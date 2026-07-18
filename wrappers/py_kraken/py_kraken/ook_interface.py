"""Pythonic facade over the native :mod:`OpenOceanKraken` extension."""

from __future__ import annotations

import json
import math
from pathlib import Path
from threading import Lock
from typing import Any, Callable

import numpy as np

from . import native
from .ook_data_model import ConfigModel, FieldData, ModeData, ModeProfileData


class OpenOceanKraken_interface:
    """High-level, independently constructible OOK solver interface."""

    _instances: dict[int, "OpenOceanKraken_interface"] = {}
    _instances_lock = Lock()

    def __init__(self, thread_num: int = 1):
        if isinstance(thread_num, bool) or not isinstance(thread_num, int) or thread_num <= 0:
            raise ValueError("thread_num must be a positive integer")
        self._thread_num = thread_num
        self._pool = native.ThreadPool(thread_num)
        self._interface = native.Interface(self._pool)
        self._interface.setNumThreads(thread_num)
        self._input_path: Path | None = None

    @classmethod
    def get_instance(cls, thread_num: int = 1) -> "OpenOceanKraken_interface":
        """Return an optional per-thread-count singleton."""

        with cls._instances_lock:
            if thread_num not in cls._instances:
                cls._instances[thread_num] = cls(thread_num)
            return cls._instances[thread_num]

    @staticmethod
    def _input_file(path: str | Path, suffix: str) -> Path:
        resolved = Path(path).expanduser().resolve()
        if not resolved.is_file():
            raise FileNotFoundError(resolved)
        if resolved.suffix.lower() != suffix:
            raise ValueError(f"expected a {suffix} input file")
        return resolved

    @staticmethod
    def _finite(value: float, name: str, *, positive: bool = False) -> float:
        result = float(value)
        if not math.isfinite(result) or (positive and result <= 0.0):
            qualifier = "positive and finite" if positive else "finite"
            raise ValueError(f"{name} must be {qualifier}")
        return result

    @staticmethod
    def _vector(value: Any, name: str) -> np.ndarray:
        result = np.asarray(value, dtype=np.float64)
        if result.ndim != 1 or result.size == 0:
            raise ValueError(f"{name} must be a non-empty one-dimensional array")
        if not np.isfinite(result).all():
            raise ValueError(f"{name} must contain only finite values")
        return np.ascontiguousarray(result)

    def _set_vector_or_range(
        self,
        setter: Callable[..., None],
        name: str,
        values: Any = None,
        *,
        start: float | None = None,
        end: float | None = None,
        count: int | None = None,
    ) -> None:
        range_arguments = (start, end, count)
        if values is not None:
            if any(argument is not None for argument in range_arguments):
                raise ValueError(f"{name}: use values or start/end/count, not both")
            setter(self._vector(values, name))
            return
        if any(argument is None for argument in range_arguments):
            raise ValueError(f"{name}: start, end, and count are all required")
        if isinstance(count, bool) or not isinstance(count, int) or count <= 0:
            raise ValueError(f"{name}: count must be a positive integer")
        range_start = self._finite(start, f"{name} start")
        range_end = self._finite(end, f"{name} end")
        setter(range_start, range_end, count)

    @staticmethod
    def _output_root(path: str | Path, suffix: str) -> str:
        result = Path(path).expanduser()
        if result.suffix.lower() == suffix:
            result = result.with_suffix("")
        result.parent.mkdir(parents=True, exist_ok=True)
        return str(result)

    def ook_load_env(self, path: str | Path) -> None:
        resolved = self._input_file(path, ".env")
        if not self._interface.from_env(str(resolved)):
            raise RuntimeError(f"OOK failed to load ENV file {resolved}")
        self._input_path = resolved

    def ook_load_json(self, path: str | Path) -> None:
        resolved = self._input_file(path, ".json")
        if not self._interface.from_json(str(resolved)):
            raise RuntimeError(f"OOK failed to load JSON file {resolved}")
        self._input_path = resolved

    def ook_export_json(self, path: str | Path) -> None:
        output = Path(path).expanduser()
        if output.suffix.lower() != ".json":
            output = output.with_suffix(".json")
        output.parent.mkdir(parents=True, exist_ok=True)
        if not self._interface.to_json(str(output)):
            raise RuntimeError(f"OOK failed to export JSON file {output}")

    def ook_get_config(self) -> ConfigModel:
        value = json.loads(self._interface.to_json_string())
        value["input"] = self._input_path
        value["threads"] = self._thread_num
        return ConfigModel.model_validate(value)

    def ook_set_frequency(self, frequency: float) -> None:
        self._interface.set_Freq(self._finite(frequency, "frequency", positive=True))

    def ook_set_frequency_vector(self, frequencies: np.ndarray) -> None:
        values = self._vector(frequencies, "frequencies")
        if np.any(values <= 0.0):
            raise ValueError("frequencies must be positive")
        self._interface.set_freqvec(values)

    def ook_set_title(self, title: str) -> None:
        if not isinstance(title, str) or not title.strip():
            raise ValueError("title must be a non-empty string")
        self._interface.set_Title(title)

    def ook_set_source_depth(self, values=None, *, start=None, end=None, count=None) -> None:
        self._set_vector_or_range(
            self._interface.set_Sz, "source depth", values, start=start, end=end, count=count
        )

    def ook_set_receiver_depth(self, values=None, *, start=None, end=None, count=None) -> None:
        self._set_vector_or_range(
            self._interface.set_Rz, "receiver depth", values, start=start, end=end, count=count
        )

    def ook_set_receiver_range(self, values=None, *, start=None, end=None, count=None) -> None:
        self._set_vector_or_range(
            self._interface.set_Rr, "receiver range", values, start=start, end=end, count=count
        )

    def ook_set_profile_range(self, values=None, *, start=None, end=None, count=None) -> None:
        self._set_vector_or_range(
            self._interface.set_RProf, "profile range", values, start=start, end=end, count=count
        )

    def ook_set_phase_speed(self, c_low: float, c_high: float) -> None:
        low = self._finite(c_low, "c_low", positive=True)
        high = self._finite(c_high, "c_high", positive=True)
        if low >= high:
            raise ValueError("c_low must be smaller than c_high")
        self._interface.set_cPhase(low, high)

    def ook_set_mode_limit(self, limit: int) -> None:
        if isinstance(limit, bool) or not isinstance(limit, int) or limit <= 0:
            raise ValueError("mode limit must be a positive integer")
        self._interface.set_MLimit(limit)

    def ook_set_ssp(self, areas: list[Any]) -> None:
        if not isinstance(areas, list) or not areas:
            raise ValueError("areas must be a non-empty list")
        self._interface.set_SSP(areas)

    def ook_set_attenuation(self, mode: Any) -> None:
        self._interface.set_AttenUnit(mode)

    def ook_set_grid_type(self, mode: Any) -> None:
        self._interface.set_GridType(mode)

    def ook_set_source_type(self, mode: Any) -> None:
        self._interface.set_SourceType(mode)

    def ook_set_run_mode(self, mode: Any) -> None:
        self._interface.set_RunMode(mode)

    def ook_set_coherence_type(self, mode: Any) -> None:
        self._interface.set_CoherenceType(mode)

    def ook_set_mode_type(self, mode: Any) -> None:
        self._interface.set_ModeType(mode)

    def ook_set_rmax(self, rmax: float) -> None:
        self._interface.set_Rmax(self._finite(rmax, "rmax", positive=True))

    @staticmethod
    def _reflection_coefficients(coefficients: np.ndarray) -> list[Any]:
        values = np.asarray(coefficients, dtype=np.float64)
        if values.ndim != 2 or values.shape[1] != 3 or values.shape[0] == 0:
            raise ValueError("reflection coefficients must have [theta, magnitude, phase] columns")
        if not np.isfinite(values).all():
            raise ValueError("reflection coefficients must be finite")
        return [native.ReflectionCoef(*row) for row in values]

    def ook_set_reflection_top(self, coefficients: np.ndarray) -> None:
        self._interface.set_ReflCoef_Top(self._reflection_coefficients(coefficients))

    def ook_set_reflection_bottom(self, coefficients: np.ndarray) -> None:
        self._interface.set_ReflCoef_Bottom(self._reflection_coefficients(coefficients))

    def ook_set_sbp(self, pattern: np.ndarray, angles: np.ndarray) -> None:
        pattern_values = self._vector(pattern, "pattern")
        angle_values = self._vector(angles, "angles")
        if pattern_values.size != angle_values.size:
            raise ValueError("pattern and angles must have the same length")
        self._interface.set_SBP(pattern_values, angle_values)

    def ook_set_velocity_enable(self, enabled: bool) -> None:
        if not isinstance(enabled, (bool, np.bool_)):
            raise ValueError("enabled must be boolean")
        self._interface.set_Velocity_enable(bool(enabled))

    def ook_run(self) -> None:
        self._interface.run()

    def ook_run_eigen(self) -> None:
        self._interface.runEigen()

    def ook_run_field(self) -> None:
        self._interface.runField()

    def ook_clear(self) -> None:
        self._interface.clearResults()

    def ook_free(self) -> None:
        self._interface.free()

    @staticmethod
    def _field(snapshot: Any) -> FieldData:
        irregular = snapshot.grid_type == native.Grid_Mode.MODE_I_Irregular
        return FieldData(
            title=snapshot.title,
            plot_type="irregular" if irregular else "rectilin",
            frequency=snapshot.frequency,
            source_depths=snapshot.source_depths,
            receiver_ranges=snapshot.receiver_ranges,
            receiver_depths=snapshot.receiver_depths,
            values=snapshot.values,
        )

    def ook_get_pressure(self) -> FieldData:
        return self._field(self._interface.get_pressure_snapshot())

    def ook_get_vertical_velocity(self) -> FieldData:
        return self._field(self._interface.get_vertical_velocity_snapshot())

    def ook_get_horizontal_velocity(self) -> FieldData:
        return self._field(self._interface.get_horizontal_velocity_snapshot())

    def ook_get_modes(self) -> ModeData:
        snapshots = self._interface.get_modes()
        profiles = [
            ModeProfileData(
                profile_index=index,
                profile_range=snapshot.profile_range,
                depth=snapshot.depth,
                wavenumbers=snapshot.wavenumbers,
                group_velocity=snapshot.group_velocity,
                mode_shapes=snapshot.mode_shapes,
            )
            for index, snapshot in enumerate(snapshots)
        ]
        config = json.loads(self._interface.to_json_string())
        return ModeData(
            title=str(config.get("Title", "")),
            frequency=float(config.get("freqinfo", {}).get("freq", 0.0)),
            profiles=profiles,
        )

    def ook_export_shd(self, path: str | Path, data_type: str = "pressure") -> None:
        types = {"pressure": 1, "vertical_velocity": 2, "horizontal_velocity": 3}
        if data_type not in types:
            raise ValueError(f"data_type must be one of {', '.join(types)}")
        self._interface.export_shd(self._output_root(path, ".shd"), types[data_type])

    def ook_export_mod(self, path: str | Path) -> None:
        self._interface.export_mod(self._output_root(path, ".mod"))

    def ook_export_result(self, path: str | Path) -> None:
        self._interface.export_result(str(Path(path).expanduser()))


__all__ = ["OpenOceanKraken_interface", "native"]
