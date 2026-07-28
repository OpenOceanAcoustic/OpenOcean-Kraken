"""Validated Python data containers used by the OOK wrappers."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


class ConfigModel(BaseModel):
    """Common high-level runner options.

    Unknown keys are retained so applications can keep their own metadata beside
    the OOK options without needing a second configuration object.
    """

    model_config = ConfigDict(extra="allow")

    input: str | Path | None = None
    output: str | Path | None = None
    threads: int = Field(default=1, ge=1)
    velocity: bool = False
    mod: bool = False
    mod_only: bool = False
    export_json: bool = False


class BiologicalAttenuationLayerModel(BaseModel):
    model_config = ConfigDict(extra="forbid")

    Z1: float
    Z2: float
    f0: float
    Q: float
    a0: float

    @model_validator(mode="after")
    def validate_layer(self) -> "BiologicalAttenuationLayerModel":
        for name in ("Z1", "Z2", "f0", "Q", "a0"):
            if not math.isfinite(getattr(self, name)):
                raise ValueError(f"{name} must be finite")
        if self.Z1 > self.Z2:
            raise ValueError("Z1 must be less than or equal to Z2")
        if self.f0 <= 0.0:
            raise ValueError("f0 must be positive")
        if self.Q <= 0.0:
            raise ValueError("Q must be positive")
        if self.a0 < 0.0:
            raise ValueError("a0 must be non-negative")
        return self


class _ArrayModel(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    @staticmethod
    def _real_array(value: Any) -> np.ndarray:
        result = np.asarray(value)
        if result.ndim != 1:
            raise ValueError("coordinate arrays must be one-dimensional")
        return np.array(result, copy=True)


class FieldData(_ArrayModel):
    """Pressure or velocity field in ``[source, range, depth]`` order."""

    title: str = ""
    plot_type: str = ""
    frequency: float
    attenuation: float = 0.0
    source_depths: np.ndarray
    receiver_ranges: np.ndarray
    receiver_depths: np.ndarray
    values: np.ndarray
    path: Path | None = None

    @field_validator("source_depths", "receiver_ranges", "receiver_depths", mode="before")
    @classmethod
    def validate_coordinates(cls, value: Any) -> np.ndarray:
        return cls._real_array(value)

    @field_validator("values", mode="before")
    @classmethod
    def validate_values(cls, value: Any) -> np.ndarray:
        result = np.asarray(value, dtype=np.complex64)
        if result.ndim != 3:
            raise ValueError("field values must have [source, range, depth] dimensions")
        return np.array(result, dtype=np.complex64, copy=True, order="C")

    @model_validator(mode="after")
    def validate_shape(self) -> "FieldData":
        source_count, range_count, depth_count = self.values.shape
        if source_count != self.source_depths.size:
            raise ValueError("field source dimension does not match source_depths")
        if range_count != self.receiver_ranges.size:
            raise ValueError("field range dimension does not match receiver_ranges")
        if self.plot_type.lower().startswith("irregular"):
            if depth_count != 1 or self.receiver_depths.size != range_count:
                raise ValueError("irregular fields require one depth value per range")
        elif depth_count != self.receiver_depths.size:
            raise ValueError("field depth dimension does not match receiver_depths")
        return self


class ModeProfileData(_ArrayModel):
    """Modes for one range-dependent environment profile."""

    profile_index: int = Field(ge=0)
    profile_range: float | None = None
    depth: np.ndarray
    wavenumbers: np.ndarray
    group_velocity: np.ndarray = Field(default_factory=lambda: np.empty(0))
    mode_shapes: np.ndarray
    media_mesh_counts: list[int] = Field(default_factory=list)
    media_materials: list[str] = Field(default_factory=list)
    halfspaces: dict[str, Any] = Field(default_factory=dict)

    @field_validator("depth", "group_velocity", mode="before")
    @classmethod
    def validate_real_arrays(cls, value: Any) -> np.ndarray:
        return cls._real_array(value)

    @field_validator("wavenumbers", mode="before")
    @classmethod
    def validate_wavenumbers(cls, value: Any) -> np.ndarray:
        result = np.asarray(value, dtype=np.complex64)
        if result.ndim != 1:
            raise ValueError("wavenumbers must be one-dimensional")
        return np.array(result, dtype=np.complex64, copy=True)

    @field_validator("mode_shapes", mode="before")
    @classmethod
    def validate_mode_shapes(cls, value: Any) -> np.ndarray:
        result = np.asarray(value, dtype=np.complex64)
        if result.ndim != 2:
            raise ValueError("mode_shapes must have [mode, depth] dimensions")
        return np.array(result, dtype=np.complex64, copy=True, order="C")

    @model_validator(mode="after")
    def validate_shape(self) -> "ModeProfileData":
        if self.mode_shapes.shape != (self.wavenumbers.size, self.depth.size):
            raise ValueError("mode_shapes dimensions must match wavenumbers and depth")
        if self.group_velocity.size not in (0, self.wavenumbers.size):
            raise ValueError("group_velocity must be empty or contain one value per mode")
        return self


class ModeData(_ArrayModel):
    """All profiles stored in an OOK MOD file or returned by the native module."""

    title: str = ""
    frequency: float
    profiles: list[ModeProfileData]
    path: Path | None = None

    @model_validator(mode="after")
    def require_profiles(self) -> "ModeData":
        if not self.profiles:
            raise ValueError("mode data must contain at least one profile")
        return self


__all__ = [
    "BiologicalAttenuationLayerModel",
    "ConfigModel",
    "FieldData",
    "ModeData",
    "ModeProfileData",
]
