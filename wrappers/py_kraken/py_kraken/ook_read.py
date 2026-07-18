"""Readers for OOK's fixed-record SHD and MOD binary formats."""

from __future__ import annotations

from pathlib import Path
import json
import struct
from typing import Any

import numpy as np

from .ook_data_model import ConfigModel, FieldData, ModeData, ModeProfileData


def _load(path: str | Path, kind: str) -> tuple[Path, bytes]:
    resolved = Path(path)
    try:
        raw = resolved.read_bytes()
    except OSError as error:
        raise ValueError(f"cannot read {kind} file {resolved}: {error}") from error
    if len(raw) < 4:
        raise ValueError(f"{kind} file {resolved} is too short")
    return resolved, raw


def _record_length(raw: bytes, kind: str) -> tuple[int, int]:
    words = struct.unpack_from("<i", raw, 0)[0]
    if words <= 0:
        raise ValueError(f"{kind} file has invalid record length {words}")
    record_bytes = 4 * words
    if record_bytes < 4 or record_bytes > len(raw):
        raise ValueError(f"{kind} file has impossible record length {words}")
    return words, record_bytes


def _record(raw: bytes, record_bytes: int, index: int, kind: str) -> memoryview:
    start = index * record_bytes
    end = start + record_bytes
    if index < 0 or end > len(raw):
        raise ValueError(f"{kind} file is truncated at record {index + 1}")
    return memoryview(raw)[start:end]


def _array(record: memoryview, dtype: str, count: int, kind: str) -> np.ndarray:
    item_size = np.dtype(dtype).itemsize
    if count < 0 or count * item_size > len(record):
        raise ValueError(f"{kind} record cannot contain {count} values of type {dtype}")
    return np.frombuffer(record, dtype=dtype, count=count).copy()


def _text(data: bytes | memoryview) -> str:
    return bytes(data).decode("latin-1", errors="replace").rstrip("\x00 ").strip()


def read_json(path: str | Path) -> ConfigModel:
    """Read an OOK JSON input while retaining extension fields."""

    resolved = Path(path)
    try:
        with resolved.open("r", encoding="utf-8") as stream:
            value = json.load(stream)
    except OSError as error:
        raise ValueError(f"cannot read JSON file {resolved}: {error}") from error
    except json.JSONDecodeError as error:
        raise ValueError(f"invalid JSON file {resolved}: {error}") from error
    if not isinstance(value, dict):
        raise ValueError("OOK JSON root must be an object")
    value.setdefault("input", resolved)
    return ConfigModel.model_validate(value)


def read_shd(path: str | Path) -> FieldData:
    """Read an OOK SHD file into an owned ``complex64`` field array."""

    resolved, raw = _load(path, "SHD")
    _, record_bytes = _record_length(raw, "SHD")
    if len(raw) < 10 * record_bytes:
        raise ValueError("SHD file is missing one or more header records")

    title = _text(memoryview(raw)[4:84])
    plot_type = _text(_record(raw, record_bytes, 1, "SHD")[:10]).lower()
    header = _record(raw, record_bytes, 2, "SHD")
    if len(header) < 36:
        raise ValueError("SHD dimension record is too short")
    nfreq, ntheta, nsx, nsy, nsz, nrz, nrr = struct.unpack_from("<7i", header, 0)
    frequency, attenuation = struct.unpack_from("<2f", header, 28)
    dimensions = (nfreq, ntheta, nsx, nsy, nsz, nrz, nrr)
    if any(value <= 0 for value in dimensions):
        raise ValueError(f"SHD file has invalid dimensions {dimensions}")

    frequencies = _array(_record(raw, record_bytes, 3, "SHD"), "<f8", nfreq, "SHD")
    if frequencies.size:
        frequency = float(frequencies[0])
    source_depths = _array(_record(raw, record_bytes, 7, "SHD"), "<f4", nsz, "SHD")
    all_receiver_depths = _array(
        _record(raw, record_bytes, 8, "SHD"), "<f4", nrz, "SHD"
    )
    receiver_ranges = _array(
        _record(raw, record_bytes, 9, "SHD"), "<f4", nrr, "SHD"
    )

    last_record_index = (len(raw) - 1) // record_bytes
    data_record_count = last_record_index - 9
    if data_record_count <= 0 or data_record_count % nsz:
        raise ValueError("SHD data-record count is inconsistent with its source count")
    depth_count = data_record_count // nsz
    irregular = plot_type.startswith("irregular")
    expected_depth_count = 1 if irregular else nrz
    if depth_count != expected_depth_count:
        raise ValueError(
            f"SHD contains {depth_count} depth records per source; expected {expected_depth_count}"
        )

    values = np.empty((nsz, nrr, depth_count), dtype=np.complex64)
    for source_index in range(nsz):
        for depth_index in range(depth_count):
            record_index = 10 + source_index * depth_count + depth_index
            start = record_index * record_bytes
            end = start + nrr * np.dtype("<c8").itemsize
            if end > len(raw):
                raise ValueError(f"SHD file is truncated at data record {record_index + 1}")
            values[source_index, :, depth_index] = np.frombuffer(
                raw, dtype="<c8", count=nrr, offset=start
            )

    receiver_depths = all_receiver_depths if irregular else all_receiver_depths[:depth_count]
    return FieldData(
        title=title,
        plot_type=plot_type,
        frequency=frequency,
        attenuation=attenuation,
        source_depths=source_depths,
        receiver_ranges=receiver_ranges,
        receiver_depths=receiver_depths,
        values=values,
        path=resolved,
    )


def _parse_halfspaces(record: memoryview) -> dict[str, Any]:
    if len(record) < 50:
        raise ValueError("MOD halfspace record is too short")

    def one(offset: int) -> dict[str, Any]:
        cp = complex(*struct.unpack_from("<2f", record, offset + 1))
        cs = complex(*struct.unpack_from("<2f", record, offset + 9))
        density, depth = struct.unpack_from("<2f", record, offset + 17)
        return {
            "boundary": _text(record[offset : offset + 1]),
            "compressional_speed": cp,
            "shear_speed": cs,
            "density": float(density),
            "depth": float(depth),
        }

    return {"top": one(0), "bottom": one(25)}


def read_mod(path: str | Path) -> ModeData:
    """Read every profile in an OOK MOD file."""

    resolved, raw = _load(path, "MOD")
    recl_words, record_bytes = _record_length(raw, "MOD")
    profiles: list[ModeProfileData] = []
    record_index = 0
    common_title = ""
    common_frequency: float | None = None

    while record_index * record_bytes < len(raw):
        header = _record(raw, record_bytes, record_index, "MOD")
        if len(header) < 100:
            raise ValueError("MOD profile header is too short")
        embedded_recl = struct.unpack_from("<i", header, 0)[0]
        if embedded_recl != recl_words:
            raise ValueError(
                f"MOD profile {len(profiles)} has record length {embedded_recl}, expected {recl_words}"
            )
        title = _text(header[4:84])
        nfreq, nmedia, depth_count, material_count = struct.unpack_from("<4i", header, 84)
        if nfreq <= 0 or nmedia <= 0 or depth_count <= 0 or material_count < 0:
            raise ValueError("MOD profile header contains invalid dimensions")

        media_record = _record(raw, record_bytes, record_index + 1, "MOD")
        if nmedia * 12 > len(media_record):
            raise ValueError("MOD media record is too short")
        mesh_counts: list[int] = []
        materials: list[str] = []
        for medium_index in range(nmedia):
            offset = medium_index * 12
            mesh_counts.append(struct.unpack_from("<i", media_record, offset)[0])
            materials.append(_text(media_record[offset + 4 : offset + 12]))

        frequency_record = _record(raw, record_bytes, record_index + 3, "MOD")
        frequencies = _array(frequency_record, "<f8", nfreq, "MOD")
        depth = _array(
            _record(raw, record_bytes, record_index + 4, "MOD"),
            "<f4",
            depth_count,
            "MOD",
        )
        mode_record_index = record_index + 5
        mode_count = struct.unpack_from(
            "<i", _record(raw, record_bytes, mode_record_index, "MOD"), 0
        )[0]
        if mode_count < 0:
            raise ValueError(f"MOD profile {len(profiles)} has invalid mode count {mode_count}")
        halfspaces = _parse_halfspaces(
            _record(raw, record_bytes, mode_record_index + 1, "MOD")
        )

        mode_shapes = np.empty((mode_count, depth_count), dtype=np.complex64)
        for mode_index in range(mode_count):
            mode_shapes[mode_index] = _array(
                _record(
                    raw,
                    record_bytes,
                    mode_record_index + 2 + mode_index,
                    "MOD",
                ),
                "<c8",
                depth_count,
                "MOD",
            )

        modes_per_record = max(1, recl_words // 2)
        wavenumbers = np.empty(mode_count, dtype=np.complex64)
        k_record_count = (mode_count + modes_per_record - 1) // modes_per_record
        first_mode = 0
        for k_index in range(k_record_count):
            take = min(modes_per_record, mode_count - first_mode)
            wavenumbers[first_mode : first_mode + take] = _array(
                _record(
                    raw,
                    record_bytes,
                    mode_record_index + 2 + mode_count + k_index,
                    "MOD",
                ),
                "<c8",
                take,
                "MOD",
            )
            first_mode += take

        frequency = float(frequencies[0])
        if common_frequency is not None and not np.isclose(common_frequency, frequency):
            raise ValueError("MOD profiles contain inconsistent frequencies")
        common_frequency = frequency
        common_title = common_title or title
        profiles.append(
            ModeProfileData(
                profile_index=len(profiles),
                profile_range=None,
                depth=depth,
                wavenumbers=wavenumbers,
                group_velocity=np.empty(0),
                mode_shapes=mode_shapes,
                media_mesh_counts=mesh_counts,
                media_materials=materials,
                halfspaces=halfspaces,
            )
        )
        record_index += 8 + mode_count + (
            (2 * mode_count - 1) // recl_words if mode_count else 0
        )

    if not profiles or common_frequency is None:
        raise ValueError("MOD file contains no profiles")
    return ModeData(
        title=common_title,
        frequency=common_frequency,
        profiles=profiles,
        path=resolved,
    )


__all__ = ["read_json", "read_mod", "read_shd"]
