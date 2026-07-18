import struct
from dataclasses import dataclass
from pathlib import Path


class ModeFileError(RuntimeError):
    pass


@dataclass(frozen=True)
class ModeSet:
    depths: list[float]
    modes: list[list[complex]]
    wavenumbers: list[complex]


def _layout(data: bytes, path: Path) -> tuple[int, int, int, int]:
    if len(data) < 100:
        raise ModeFileError(f"mode file is too short: {path}")
    (record_length_words,) = struct.unpack_from("<i", data, 0)
    if record_length_words < 25:
        raise ModeFileError(f"invalid MOD record length: {record_length_words}")
    record_bytes = 4 * record_length_words
    _, _, tabulated_count, material_count = struct.unpack_from("<iiii", data, 84)
    if tabulated_count < 1 or material_count < 1:
        raise ModeFileError("invalid MOD mode-grid dimensions")
    mode_count_offset = 5 * record_bytes
    if mode_count_offset + 4 > len(data):
        raise ModeFileError("MOD file does not contain the first mode-count record")
    (mode_count,) = struct.unpack_from("<i", data, mode_count_offset)
    if mode_count < 0:
        raise ModeFileError(f"invalid MOD mode count: {mode_count}")
    return record_bytes, tabulated_count, material_count, mode_count


def read_first_mode_set(path: Path) -> ModeSet:
    path = Path(path)
    data = path.read_bytes()
    record_bytes, tabulated_count, material_count, mode_count = _layout(data, path)

    depth_start = 4 * record_bytes
    depth_end = depth_start + 4 * tabulated_count
    if depth_end > len(data):
        raise ModeFileError("MOD file is truncated in the depth record")
    depths = list(struct.unpack_from(f"<{tabulated_count}f", data, depth_start))

    modes = []
    mode_values_bytes = 8 * material_count
    if mode_values_bytes > record_bytes:
        raise ModeFileError("MOD mode vector does not fit its direct-access record")
    for mode_index in range(mode_count):
        start = (7 + mode_index) * record_bytes
        end = start + mode_values_bytes
        if end > len(data):
            raise ModeFileError("MOD file is truncated in the mode records")
        values = []
        for value_index in range(material_count):
            real, imaginary = struct.unpack_from("<ff", data, start + 8 * value_index)
            values.append(complex(real, imaginary))
        modes.append(values)

    return ModeSet(depths, modes, read_first_wavenumber_set(path))


def read_first_wavenumber_set(path: Path) -> list[complex]:
    path = Path(path)
    data = path.read_bytes()
    record_bytes, _, _, mode_count = _layout(data, path)
    if mode_count == 0:
        return []

    first_wavenumber_record = mode_count + 8
    start = (first_wavenumber_record - 1) * record_bytes
    required_records = 1 + (2 * mode_count - 1) // (record_bytes // 4)
    end = start + required_records * record_bytes
    if end > len(data):
        raise ModeFileError("MOD file is truncated in the wavenumber records")

    result = []
    cursor = start
    while len(result) < mode_count:
        real, imaginary = struct.unpack_from("<ff", data, cursor)
        result.append(complex(real, imaginary))
        cursor += 8
    return result


def read_wavenumber_sets(path: Path) -> list[list[complex]]:
    """Read every one-frequency profile's wavenumbers from a KRAKENC MOD file."""
    path = Path(path)
    data = path.read_bytes()
    if len(data) < 100:
        raise ModeFileError(f"mode file is too short: {path}")
    (record_words,) = struct.unpack_from("<i", data, 0)
    if record_words < 25:
        raise ModeFileError(f"invalid MOD record length: {record_words}")
    record_bytes = 4 * record_words
    complex_per_record = record_words // 2
    profile_record = 0
    result = []
    while profile_record * record_bytes < len(data):
        header = profile_record * record_bytes
        if header + 100 > len(data):
            raise ModeFileError("MOD file is truncated in a profile header")
        if struct.unpack_from("<i", data, header)[0] != record_words:
            raise ModeFileError("MOD profile record length mismatch")
        mode_offset = (profile_record + 5) * record_bytes
        if mode_offset + 4 > len(data):
            raise ModeFileError("MOD file is truncated before a mode count")
        (mode_count,) = struct.unpack_from("<i", data, mode_offset)
        if mode_count < 1:
            raise ModeFileError("MOD profile must contain at least one mode")
        wave_start = profile_record + 7 + mode_count
        values = []
        for mode in range(mode_count):
            wave_record = wave_start + mode // complex_per_record
            wave_offset = (wave_record * record_bytes +
                           8 * (mode % complex_per_record))
            if wave_offset + 8 > len(data):
                raise ModeFileError("MOD file is truncated in profile wavenumbers")
            real, imaginary = struct.unpack_from("<ff", data, wave_offset)
            values.append(complex(real, imaginary))
        result.append(values)
        wave_records = (mode_count + complex_per_record - 1) // complex_per_record
        profile_record += 7 + mode_count + wave_records
    return result
