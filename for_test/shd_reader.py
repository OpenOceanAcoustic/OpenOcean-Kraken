import struct
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ShadeFile:
    title: str
    frequency: float
    source_depths: list[float]
    receiver_depths: list[float]
    ranges_metres: list[float]
    pressure: list[complex]


def read_shade_file(path: Path) -> ShadeFile:
    path = Path(path)
    data = path.read_bytes()
    if len(data) < 4:
        raise ValueError(f"SHD file is too short: {path}")
    (record_words,) = struct.unpack_from("<i", data, 0)
    record_bytes = 4 * record_words
    if record_words < 41 or len(data) < 10 * record_bytes:
        raise ValueError("invalid SHD record layout")
    title = data[4:84].decode("ascii", errors="replace").rstrip()
    header = 2 * record_bytes
    nfreq, ntheta, nsx, nsy, nsz, nrz, nrr = struct.unpack_from("<7i", data, header)
    (frequency,) = struct.unpack_from("<d", data, 3 * record_bytes)

    def doubles(record: int, count: int) -> list[float]:
        return list(struct.unpack_from(f"<{count}d", data, record * record_bytes))

    def floats(record: int, count: int) -> list[float]:
        return list(struct.unpack_from(f"<{count}f", data, record * record_bytes))

    source_depths = floats(7, nsz)
    receiver_depths = floats(8, nrz)
    ranges = doubles(9, nrr)
    pressure = []
    record = 10
    for _ in range(nfreq * ntheta * nsx * nsy * nsz * nrz):
        start = record * record_bytes
        for index in range(nrr):
            real, imaginary = struct.unpack_from("<ff", data, start + 8 * index)
            pressure.append(complex(real, imaginary))
        record += 1
    return ShadeFile(title, frequency, source_depths, receiver_depths, ranges, pressure)
