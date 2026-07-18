"""Matplotlib visualizations for OOK fields and modes."""

from __future__ import annotations

from typing import Iterable

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import numpy as np

from .ook_data_model import FieldData, ModeData


def plot_field(
    field: FieldData,
    source_index: int = 0,
    *,
    component: str = "pressure",
    axes: Axes | None = None,
    color_limits: tuple[float, float] | None = None,
) -> Axes:
    """Plot transmission loss for one source depth and return its axes."""

    if source_index < 0 or source_index >= field.source_depths.size:
        raise IndexError("source_index is outside the field")
    if axes is None:
        _, axes = plt.subplots()
    magnitude = np.abs(field.values[source_index])
    tiny = np.finfo(np.float32).tiny
    transmission_loss = -20.0 * np.log10(np.maximum(magnitude, tiny))
    label = {
        "pressure": "Transmission loss (dB)",
        "vertical_velocity": "Vertical-velocity level (dB)",
        "horizontal_velocity": "Horizontal-velocity level (dB)",
    }.get(component)
    if label is None:
        raise ValueError("unknown field component")

    if field.plot_type.lower().startswith("irregular"):
        artist = axes.scatter(
            field.receiver_ranges / 1000.0,
            field.receiver_depths,
            c=transmission_loss[:, 0],
            s=8,
        )
    else:
        artist = axes.pcolormesh(
            field.receiver_ranges / 1000.0,
            field.receiver_depths,
            transmission_loss.T,
            shading="auto",
        )
    if color_limits is not None:
        artist.set_clim(*color_limits)
    axes.figure.colorbar(artist, ax=axes, label=label)
    axes.set_xlabel("Range (km)")
    axes.set_ylabel("Depth (m)")
    axes.set_title(
        f"{field.title} — {field.frequency:g} Hz, source {field.source_depths[source_index]:g} m"
    )
    axes.invert_yaxis()
    return axes


def plot_modes(
    modes: ModeData,
    profile_index: int = 0,
    *,
    mode_indices: Iterable[int] | None = None,
) -> Figure:
    """Plot modal wavenumbers and normalized shapes for one profile."""

    if profile_index < 0 or profile_index >= len(modes.profiles):
        raise IndexError("profile_index is outside the mode data")
    profile = modes.profiles[profile_index]
    if mode_indices is None:
        selected = list(range(min(6, profile.wavenumbers.size)))
    else:
        selected = [int(index) for index in mode_indices]
    if not selected or any(index < 0 or index >= profile.wavenumbers.size for index in selected):
        raise IndexError("mode_indices contain an invalid mode")

    figure, (wave_axes, shape_axes) = plt.subplots(1, 2, figsize=(11, 5))
    mode_numbers = np.arange(1, profile.wavenumbers.size + 1)
    wave_axes.plot(mode_numbers, profile.wavenumbers.real, label="Re(k)")
    wave_axes.plot(mode_numbers, profile.wavenumbers.imag, label="Im(k)")
    wave_axes.set_xlabel("Mode number")
    wave_axes.set_ylabel("Wavenumber (1/m)")
    wave_axes.legend()
    wave_axes.grid(True, alpha=0.25)

    for mode_index in selected:
        shape = profile.mode_shapes[mode_index]
        scale = np.max(np.abs(shape))
        normalized = shape.real / scale if scale > 0.0 else shape.real
        shape_axes.plot(normalized, profile.depth, label=f"Mode {mode_index + 1}")
    shape_axes.set_xlabel("Normalized real mode shape")
    shape_axes.set_ylabel("Depth (m)")
    shape_axes.invert_yaxis()
    shape_axes.legend()
    shape_axes.grid(True, alpha=0.25)
    figure.suptitle(f"{modes.title} — {modes.frequency:g} Hz, profile {profile_index}")
    figure.tight_layout()
    return figure


__all__ = ["plot_field", "plot_modes"]
