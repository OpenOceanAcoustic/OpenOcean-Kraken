"""Python interface for OpenOceanKraken."""

import sys

from .ook_data_model import (
    BiologicalAttenuationLayerModel,
    ConfigModel,
    FieldData,
    ModeData,
    ModeProfileData,
)
from .ook_read import read_json, read_mod, read_shd

try:
    if "OpenOceanKraken" in sys.modules:
        native = sys.modules["OpenOceanKraken"]
    else:
        from . import OpenOceanKraken as native

        sys.modules["OpenOceanKraken"] = native
except ImportError:
    native = None

if native is not None:
    from .ook_interface import OpenOceanKraken_interface
else:
    OpenOceanKraken_interface = None

from .ook_plot import plot_field, plot_modes

__all__ = [
    "BiologicalAttenuationLayerModel",
    "ConfigModel",
    "FieldData",
    "ModeData",
    "ModeProfileData",
    "OpenOceanKraken_interface",
    "native",
    "plot_field",
    "plot_modes",
    "read_json",
    "read_mod",
    "read_shd",
]
