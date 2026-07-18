from pathlib import Path
import sys
import unittest

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE_ROOT))

from py_kraken import OpenOceanKraken_interface, plot_field, plot_modes


class PlotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        fixture = (Path(__file__).resolve().parents[3] / "../test/MunkK.env").resolve()
        interface = OpenOceanKraken_interface(thread_num=1)
        interface.ook_load_env(fixture)
        interface.ook_run()
        cls.field = interface.ook_get_pressure()
        cls.modes = interface.ook_get_modes()

    def tearDown(self):
        plt.close("all")

    def test_field_plot_uses_range_km_and_inverted_depth(self):
        axes = plot_field(self.field)
        self.assertEqual(axes.get_xlabel(), "Range (km)")
        self.assertEqual(axes.get_ylabel(), "Depth (m)")
        self.assertTrue(axes.yaxis_inverted())

    def test_mode_plot_has_wavenumber_and_shape_axes(self):
        figure = plot_modes(self.modes, profile_index=0)
        self.assertGreaterEqual(len(figure.axes), 2)


if __name__ == "__main__":
    unittest.main()
