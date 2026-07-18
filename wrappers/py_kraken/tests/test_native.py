from pathlib import Path
import sys
import unittest

import numpy as np


PACKAGE_DIR = Path(__file__).resolve().parents[1] / "py_kraken"
sys.path.insert(0, str(PACKAGE_DIR))

import OpenOceanKraken as native


class NativeBindingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = Path(__file__).resolve().parents[3] / "../test/MunkK.env"
        cls.fixture = cls.fixture.resolve()

    def test_required_symbols_are_exported(self):
        for name in (
            "Interface",
            "ThreadPool",
            "Run_Mode",
            "SSP_Mode",
            "Media_Mode",
            "BC_Mode",
            "Grid_Mode",
            "CoherenceType",
            "ModeType",
            "SSPLayer",
            "Range_Independent_Area",
        ):
            self.assertTrue(hasattr(native, name), name)

    def test_pressure_is_owned_complex64_source_range_depth_array(self):
        interface = native.Interface()
        pool = native.ThreadPool(1)
        interface.setThreadPool(pool)
        interface.setNumThreads(1)
        self.assertTrue(interface.from_env(str(self.fixture)))
        interface.set_Velocity_enable(True)
        interface.run()

        pressure = interface.get_pressure()
        self.assertEqual(pressure.dtype, np.complex64)
        self.assertEqual(pressure.ndim, 3)
        snapshot = interface.get_pressure_snapshot()
        self.assertEqual(
            pressure.shape,
            (snapshot.source_count, snapshot.range_count, snapshot.depth_count),
        )
        np.testing.assert_array_equal(pressure, snapshot.values)

        saved = pressure.copy()
        interface.clearResults()
        np.testing.assert_array_equal(pressure, saved)

    def test_modes_are_owned_profile_values(self):
        interface = native.Interface()
        pool = native.ThreadPool(1)
        interface.setThreadPool(pool)
        interface.setNumThreads(1)
        self.assertTrue(interface.from_env(str(self.fixture)))
        interface.runEigen()
        profiles = interface.get_modes()
        self.assertGreater(len(profiles), 0)
        self.assertGreater(profiles[0].wavenumbers.size, 0)
        self.assertEqual(
            profiles[0].mode_shapes.shape[0],
            profiles[0].wavenumbers.size,
        )


if __name__ == "__main__":
    unittest.main()
