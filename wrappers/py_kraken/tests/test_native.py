from pathlib import Path
import json
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

    def test_biological_attenuation_binding_and_value_error(self):
        layer = native.BiologicalAttenuationLayer()
        layer.Z1 = 10.0
        layer.Z2 = 30.0
        layer.f0 = 1000.0
        layer.Q = 5.0
        layer.a0 = 0.04

        mode = native.Atten_Mode()
        mode.absModel = native.OceanAbsorptionModel.Biological
        mode.biologicalLayers = [layer]

        interface = native.Interface()
        interface.set_AttenUnit(mode)
        payload = json.loads(interface.to_json_string())
        self.assertEqual(
            payload["AttenUnit"]["OceanAbsorptionModel"], "Biological")
        self.assertEqual(
            payload["AttenUnit"]["BiologicalLayers"][0]["f0"], 1000.0)

        for description, changes in (
            ("non-finite Z1", {"Z1": float("nan")}),
            ("Z1 greater than Z2", {"Z1": 30.0, "Z2": 10.0}),
            ("non-positive f0", {"f0": 0.0}),
            ("non-positive Q", {"Q": 0.0}),
            ("negative a0", {"a0": -0.01}),
        ):
            with self.subTest(description=description):
                invalid = native.BiologicalAttenuationLayer()
                invalid.Z1 = 10.0
                invalid.Z2 = 30.0
                invalid.f0 = 1000.0
                invalid.Q = 5.0
                invalid.a0 = 0.04
                for name, value in changes.items():
                    setattr(invalid, name, value)
                mode.biologicalLayers = [invalid]
                with self.assertRaises(ValueError):
                    interface.set_AttenUnit(mode)

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
