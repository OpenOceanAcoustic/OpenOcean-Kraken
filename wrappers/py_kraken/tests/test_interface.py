from pathlib import Path
import json
import sys
import unittest

import numpy as np


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE_ROOT))

import py_kraken
from py_kraken import BiologicalAttenuationLayerModel, OpenOceanKraken_interface


class HighLevelInterfaceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = (Path(__file__).resolve().parents[3] / "../test/MunkK.env").resolve()

    def test_instances_are_independent_and_singleton_is_optional(self):
        first = OpenOceanKraken_interface(thread_num=1)
        second = OpenOceanKraken_interface(thread_num=1)
        self.assertIsNot(first, second)
        self.assertIs(
            OpenOceanKraken_interface.get_instance(1),
            OpenOceanKraken_interface.get_instance(1),
        )

    def test_file_and_vector_range_validation(self):
        interface = OpenOceanKraken_interface(thread_num=1)
        with self.assertRaises(FileNotFoundError):
            interface.ook_load_env("missing.env")
        with self.assertRaises(ValueError):
            interface.ook_set_receiver_range(start=0.0, end=1000.0, count=None)
        with self.assertRaises(ValueError):
            interface.ook_set_receiver_range(values=[0.0], start=0.0, end=1.0, count=2)

    def test_high_level_biological_attenuation_preserves_order(self):
        interface = OpenOceanKraken_interface(thread_num=1)
        interface.ook_set_biological_attenuation([
            {"Z1": 10.0, "Z2": 30.0, "f0": 1000.0, "Q": 5.0, "a0": 0.04},
            {"Z1": 40.0, "Z2": 60.0, "f0": 1200.0, "Q": 4.0, "a0": 0.02},
        ])
        payload = json.loads(interface._interface.to_json_string())
        self.assertEqual(
            [item["f0"] for item in payload["AttenUnit"]["BiologicalLayers"]],
            [1000.0, 1200.0])

        interface.ook_set_biological_attenuation([])
        payload = json.loads(interface._interface.to_json_string())
        self.assertEqual(payload["AttenUnit"]["BiologicalLayers"], [])

        with self.assertRaises(ValueError):
            interface.ook_set_biological_attenuation([
                {"Z1": 30.0, "Z2": 10.0, "f0": 1000.0, "Q": 5.0, "a0": 0.04}
            ])

    def test_biological_layer_model_is_exported_from_package(self):
        self.assertIn("BiologicalAttenuationLayerModel", py_kraken.__all__)
        self.assertIs(
            py_kraken.BiologicalAttenuationLayerModel,
            BiologicalAttenuationLayerModel,
        )

    def test_biological_layer_model_rejects_extra_fields(self):
        with self.assertRaises(ValueError):
            BiologicalAttenuationLayerModel.model_validate({
                "Z1": 10.0,
                "Z2": 30.0,
                "f0": 1000.0,
                "Q": 5.0,
                "a0": 0.04,
                "unexpected": 1.0,
            })

    def test_biological_layer_model_rejects_nonfinite_fields(self):
        for field, value in (
            ("Z1", float("nan")),
            ("Z2", float("inf")),
            ("f0", float("nan")),
            ("Q", float("inf")),
            ("a0", float("nan")),
        ):
            with self.subTest(field=field):
                layer = {
                    "Z1": 10.0,
                    "Z2": 30.0,
                    "f0": 1000.0,
                    "Q": 5.0,
                    "a0": 0.04,
                }
                layer[field] = value
                with self.assertRaises(ValueError):
                    BiologicalAttenuationLayerModel.model_validate(layer)

    def test_biological_layer_model_rejects_invalid_attenuation_values(self):
        for field, value in (
            ("f0", 0.0),
            ("Q", 0.0),
            ("a0", -0.01),
        ):
            with self.subTest(field=field):
                layer = {
                    "Z1": 10.0,
                    "Z2": 30.0,
                    "f0": 1000.0,
                    "Q": 5.0,
                    "a0": 0.04,
                }
                layer[field] = value
                with self.assertRaises(ValueError):
                    BiologicalAttenuationLayerModel.model_validate(layer)

    def test_high_level_biological_attenuation_accepts_model_and_defaults_unit(self):
        interface = OpenOceanKraken_interface(thread_num=1)
        layer = BiologicalAttenuationLayerModel(
            Z1=10.0,
            Z2=30.0,
            f0=1000.0,
            Q=5.0,
            a0=0.04,
        )
        interface.ook_set_biological_attenuation([layer])
        payload = json.loads(interface._interface.to_json_string())
        self.assertEqual(payload["AttenUnit"]["AttenuationUnit"], "dB/lambda")
        self.assertEqual(
            payload["AttenUnit"]["BiologicalLayers"],
            [{"Z1": 10.0, "Z2": 30.0, "f0": 1000.0, "Q": 5.0, "a0": 0.04}],
        )

    def test_high_level_biological_attenuation_requires_list(self):
        interface = OpenOceanKraken_interface(thread_num=1)
        with self.assertRaises(ValueError):
            interface.ook_set_biological_attenuation(({
                "Z1": 10.0,
                "Z2": 30.0,
                "f0": 1000.0,
                "Q": 5.0,
                "a0": 0.04,
            },))

    def test_high_level_biological_attenuation_rejects_more_than_200_layers(self):
        interface = OpenOceanKraken_interface(thread_num=1)
        layers = [
            {"Z1": 10.0, "Z2": 30.0, "f0": 1000.0, "Q": 5.0, "a0": 0.04}
            for _ in range(201)
        ]
        with self.assertRaisesRegex(
            ValueError,
            "layers must contain at most 200 entries",
        ):
            interface.ook_set_biological_attenuation(layers)

    def test_run_returns_owned_pressure_velocity_and_modes(self):
        interface = OpenOceanKraken_interface(thread_num=1)
        interface.ook_load_env(self.fixture)
        interface.ook_set_velocity_enable(True)
        interface.ook_run()

        pressure = interface.ook_get_pressure()
        vertical = interface.ook_get_vertical_velocity()
        horizontal = interface.ook_get_horizontal_velocity()
        modes = interface.ook_get_modes()
        self.assertEqual(pressure.values.ndim, 3)
        self.assertEqual(vertical.values.shape, pressure.values.shape)
        self.assertEqual(horizontal.values.shape, pressure.values.shape)
        self.assertGreater(len(modes.profiles), 0)

        saved = pressure.values.copy()
        config = interface.ook_get_config()
        self.assertEqual(Path(config.input), self.fixture)
        interface.ook_clear()
        np.testing.assert_array_equal(pressure.values, saved)


if __name__ == "__main__":
    unittest.main()
