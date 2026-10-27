import unittest
import numpy as np
from phonon_absorption import (
    Channel,
    VibronicChannel,
    model_ii,
    spectrum,
    spectrum_independent,
    prepare_spectrum,
    model_ii_components,
)
from phonon_absorption.thermal import distribution
from phonon_absorption.parameters import ANCHORS
from phonon_absorption import _legacy_model_ii_independent as original
from phonon_absorption.fitting import FitSpec, fit_fixed_energies, predict_fit
from phonon_absorption.demo import synthetic_case


class OptimizedPipeline(unittest.TestCase):
    def test_adaptive_tail_and_strong_coupling(self):
        for s, temperature in [(0, 0), (2.3, 0), (7, 0.0232), (144, 0.08)]:
            adaptive = distribution(s, 0.032, temperature, cutoff=None)
            self.assertLessEqual(adaptive.omitted_mass, 1e-12)
            self.assertAlmostEqual(adaptive.weights.sum(), 1, places=10)
        manual = distribution(144, 0.032, 0.08, cutoff=20)
        self.assertGreater(manual.omitted_mass, 0.9)
        with self.assertRaises(ValueError):
            distribution(144, 0.032, 0.08, cutoff=None, max_terms=10)

    def test_large_low_temperature_support_regression(self):
        # This request triggered Boost's high-order gamma overflow in CDFs
        # on the validated SciPy version. Its physical distribution is finite.
        values = distribution(0.01, 0.032, 0.001, cutoff=1200)
        self.assertTrue(np.isfinite(values.weights).all())
        self.assertAlmostEqual(float(values.weights.sum()), 1.0, places=12)
        self.assertAlmostEqual(float(values.weights[1200]), np.exp(-0.01), places=12)

    def test_independent_frequency_and_linewidth_source_equivalence(self):
        x = np.linspace(0.85, 1.45, 157)
        modes = [
            VibronicChannel(Channel(1.281, 0.85, 0.45**2, 0.6), 0.033, 0.012),
            VibronicChannel(Channel(0.897, 1.34, 0.51**2, 0.48), 0.0295, 0.0201),
        ]
        # Match the source's [-20,19] truncation, separately from cutoff convergence.
        expanded = spectrum_independent(x, modes, "II", 0.0232, cutoff=20)
        correction = np.zeros_like(x)
        for item in modes:
            c = item.channel
            tail = distribution(
                4 * c.coupling**2, item.phonon_energy_eV, 0.0232, 20
            ).weights[-1]
            ct, st = np.cos(c.mixing_angle_rad / 2), np.sin(c.mixing_angle_rad / 2)
            coeff = ct**8 + st**8 + 2 * ct**4 * st**4
            correction += (
                c.amplitude
                * tail
                * coeff
                * item.linewidth_eV
                / (
                    (x - c.exciton_energy_eV - 20 * item.phonon_energy_eV) ** 2
                    + item.linewidth_eV**2
                )
            )
        expected = original.Alpha_fit_gamma12(
            x,
            0.033,
            0.0295,
            0.012,
            0.0201,
            0.0232,
            0.6,
            0.85,
            1.281,
            0.45,
            0.48,
            1.34,
            0.897,
            0.51,
        )
        np.testing.assert_allclose(
            expanded - correction, expected, rtol=3e-13, atol=3e-13
        )

    def test_components_chunking_and_preparation(self):
        x = np.linspace(0.85, 1.45, 400).reshape(20, 20)
        params = ANCHORS["II"]
        modes = [
            VibronicChannel(c, params["phonon_energy_eV"], params["linewidth_eV"])
            for c in params["channels"]
        ]
        prepared = prepare_spectrum(modes, "II", params["kBT_eV"], cutoff=None)
        np.testing.assert_allclose(
            prepared.evaluate(x, max_working_mb=0.001), prepared.evaluate(x), rtol=1e-14
        )
        np.testing.assert_allclose(
            prepared.evaluate(x, components=True).sum(axis=0),
            prepared.evaluate(x),
            rtol=1e-14,
        )
        self.assertEqual(np.shape(prepared.evaluate(1.281)), ())
        c = params["channels"][0]
        pieces = model_ii_components(
            x,
            c,
            params["phonon_energy_eV"],
            params["linewidth_eV"],
            params["kBT_eV"],
            cutoff=None,
        )
        np.testing.assert_allclose(
            sum(pieces.values()),
            model_ii(
                x,
                c,
                params["phonon_energy_eV"],
                params["linewidth_eV"],
                params["kBT_eV"],
                cutoff=None,
            ),
            rtol=2e-14,
        )
        with self.assertRaises(ValueError):
            prepared.kernels[0][1][0] = 3.0

    def test_synthetic_fit_and_metrics(self):
        x, y, sigma, truth = synthetic_case("II")
        fit = fit_fixed_energies(
            x, y, FitSpec("II"), sigma=sigma, starts=1, max_nfev=150
        )
        self.assertTrue(fit["optimizer"]["success"])
        self.assertTrue(all(v >= 0 for v in fit["amplitudes"]))
        self.assertLess(
            np.sqrt(np.mean((np.asarray(fit["prediction"]) - truth) ** 2)), 0.025
        )
        np.testing.assert_allclose(predict_fit(x, fit), fit["prediction"], rtol=1e-12)
        expected_chi_square = float(
            np.sum(((y - np.asarray(fit["prediction"])) / sigma) ** 2)
        )
        self.assertAlmostEqual(
            fit["metrics"]["chi_square"], expected_chi_square, places=8
        )
        unweighted = fit_fixed_energies(x, y, FitSpec("II"), starts=1, max_nfev=150)
        self.assertIsNone(unweighted["metrics"]["chi_square"])

    def test_fit_input_and_named_configuration(self):
        x, y, _, _ = synthetic_case("I")
        with self.assertRaises(ValueError):
            fit_fixed_energies(x, y, initial={"bad_name": 1})
        with self.assertRaises(ValueError):
            fit_fixed_energies(x, y, sigma=0)
        with self.assertRaises(ValueError):
            fit_fixed_energies(x, y, bounds={"linewidth_eV": [0, 0.1]})
        spec = FitSpec("II", (1.281, 0.897), False, False, (0.033, 0.0295))
        names, p0, _, _ = spec.defaults()
        self.assertNotIn("phonon_energy_0_eV", names)
        self.assertIn("linewidth_1_eV", names)
        _, modes = spec.decode(p0)
        self.assertEqual([m.phonon_energy_eV for m in modes], [0.033, 0.0295])


if __name__ == "__main__":
    unittest.main()
