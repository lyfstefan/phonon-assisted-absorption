import json
import platform
from pathlib import Path
import numpy as np
import scipy
from phonon_absorption import Channel, model_i, model_ii, spectrum, sideband_weights
from phonon_absorption.parameters import ANCHORS
from phonon_absorption import _legacy_model_i as old_i, _legacy_model_ii as old_ii
from phonon_absorption.thermal import distribution
from phonon_absorption.fitting import FitSpec, fit_fixed_energies
from phonon_absorption.demo import synthetic_case


def main():
    grid = np.linspace(0.85, 1.45, 512)
    report = {
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "grid": {"points": 512, "min_eV": 0.85, "max_eV": 1.45},
        "models": {},
        "low_temperature_checks": [],
        "limit_checks": [],
    }
    rel = lambda a, b: float(np.linalg.norm(a - b) / np.linalg.norm(a))
    for name in ["I", "II"]:
        params = ANCHORS[name]
        public = spectrum(grid, model=name, **params)
        legacy = np.zeros_like(grid)
        for channel in params["channels"]:
            args = (
                channel.coupling,
                params["phonon_energy_eV"],
                channel.exciton_energy_eV,
                params["linewidth_eV"],
                np.sqrt(channel.amplitude),
                params["kBT_eV"],
            )
            legacy += (
                old_i.alpha(grid, *args)
                if name == "I"
                else old_ii.alpha(grid, channel.mixing_angle_rad, *args)
            )
        # Match the exact asymmetric cutoffs of the historical scripts in the
        # stable evaluator, so equivalence and cutoff error are not conflated.
        matched = np.zeros_like(grid)
        cut = 50 if name == "I" else 20
        for c in params["channels"]:
            fn = model_i if name == "I" else model_ii
            full = fn(
                grid,
                c,
                params["phonon_energy_eV"],
                params["linewidth_eV"],
                params["kBT_eV"],
                cut,
            )
            S = c.coupling**2 * (1 if name == "I" else 4)
            m, w = sideband_weights(
                S, params["phonon_energy_eV"], params["kBT_eV"], cut
            )
            coeff = 1.0
            if name == "II":
                ct, st = np.cos(c.mixing_angle_rad / 2), np.sin(c.mixing_angle_rad / 2)
                coeff = ct**8 + st**8 + (-1.0) ** cut * 2 * ct**4 * st**4
            lor = params["linewidth_eV"] / (
                (grid - c.exciton_energy_eV - cut * params["phonon_energy_eV"]) ** 2
                + params["linewidth_eV"] ** 2
            )
            matched += full - c.amplitude * coeff * w[-1] * lor
        expanded = spectrum(grid, model=name, cutoff=120, **params)
        error = rel(legacy, matched)
        assert error < 1e-12
        assert rel(expanded, public) < 1e-10
        assert np.isfinite(public).all() and np.min(public) >= 0
        report["models"][name] = {
            "matched_historical_cutoff_relative_L2": error,
            "historical_vs_default_cutoff_relative_L2": rel(public, legacy),
            "cutoff80_vs120_relative_L2": rel(expanded, public),
            "min_spectrum": float(public.min()),
            "max_spectrum": float(public.max()),
        }
    for g in [0, 0.4, 0.85, 1.34]:
        for kBT in [0, 0.001, 0.0232, 0.05]:
            a = model_i(grid, Channel(1.281, 2 * g, 1), 0.032, 0.0231, kBT)
            b = model_ii(grid, Channel(1.281, g, 1, 0), 0.032, 0.0231, kBT)
            err = rel(a, b)
            assert err < 1e-13
            report["limit_checks"].append(
                {"g_norm": g, "kBT_eV": kBT, "theta_zero_gcv_2g_relative_L2": err}
            )
    for temp in [0, 0.001, 20 * 8.617333262145e-5, 0.0232]:
        for name in ["I", "II"]:
            params = {**ANCHORS[name], "kBT_eV": temp}
            y = spectrum(grid, model=name, **params)
            assert np.isfinite(y).all() and np.min(y) >= 0
            report["low_temperature_checks"].append(
                {"model": name, "kBT_eV": temp, "finite_nonnegative": True}
            )
    # Thermal weights satisfy detailed balance P(-m)/P(m)=exp(-m Omega/kBT).
    m, w = sideband_weights(3, 0.032, 0.0232, 80)
    balance = float(
        np.max(
            abs(
                w[80 - np.arange(1, 8)] / w[80 + np.arange(1, 8)]
                - np.exp(-np.arange(1, 8) * 0.032 / 0.0232)
            )
        )
    )
    assert balance < 1e-13
    report["thermal_weight_sum"] = float(w.sum())
    report["max_detailed_balance_error"] = balance
    report["adaptive_checks"] = []
    xcheck = np.linspace(0.85, 1.45, 128)
    for factor in [0, 0.01, 1, 7, 25, 144]:
        for temperature in [0, 0.001, 0.0232, 0.1]:
            weights = distribution(factor, 0.032, temperature, cutoff=None)
            assert (
                weights.omitted_mass <= 1e-12 and abs(weights.weights.sum() - 1) < 1e-10
            )
            for model in ["I", "II"]:
                channel = Channel(
                    1.281, np.sqrt(factor) / (1 if model == "I" else 2), 1.0, 0.4
                )
                fn = model_i if model == "I" else model_ii
                adaptive = fn(xcheck, channel, 0.032, 0.0231, temperature, cutoff=None)
                expanded = fn(xcheck, channel, 0.032, 0.0231, temperature, cutoff=1200)
                error = rel(expanded, adaptive)
                assert error < 1e-9
                report["adaptive_checks"].append(
                    {
                        "model": model,
                        "huang_rhys": factor,
                        "kBT_eV": temperature,
                        **weights.diagnostics(),
                        "relative_l2_vs_cutoff1200": error,
                    }
                )
    report["synthetic_fitting_checks"] = {}
    for name in ["I", "II"]:
        x, y, sigma, truth = synthetic_case(name)
        energies = tuple(c.exciton_energy_eV for c in ANCHORS[name]["channels"])
        fit = fit_fixed_energies(x, y, FitSpec(name, energies), sigma=sigma, starts=3)
        rmse = float(np.sqrt(np.mean((np.asarray(fit["prediction"]) - truth) ** 2)))
        assert fit["optimizer"]["success"] and rmse < 0.025
        report["synthetic_fitting_checks"][name] = {
            "data_origin": "seeded synthetic theory data only",
            "seed": 2026,
            "noise_sigma": 0.02,
            "n_observations": 120,
            "mse": fit["metrics"]["mse"],
            "rmse_vs_generating_theory": rmse,
            "nonnegative_amplitudes": bool(min(fit["amplitudes"]) >= 0),
            "optimizer_success": fit["optimizer"]["success"],
            "truncation_recheck": fit["truncation"],
        }
    report["scope"] = (
        "Source/limit comparisons, bounded adaptive-tail probes and synthetic fitting checks. No experimental point release, universal stability proof or original-paper-fit recovery claim."
    )
    report["all_checks_passed"] = True
    Path(__file__).with_name("report.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps({"all_checks_passed": True, "models": report["models"]}, indent=2))


if __name__ == "__main__":
    main()
