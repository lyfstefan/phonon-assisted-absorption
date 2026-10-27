"""Bounded measured costs; no fitting or general-domain speedup claim."""

from pathlib import Path
import json, time, platform
import numpy as np
import scipy
from phonon_absorption import VibronicChannel, prepare_spectrum, spectrum
from phonon_absorption.parameters import ANCHORS
from phonon_absorption import _legacy_model_i as old_i, _legacy_model_ii as old_ii


def main():
    x = np.linspace(0.85, 1.45, 512)
    report = {
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "scope": "512 energies, two channels, rounded anchors, 40 warm repetitions; excludes imports; same parameters repeated; not whole-domain or fitting speedup",
        "models": {},
    }

    def measure(fn):
        fn()
        times = []
        for _ in range(40):
            t = time.perf_counter()
            fn()
            times.append(time.perf_counter() - t)
        return {
            "median_seconds": float(np.median(times)),
            "p95_seconds": float(np.percentile(times, 95)),
        }

    for model in ["I", "II"]:
        p = ANCHORS[model]
        modes = [
            VibronicChannel(c, p["phonon_energy_eV"], p["linewidth_eV"])
            for c in p["channels"]
        ]
        prepared = prepare_spectrum(modes, model, p["kBT_eV"], cutoff=None)

        def original():
            total = np.zeros_like(x)
            for c in p["channels"]:
                args = (
                    c.coupling,
                    p["phonon_energy_eV"],
                    c.exciton_energy_eV,
                    p["linewidth_eV"],
                    np.sqrt(c.amplitude),
                    p["kBT_eV"],
                )
                total += (
                    old_i.alpha(x, *args)
                    if model == "I"
                    else old_ii.alpha(x, c.mixing_angle_rad, *args)
                )
            return total

        report["models"][model] = {
            "source_functions_native_cutoff": measure(original),
            "stable_fixed80": measure(lambda: spectrum(x, model=model, **p)),
            "stable_adaptive_including_preparation": measure(
                lambda: spectrum(x, model=model, cutoff=None, **p)
            ),
            "prepared_adaptive_evaluation_only": measure(lambda: prepared.evaluate(x)),
            "adaptive_diagnostics": list(prepared.diagnostics),
        }
    Path(__file__).with_name("report.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
