"""Reproducible synthetic fitting and theory examples, with no experimental data."""

from pathlib import Path
import numpy as np
from .models import Channel, VibronicChannel, spectrum, model_ii_components
from .parameters import ANCHORS
from .fitting import FitSpec, fit_fixed_energies
from .output import write_json
from .publication_style import apply_style, width, save, C, DASH


def synthetic_case(model, seed=2026):
    energy = np.r_[np.linspace(0.87, 1.05, 60), np.linspace(1.23, 1.43, 60)]
    parameters = ANCHORS[model]
    truth = spectrum(energy, model=model, cutoff=None, **parameters)
    sigma = np.full_like(energy, 0.02)
    observed = truth + np.random.default_rng(seed).normal(0, sigma)
    return energy, observed, sigma, truth


def run(output, font_path=None):
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from .plotting import plot_fit

    apply_style(font_path=font_path)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    summaries = {}
    for model in ["I", "II"]:
        x, y, sigma, truth = synthetic_case(model)
        energies = tuple(c.exciton_energy_eV for c in ANCHORS[model]["channels"])
        fit = fit_fixed_energies(x, y, FitSpec(model, energies), sigma=sigma, starts=3)
        fit["data_origin"] = (
            "Synthetic theory observations with seeded Gaussian noise; no experimental data"
        )
        write_json(output / f"synthetic_fit_{model}.json", fit)
        np.savetxt(
            output / f"synthetic_observations_{model}.csv",
            np.column_stack([x, y, sigma, truth]),
            delimiter=",",
            header="energy_eV,absorption,sigma,generating_theory",
            comments="",
        )
        fig = plot_fit(
            x,
            y,
            fit,
            "Synthetic observations",
            f"Synthetic validation: Model {model}",
            font_path=font_path,
        )
        apply_style(font_path=font_path)
        save(fig, output / f"synthetic_fit_{model}")
        fig.savefig(output / f"synthetic_fit_{model}.svg")
        plt.close(fig)
        summaries[model] = {
            "mse": fit["metrics"]["mse"],
            "success": fit["optimizer"]["success"],
            "rmse_against_generating_theory": float(
                np.sqrt(np.mean((np.array(fit["prediction"]) - truth) ** 2))
            ),
            "noise_sigma": 0.02,
            "n_observations": len(x),
        }
    # Explicit theory temperature dependence, using the same table anchors.
    fig, axes = plt.subplots(2, 2, figsize=(width(2), 5.2), constrained_layout=True)
    for row, model in enumerate(["I", "II"]):
        for col, (lo, hi) in enumerate([(0.84, 1.08), (1.2, 1.48)]):
            energy = np.linspace(lo, hi, 500)
            for ti, temperature in enumerate([0, 20, 100, 270]):
                params = {**ANCHORS[model], "kBT_eV": temperature * 8.617333262145e-5}
                axes[row, col].plot(
                    energy,
                    spectrum(energy, model=model, cutoff=None, **params),
                    label=f"{temperature} K",
                    color=[C["blue"], C["orange"], C["green"], C["red"]][ti],
                    linestyle=DASH[ti],
                )
            axes[row, col].set(
                title=f"Model {model}",
                xlabel="Photon energy (eV)",
                ylabel="Absorption (model units)",
            )
            axes[row, col].legend(frameon=False)
    fig.suptitle("Theory: temperature dependence")
    apply_style(font_path=font_path)
    save(fig, output / "temperature_spectra")
    fig.savefig(output / "temperature_spectra.svg")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(width(2), 3.1), constrained_layout=True)
    energy = np.linspace(0.8, 1.6, 1200)
    params = ANCHORS["II"]
    total = np.zeros_like(energy)
    for ci, (label, ch) in enumerate(
        zip(["KK channel", r"$\Gamma$-K channel"], params["channels"])
    ):
        pieces = model_ii_components(
            energy,
            ch,
            params["phonon_energy_eV"],
            params["linewidth_eV"],
            params["kBT_eV"],
            cutoff=None,
        )
        channel = sum(pieces.values())
        total += channel
        ax.plot(
            energy,
            channel,
            linestyle=DASH[ci + 1],
            color=[C["blue"], C["orange"]][ci],
            label=label,
        )
    ax.plot(energy, total, c=C["black"], label="Model II total")
    ax.set(
        xlabel="Photon energy (eV)",
        ylabel="Absorption (model units)",
        title="Theory-only exciton-channel decomposition",
    )
    ax.legend(frameon=False)
    apply_style(font_path=font_path)
    save(fig, output / "channel_decomposition")
    fig.savefig(output / "channel_decomposition.svg")
    plt.close(fig)
    write_json(
        output / "synthetic_validation.json",
        {"seed": 2026, "data_origin": "synthetic only", "models": summaries},
    )
    return summaries
