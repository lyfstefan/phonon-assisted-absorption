"""Reproducible synthetic fitting and theory examples, with no experimental data."""

from pathlib import Path
import numpy as np
from .models import Channel, VibronicChannel, spectrum, model_ii_components
from .parameters import ANCHORS
from .fitting import FitSpec, fit_fixed_energies
from .output import write_json


def synthetic_case(model, seed=2026):
    energy = np.r_[np.linspace(0.87, 1.05, 60), np.linspace(1.23, 1.43, 60)]
    parameters = ANCHORS[model]
    truth = spectrum(energy, model=model, cutoff=None, **parameters)
    sigma = np.full_like(energy, 0.02)
    observed = truth + np.random.default_rng(seed).normal(0, sigma)
    return energy, observed, sigma, truth


def run(output):
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from .plotting import plot_fit

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
            x, y, fit, "Synthetic observations", f"Synthetic validation: Model {model}"
        )
        for suffix in ["png", "svg"]:
            fig.savefig(output / f"synthetic_fit_{model}.{suffix}", dpi=170)
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
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for row, model in enumerate(["I", "II"]):
        for col, (lo, hi) in enumerate([(0.84, 1.08), (1.2, 1.48)]):
            energy = np.linspace(lo, hi, 500)
            for temperature in [0, 20, 100, 270]:
                params = {**ANCHORS[model], "kBT_eV": temperature * 8.617333262145e-5}
                axes[row, col].plot(
                    energy,
                    spectrum(energy, model=model, cutoff=None, **params),
                    label=f"{temperature} K",
                )
            axes[row, col].set(
                title=f"Model {model}",
                xlabel="Photon energy (eV)",
                ylabel="Absorption (model units)",
            )
            axes[row, col].legend(fontsize=8)
    fig.suptitle(
        "Theory-only temperature dependence at fixed nonthermal table parameters"
    )
    for suffix in ["png", "svg"]:
        fig.savefig(output / f"temperature_spectra.{suffix}", dpi=170)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 4), constrained_layout=True)
    energy = np.linspace(0.8, 1.6, 1200)
    params = ANCHORS["II"]
    total = np.zeros_like(energy)
    for label, ch in zip(["KK channel", "Gamma-K channel"], params["channels"]):
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
        ax.plot(energy, channel, "--", label=label)
    ax.plot(energy, total, c="black", label="Model II total")
    ax.set(
        xlabel="Photon energy (eV)",
        ylabel="Absorption (model units)",
        title="Theory-only exciton-channel decomposition",
    )
    ax.legend()
    for suffix in ["png", "svg"]:
        fig.savefig(output / f"channel_decomposition.{suffix}", dpi=170)
    plt.close(fig)
    write_json(
        output / "synthetic_validation.json",
        {"seed": 2026, "data_origin": "synthetic only", "models": summaries},
    )
    return summaries
