"""Figures with explicit data labels and no TeX/global plotting configuration."""

import numpy as np
import matplotlib.pyplot as plt
from .fitting import predict_fit


def energy_windows(energy):
    x = np.sort(np.unique(np.asarray(energy)))
    steps = np.diff(x)
    if len(steps) and steps.max() > 5 * np.median(steps):
        index = int(np.argmax(steps)) + 1
        return [(x[0], x[index - 1]), (x[index], x[-1])]
    return [(x[0], x[-1])]


def plot_fit(
    energy,
    observations,
    fit,
    data_label="Supplied observations",
    title="Fixed-energy spectral fit",
):
    windows = energy_windows(energy)
    fig, axes = plt.subplots(
        2,
        len(windows),
        figsize=(5 * len(windows), 6),
        sharex="col",
        squeeze=False,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )
    x, y = np.asarray(energy), np.asarray(observations)
    pred = np.asarray(fit["prediction"])
    for col, (lo, hi) in enumerate(windows):
        use = (x >= lo) & (x <= hi)
        grid = np.linspace(lo, hi, 400)
        axes[0, col].scatter(
            x[use], y[use], s=12, c="black", alpha=0.6, label=data_label
        )
        axes[0, col].plot(
            grid,
            predict_fit(grid, fit),
            color="tab:red",
            lw=1.8,
            label=f"Model {fit['model']} fit",
        )
        axes[1, col].scatter(x[use], (y - pred)[use], s=12, color="tab:blue")
        axes[1, col].axhline(0, color="gray", lw=0.8)
        axes[0, col].set(ylabel="Absorption (model units)")
        axes[1, col].set(xlabel="Photon energy (eV)", ylabel="Residual")
        axes[0, col].legend(fontsize=8)
    fig.suptitle(
        f"{title}\nMSE={fit['metrics']['mse']:.6g}; fixed exciton energies; nonnegative profiled amplitudes",
        fontsize=11,
    )
    return fig
