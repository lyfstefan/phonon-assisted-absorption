"""Figures with explicit data labels and no TeX/global plotting configuration."""

import numpy as np
import matplotlib.pyplot as plt
from .fitting import predict_fit
from .publication_style import apply_style, width, C


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
    font_path=None,
):
    apply_style(font_path=font_path)
    windows = energy_windows(energy)
    fig, axes = plt.subplots(
        2,
        len(windows),
        figsize=(width(2 if len(windows) > 1 else 1), 4.4),
        sharex="col",
        squeeze=False,
        constrained_layout=False,
        gridspec_kw={"height_ratios": [3, 1]},
    )
    fig.subplots_adjust(
        left=0.11, right=0.985, bottom=0.13, top=0.77, wspace=0.31, hspace=0.09
    )
    x, y = np.asarray(energy), np.asarray(observations)
    pred = np.asarray(fit["prediction"])
    for col, (lo, hi) in enumerate(windows):
        use = (x >= lo) & (x <= hi)
        grid = np.linspace(lo, hi, 400)
        axes[0, col].scatter(
            x[use], y[use], s=10, c=C["black"], alpha=0.8, label=data_label
        )
        axes[0, col].plot(
            grid,
            predict_fit(grid, fit),
            color=C["red"],
            lw=1.2,
            label=f"Model {fit['model']} fit",
        )
        axes[1, col].scatter(x[use], (y - pred)[use], s=9, color=C["blue"])
        axes[1, col].axhline(0, color=C["gray"], lw=0.8)
        axes[0, col].set(ylabel="Absorption (model units)")
        axes[1, col].set(xlabel="Photon energy (eV)", ylabel="Residual")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.53, 0.89),
        ncol=2 if len(windows) > 1 else 1,
        frameon=False,
    )
    fig.suptitle(
        f"{title}\nMSE = {fit['metrics']['mse']:.6g}",
    )
    return fig
