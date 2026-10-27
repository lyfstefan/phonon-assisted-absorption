"""Rounded table parameters, theoretical curves only; no experimental points."""

import argparse
from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from phonon_absorption import spectrum
from phonon_absorption.parameters import ANCHORS
from phonon_absorption.publication_style import apply_style, width, save, C, DASH


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("examples/output"))
    parser.add_argument("--font-path", type=Path, help="optional Roman font file")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    apply_style(font_path=args.font_path)
    energy = np.linspace(0.8, 1.6, 1600)
    fig, ax = plt.subplots(figsize=(width(1), 2.65), constrained_layout=True)
    columns = [energy]
    for index, model in enumerate(["I", "II"]):
        y = spectrum(energy, model=model, **ANCHORS[model])
        columns.append(y)
        ax.plot(
            energy,
            y,
            label=f"Model {model}",
            color=[C["blue"], C["red"]][index],
            linestyle=DASH[index],
        )
    ax.set(
        xlabel="Photon energy (eV)",
        ylabel="Absorption (model units)",
    )
    fig.suptitle("Theory at rounded parameters")
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2)
    save(fig, args.output / "theory_spectra")
    fig.savefig(args.output / "theory_spectra.svg")
    plt.close(fig)
    np.savetxt(
        args.output / "theory_spectra.csv",
        np.column_stack(columns),
        delimiter=",",
        header="energy_eV,model_I,model_II",
        comments="",
    )
    print(args.output / "theory_spectra.png")


if __name__ == "__main__":
    main()
