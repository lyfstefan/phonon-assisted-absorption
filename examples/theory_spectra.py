"""Rounded table parameters, theoretical curves only; no experimental points."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from phonon_absorption import spectrum
from phonon_absorption.parameters import ANCHORS


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=Path("examples/output"))
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    energy=np.linspace(.8,1.6,1600)
    fig,ax=plt.subplots(figsize=(8,4.5),constrained_layout=True)
    columns=[energy]
    for model in ["I","II"]:
        y=spectrum(energy,model=model,**ANCHORS[model])
        columns.append(y);ax.plot(energy,y,label=f"Model {model}")
    ax.set(xlabel="Photon energy (eV)",ylabel="Absorption (model units)",title="Phonon-assisted absorption: theory at rounded table parameters")
    ax.legend();fig.savefig(args.output/"theory_spectra.png",dpi=180);plt.close(fig)
    np.savetxt(args.output/"theory_spectra.csv",np.column_stack(columns),delimiter=",",header="energy_eV,model_I,model_II",comments="")
    print(args.output/"theory_spectra.png")


if __name__=="__main__":main()
