"""Explicit fit/demo commands; no source-script import or optimization on import."""

import argparse, json
from pathlib import Path
import numpy as np
from .fitting import FitSpec, fit_fixed_energies
from .output import write_json


def main():
    parser = argparse.ArgumentParser(
        description="Polaron absorption models and bounded spectral fitting"
    )
    parser.add_argument("command", choices=["demo", "fit"])
    parser.add_argument(
        "--data",
        type=Path,
        help="CSV with energy_eV, absorption and optional sigma columns",
    )
    parser.add_argument(
        "--config", type=Path, help="JSON with spec and optional fit settings"
    )
    parser.add_argument("--output", type=Path, default=Path("examples/output"))
    parser.add_argument(
        "--font-path", type=Path, help="legacy option; LaTeX selects Computer Modern"
    )
    args = parser.parse_args()
    if args.command == "demo":
        from .demo import run

        print(json.dumps(run(args.output, font_path=args.font_path), indent=2))
        return
    if args.data is None or args.config is None:
        parser.error("fit requires --data and --config")
    data = np.atleast_1d(np.genfromtxt(args.data, delimiter=",", names=True))
    if data.dtype.names is None or not {"energy_eV", "absorption"} <= set(
        data.dtype.names
    ):
        parser.error("CSV requires energy_eV and absorption columns")
    config = json.loads(args.config.read_text())
    allowed = {
        "spec",
        "initial",
        "bounds",
        "starts",
        "seed",
        "max_nfev",
        "cutoff",
        "tail_tolerance",
    }
    if set(config) - allowed:
        parser.error("unknown configuration keys: " + str(set(config) - allowed))
    spec = FitSpec(**config.pop("spec", {}))
    fit = fit_fixed_energies(
        data["energy_eV"],
        data["absorption"],
        spec,
        sigma=data["sigma"] if "sigma" in data.dtype.names else None,
        **config,
    )
    args.output.mkdir(parents=True, exist_ok=True)
    write_json(args.output / "fit.json", fit)
    np.savetxt(
        args.output / "fit_at_observations.csv",
        np.column_stack([data["energy_eV"], data["absorption"], fit["prediction"]]),
        delimiter=",",
        header="energy_eV,observed,fit",
        comments="",
    )
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from .plotting import plot_fit

    fig = plot_fit(data["energy_eV"], data["absorption"], fit, font_path=args.font_path)
    from .publication_style import save

    save(fig, args.output / "fit")
    plt.close(fig)
    print(
        json.dumps(
            {"metrics": fit["metrics"], "success": fit["optimizer"]["success"]},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
