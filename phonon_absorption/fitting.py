"""Bounded fixed-exciton-energy fitting with profiled nonnegative amplitudes.

Nonlinear spectral parameters use least_squares; at every evaluation the
linear |psi|^2 amplitudes are solved by NNLS, not optimized as squared wf.
This is an explicitly changed optimizer, not the historical curve_fit path.
"""

from dataclasses import dataclass, asdict
import numpy as np
from scipy.optimize import least_squares, nnls
from .models import Channel, VibronicChannel, prepare_spectrum, spectrum_independent


@dataclass(frozen=True)
class FitSpec:
    model: str = "II"
    exciton_energies_eV: tuple = (1.281, 0.898)
    shared_phonon_energy: bool = True
    shared_linewidth: bool = True
    fixed_phonon_energies_eV: tuple | None = None

    def __post_init__(self):
        object.__setattr__(self, "exciton_energies_eV", tuple(self.exciton_energies_eV))
        if self.model not in ("I", "II") or not self.exciton_energies_eV:
            raise ValueError("choose model I/II and at least one fixed exciton energy")
        if not isinstance(self.shared_phonon_energy, bool) or not isinstance(
            self.shared_linewidth, bool
        ):
            raise ValueError("shared parameter flags must be booleans")
        if not np.isfinite(self.exciton_energies_eV).all():
            raise ValueError("exciton energies must be finite")
        if self.fixed_phonon_energies_eV is not None:
            values = tuple(self.fixed_phonon_energies_eV)
            if (
                len(values) != len(self.exciton_energies_eV)
                or not np.isfinite(values).all()
                or min(values) <= 0
            ):
                raise ValueError("provide one positive fixed phonon energy per exciton")
            object.__setattr__(self, "fixed_phonon_energies_eV", values)

    def defaults(self):
        """Named computational bounds; not experimentally established priors."""
        names, initial, lower, upper = [], [], [], []

        def add(name, value, lo, hi):
            names.append(name)
            initial.append(value)
            lower.append(lo)
            upper.append(hi)

        n = len(self.exciton_energies_eV)
        if self.fixed_phonon_energies_eV is None:
            for i in range(1 if self.shared_phonon_energy else n):
                add(
                    (
                        "phonon_energy_eV"
                        if self.shared_phonon_energy
                        else f"phonon_energy_{i}_eV"
                    ),
                    0.0346 if self.model == "I" else 0.032,
                    0.005,
                    0.1,
                )
        for i in range(1 if self.shared_linewidth else n):
            add(
                "linewidth_eV" if self.shared_linewidth else f"linewidth_{i}_eV",
                0.0241 if self.model == "I" else 0.0231,
                0.0005,
                0.1,
            )
        add("kBT_eV", 0.0284 if self.model == "I" else 0.0232, 0.0001, 0.15)
        for i in range(n):
            coupling = ([1.5755, 2.6136] if self.model == "I" else [0.85, 1.34])[i % 2]
            add(f"coupling_{i}", coupling, 0.0, 5.0)
        if self.model == "II":
            for i in range(n):
                add(f"mixing_angle_{i}_rad", [0.6, 0.48][i % 2], 0.0, np.pi / 2)
        return names, np.array(initial), np.array(lower), np.array(upper)

    def decode(self, values, amplitudes=None):
        names, _, _, _ = self.defaults()
        parameters = dict(zip(names, values))
        modes = []
        for i, energy in enumerate(self.exciton_energies_eV):
            omega = (
                self.fixed_phonon_energies_eV[i]
                if self.fixed_phonon_energies_eV is not None
                else parameters[
                    (
                        "phonon_energy_eV"
                        if self.shared_phonon_energy
                        else f"phonon_energy_{i}_eV"
                    )
                ]
            )
            linewidth = parameters[
                "linewidth_eV" if self.shared_linewidth else f"linewidth_{i}_eV"
            ]
            channel = Channel(
                energy,
                parameters[f"coupling_{i}"],
                1.0 if amplitudes is None else float(amplitudes[i]),
                parameters.get(f"mixing_angle_{i}_rad", 0.0),
            )
            modes.append(VibronicChannel(channel, omega, linewidth))
        return parameters, modes


def fit_fixed_energies(
    energy_eV,
    observed,
    spec=None,
    sigma=None,
    initial=None,
    bounds=None,
    starts=3,
    seed=2026,
    max_nfev=300,
    cutoff=None,
    tail_tolerance=1e-11,
):
    """Fit supplied data without inventing samples between energy windows.

    Without supplied measurement uncertainties, report RSS/MSE only; no
    chi-square p-values or uncertainty estimates are manufactured. A weighted
    chi-square is reported only when positive sigma values are supplied.
    Local multistart convergence does not certify a unique physical fit.
    """
    spec = spec or FitSpec()
    x, y = np.asarray(energy_eV, dtype=float), np.asarray(observed, dtype=float)
    names, p0, lo, hi = spec.defaults()
    if (
        x.ndim != 1
        or x.shape != y.shape
        or not np.isfinite(x).all()
        or not np.isfinite(y).all()
    ):
        raise ValueError(
            "finite one-dimensional energy and observation arrays of equal size required"
        )
    if len(x) <= len(names) + len(spec.exciton_energies_eV):
        raise ValueError("more observations than fitted parameters are required")
    if (
        not isinstance(starts, int)
        or starts < 1
        or not isinstance(max_nfev, int)
        or max_nfev < 1
    ):
        raise ValueError("starts and max_nfev must be positive integers")
    if initial is not None:
        if not isinstance(initial, dict) or set(initial) - set(names):
            raise ValueError("initial must map known parameter names to values")
        p0 = np.array([initial.get(name, value) for name, value in zip(names, p0)])
    if bounds is not None:
        if not isinstance(bounds, dict) or set(bounds) - set(names):
            raise ValueError("bounds must map known parameter names to [lower, upper]")
        for i, name in enumerate(names):
            if name in bounds:
                lo[i], hi[i] = bounds[name]
    if (
        not np.isfinite(np.r_[p0, lo, hi]).all()
        or np.any(lo >= hi)
        or np.any(p0 < lo)
        or np.any(p0 > hi)
    ):
        raise ValueError("finite ordered bounds and a feasible initial point required")
    # Bounds must keep every trial inside the forward model's physical domain.
    for i, name in enumerate(names):
        if name.endswith("_eV") and lo[i] <= 0:
            raise ValueError(
                "fitted energies/widths/thermal energy require positive lower bounds"
            )
        if name.startswith("coupling_") and lo[i] < 0:
            raise ValueError("coupling bounds must be nonnegative")
        if name.startswith("mixing_angle_") and (lo[i] < 0 or hi[i] > np.pi / 2):
            raise ValueError("fit angles are restricted to [0, pi/2]")
    uncertainty = (
        np.ones_like(y)
        if sigma is None
        else np.broadcast_to(np.asarray(sigma, dtype=float), y.shape)
    )
    if not np.isfinite(uncertainty).all() or np.any(uncertainty <= 0):
        raise ValueError("sigma must be finite and positive")
    calls = 0

    def projected(values):
        nonlocal calls
        calls += 1
        parameters, modes = spec.decode(values)
        basis = (
            prepare_spectrum(
                modes, spec.model, parameters["kBT_eV"], cutoff, tail_tolerance
            )
            .evaluate(x, components=True)
            .T
        )
        amplitudes, _ = nnls(basis / uncertainty[:, None], y / uncertainty)
        prediction = basis @ amplitudes
        return (prediction - y) / uncertainty, amplitudes, prediction

    rng = np.random.default_rng(seed)
    trials, solutions = [], []
    for index in range(starts):
        start = (
            p0
            if index == 0
            else np.clip(
                p0
                + rng.normal(size=len(p0))
                * 0.25
                * np.maximum(abs(p0), 0.1 * (hi - lo)),
                lo + 1e-8 * (hi - lo),
                hi - 1e-8 * (hi - lo),
            )
        )
        result = least_squares(
            lambda p: projected(p)[0],
            start,
            bounds=(lo, hi),
            x_scale="jac",
            max_nfev=max_nfev,
            ftol=1e-9,
            xtol=1e-9,
            gtol=1e-9,
        )
        residual, amplitudes, prediction = projected(result.x)
        objective = float(residual @ residual)
        trials.append(
            {
                "start": index,
                "initial": dict(zip(names, map(float, start))),
                "success": bool(result.success),
                "status": int(result.status),
                "message": str(result.message),
                "nfev": int(result.nfev),
                "weighted_rss": objective,
            }
        )
        solutions.append((objective, result, amplitudes, prediction))
    objective, result, amplitudes, prediction = min(
        solutions, key=lambda value: value[0]
    )
    parameters, modes = spec.decode(result.x, amplitudes)
    singular = np.linalg.svd(result.jac, compute_uv=False)
    threshold = np.finfo(float).eps * max(result.jac.shape) * singular[0]
    rank = int(np.sum(singular > threshold))
    rss = float(np.sum((prediction - y) ** 2))
    modes_prepared = prepare_spectrum(
        modes, spec.model, parameters["kBT_eV"], cutoff, tail_tolerance
    )
    check = prepare_spectrum(
        modes, spec.model, parameters["kBT_eV"], None, 1e-13
    ).evaluate(x)
    fit = {
        "model": spec.model,
        "spec": asdict(spec),
        "parameters": {k: float(v) for k, v in parameters.items()},
        "channels": [asdict(mode) for mode in modes],
        "amplitudes": amplitudes.tolist(),
        "prediction": prediction.tolist(),
        "residual_observed_minus_fit": (y - prediction).tolist(),
        "metrics": {
            "n_observations": len(x),
            "rss": rss,
            "mse": rss / len(x),
            "weighted_rss": objective,
            "chi_square": objective if sigma is not None else None,
            "nominal_dof": len(x) - len(names) - len(amplitudes),
        },
        "optimizer": {
            "method": "bounded least_squares with NNLS-profiled amplitudes",
            "success": bool(result.success),
            "multistart_seed": seed,
            "trials": trials,
            "forward_evaluations": calls,
            "projected_jacobian_rank": rank,
            "projected_jacobian_columns": len(names),
            "projected_jacobian_singular_values": singular.tolist(),
            "nonlinear_bounds": {
                name: [float(a), float(b)] for name, a, b in zip(names, lo, hi)
            },
        },
        "truncation": {
            "cutoff": cutoff,
            "tail_tolerance": tail_tolerance,
            "channels": list(modes_prepared.diagnostics),
            "relative_l2_change_at_tighter_adaptive_tail": float(
                np.linalg.norm(check - prediction) / max(np.linalg.norm(check), 1e-300)
            ),
        },
        "scope": "New local fit using supplied data and explicit computational bounds; no global-optimum, parameter-identifiability or original-paper-fit claim.",
    }
    return fit


def predict_fit(energy_eV, fit):
    modes = [
        VibronicChannel(
            Channel(**mode["channel"]), mode["phonon_energy_eV"], mode["linewidth_eV"]
        )
        for mode in fit["channels"]
    ]
    return spectrum_independent(
        energy_eV,
        modes,
        fit["model"],
        fit["parameters"]["kBT_eV"],
        cutoff=fit["truncation"]["cutoff"],
        tail_tolerance=fit["truncation"]["tail_tolerance"],
    )
