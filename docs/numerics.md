# Fitting and numerical details

## Fitting supplied observations

`fit_fixed_energies` takes a one-dimensional energy array, observations,
`FitSpec`, and optional positive standard deviations `sigma`. It evaluates
only the supplied energies, preserving their order. Plotting separates the
two energy windows when it detects the large gap; no missing observations
are synthesized between them.

At every nonlinear trial, unit-amplitude channel spectra form a design
matrix. `scipy.optimize.nnls` solves the linear nonnegative squared
amplitudes. `least_squares` then optimizes phonon energies, linewidths,
thermal energy, coupling magnitudes, and Model II mixing angles. The first
start uses the configured initial values; later starts are seeded
perturbations. The returned solution has the lowest objective among the
tried starts, with success and termination messages reported separately.
Convergence is local and does not prove parameter identifiability.

Configurations accept `spec`, `initial`, `bounds`, `starts`, `seed`,
`max_nfev`, `cutoff`, and `tail_tolerance`. `initial` and `bounds` use named
parameters, for example:

```json
{
  "spec": {"model": "II", "exciton_energies_eV": [1.280999, 0.898]},
  "initial": {"phonon_energy_eV": 0.032},
  "bounds": {"linewidth_eV": [0.001, 0.06]},
  "starts": 5,
  "seed": 2026,
  "cutoff": null
}
```

Shared parameters are `phonon_energy_eV` and `linewidth_eV`; independent
ones are `phonon_energy_0_eV`, `linewidth_0_eV`, etc. Couplings and mixing
angles are named `coupling_0` and `mixing_angle_0_rad`. Thermal energy is
`kBT_eV`. Use `FitSpec.defaults()` to inspect exact names and starts.
The default computational bounds are phonons [0.005, 0.1] eV,
linewidths [0.0005, 0.1] eV, kBT [0.0001, 0.15] eV, couplings [0, 5], and
angles [0, pi/2]. These are search bounds, not experimentally justified
priors. Initial angles are away from the zero-derivative endpoint.

Outputs include parameters, squared amplitudes, predictions, residuals,
RSS/MSE, trial statuses, forward-call count, bounds, and projected Jacobian
singular values. Singular values do not establish physical identifiability.
`chi_square` is null without supplied sigma. With sigma, weighted RSS is
reported as chi-square, meaningful as such only when sigma represents
known measurement standard deviations and the assumed error model is
appropriate. No p-values or parameter error bars are inferred; the reported
degrees of freedom are explicitly nominal.

## Sideband support and stability

The thermal law is Skellam with means `S*(nB+1)` and `S*nB`, with
`S=g_cv**2` in Model I and `S=4*g_norm**2` in Model II. At zero temperature
it reduces to Poisson; at zero coupling, to a delta at m=0. Truncated
weights are never renormalized.

With `cutoff=None`, an initial mean/variance bracket expands until the
computed CDF/SF omitted tail meets `tail_tolerance`. No Gaussian weights
are used. If SciPy's Skellam routine overflows, a conservative union bound
from the two independent Poisson variables estimates the tail. PMF
fallbacks use a log-space hypergeometric expression or scaled Bessel
functions. Diagnostics name the method used. Numerical tail estimates are
not certified arithmetic. An explicit term-count guard rejects excessively
large requested supports rather than silently truncating them.

Each prepared channel reports retained weight, omitted thermal mass, and
an absolute spectral truncation-bound estimate based on the Lorentzian
maximum and absolute coefficients. These estimates cover omitted
sidebands, not all floating-point errors. Fits recheck their predictions
with a tighter adaptive tail of 1e-13. The default public forward interface
retains `cutoff=80` for compatibility; fitting defaults to adaptive support.

`prepare_spectrum` stores fixed weights/centers for reuse on many energy
grids. Its `components=True` evaluation returns independent exciton
channels. `model_ii_components` instead returns the three algebraic terms;
the interference term may be signed. Evaluation uses blocks controlled by
`max_working_mb`, an estimate of temporary array space rather than a total
memory limit.

## Measured performance

Run `python -m benchmarks.run` to regenerate
[benchmarks/report.json](../benchmarks/report.json). The recorded run uses
Python 3.11.5, NumPy 1.26.4, SciPy 1.15.1, 512 energies, two channels,
rounded anchor parameters, and 40 warm repetitions. Median times:

| Evaluation | Model I (ms) | Model II (ms) |
|---|---:|---:|
| Historical functions, native cutoffs | 0.878 | 0.785 |
| Stable implementation, fixed cutoff 80 | 0.822 | 0.867 |
| Adaptive support including preparation | 0.770 | 0.832 |
| Prepared adaptive kernels, evaluation only | 0.107 | 0.132 |

Historical supports differ from the expanded supports. Prepared timings
exclude kernel construction and apply only when parameters remain fixed.
There is no uniform single-call speedup: the main measured gain is reuse
for repeated spectra. These timings exclude imports and do not measure
end-to-end nonlinear fitting or establish a whole-domain speedup.
