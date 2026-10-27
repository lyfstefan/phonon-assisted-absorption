# Research source and version choices

The author identified `Version_2_removing_middle` as the phonon project
folder. All three Python fitting variants in that folder were read without
importing or executing them:

| Source | Extracted reference module | Supported assumptions |
|---|---|---|
| `fitting_both_KK_and_KGamma_single.py` | `_legacy_model_i.py` | Model I |
| `fitting_both_KK_and_KGamma_single_interlayer_same_freq.py` | `_legacy_model_ii.py` | Model II, shared frequency/linewidth |
| `fitting_both_KK_and_KGamma_single_interlayer_diff_freq.py` | `_legacy_model_ii_independent.py` | Model II, independent frequencies and optional independent linewidths |

The last file has a top-level fitting call. Its nested
`Alpha_fit_gamma12` function takes all inputs explicitly and captures no
local variables; it was extracted as a standalone reference definition.
Reference modules retain selected function source text, removing only
trailing whitespace, and add required imports/provenance docstrings.
[source_manifest.json](source_manifest.json) records the source hashes and
selected functions. No embedded observations, fitting/plotting calls,
timestamps, or global Matplotlib configuration were copied into them.
All original Physics files remain unchanged.

No supplied evidence establishes which variant generated the final paper
figures. These remain candidate historical implementations. The separate
Version 1, FFT/multimode routes, and full paper-figure reproduction are not
part of this release.

## Changes made for version 0.2

The public forward model replaces `exp * iv` thermal factors with equivalent
Skellam weights, including analytic zero-temperature and zero-coupling
limits. It adds adaptive truncation, diagnostics, blocked evaluation, and
reusable prepared kernels. Tests compare the formulas using the exact
historical asymmetric sideband supports before testing convergence with
larger supports. The archived reference functions can overflow and are not
the recommended public API.

The new optimizer uses bounded nonlinear least squares and solves the
nonnegative linear squared amplitudes at every evaluation. This replaces
the historical `curve_fit` procedure and its squared-`wf` parameters.
Multiple seeded starting points and explicit named configurations make
new fits reproducible. Exciton energies are fixed. No claim is made of a
global optimum, unique physical parameters, or historical fit recovery.

Public demonstrations use independently generated synthetic observations
or theory curves. Local fitting of the scripts' embedded comparison points
produced only the aggregate public `validation/local_refit_summary.json`;
point arrays, fitted values at those points, residuals, and experimental
overlay figures stay outside this repository. Those points were used for
comparison with published experiments and are not the author's raw
measurements. Their exact digitization and redistribution status remain
unestablished. The unresolved original rounded-parameter MSE discrepancy
is separate from the newly optimized fit statistics.
