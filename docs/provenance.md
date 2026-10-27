# Research source and numerical changes

The source audit supplied with this project was read before extraction.
Two files from the author's `Version_2_removing_middle` directory were
selected as forward-model candidates:

1. `fitting_both_KK_and_KGamma_single.py` -> `_legacy_model_i.py`.
2. `fitting_both_KK_and_KGamma_single_interlayer_same_freq.py` ->
   `_legacy_model_ii.py`.

The new modules contain the selected function definitions (only trailing
whitespace is removed),
plus the necessary imports and a provenance docstring. Source SHA-256
fingerprints and selected names are in `source_manifest.json`.
Top-level data arrays, fitting/plotting calls, timestamps and Matplotlib
global configuration were not copied or executed. The original Physics
files were read only and were not edited. Their private absolute paths
are not needed by this repository.

These are candidate versions, not asserted to be the final scripts used to
generate the published figures. The `diff_freq` extension has different
channel assumptions and was not selected for the shared-frequency example.
Version 1 has duplicate definitions and top-level optimization; it was not
imported. The million-point FFT/multimode route was not imported or claimed
to be validated.

## Adapted interface

`models.py` is a new, documented interface adapted from these existing
formulas. The physical expressions and table anchors were checked against
arXiv:2602.02803v2, sections II/III and the local thesis theory and parameter
table. Finite-temperature sideband weights are evaluated through the
equivalent Skellam law rather than the original `exp * iv` expression.
Zero temperature is handled with its one-sided Poisson limit, and zero
coupling with a delta weight at m=0. These are numerical/analytic changes
documented as such, not verbatim recovery of historical implementation.

The wrapper exposes a multiplicative squared amplitude rather than an
ambiguously named `wf`. Default cutoffs are symmetric and larger. Source
equivalence is therefore checked with the original asymmetric cutoffs
before expanded-cutoff spectra are compared. The reference functions may
overflow at low temperature and should not be used as the public API.

The validation report contains no benchmark speedup claim, no refitted
parameters, and no experimental arrays. It verifies only the stated grid,
anchors, thermal probes and limit cases. A globally converged fit, complete
figure reproduction, multi-mode extension, and experimental-data release
remain outside this release.
