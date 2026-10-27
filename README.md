# Phonon-assisted absorption in TMD heterostructures

Absorption spectra and spectral fitting for phonon-coupled excitons in
transition metal dichalcogenide heterostructures. The calculations accompany
Y. Liu, R. Dawson, N. Gabor and V. Aji, *Physical Review B* **113**, 235306 (2026).

![Temperature dependence of the absorption spectra](docs/temperature_spectra.png)

## Models and units

Model I describes intraband exciton–phonon coupling. For independent exciton
channels, its spectrum has the form

$$
\alpha(E)=\sum_a A_a\sum_m W_{a,m}(T)
\frac{\Gamma_a}{(E-E_a-m\Omega_a)^2+\Gamma_a^2}.
$$

Here $E_a$ is the exciton energy, $\Omega_a$ the phonon energy, $\Gamma_a$ the
Lorentzian half width at half maximum, and $A_a$ the squared optical amplitude.
The Lorentzian convention does not include a factor of $1/\pi$.
The thermal sideband weights reduce to a Poisson distribution at zero
temperature and to a single zero-phonon term at zero coupling.

Model II includes coherent interband contributions. Its coupling parameters
are $g_{\mathrm{inter}}=g\sin\theta$ and
$g_{\mathrm{diag}}=g\cos\theta$. The intraband limit is
$\theta=0$ with Model I coupling $g_{cv}=2g$.

All energies, linewidths and `kBT_eV` are in eV. Couplings and mixing angles
are dimensionless, with angles in radians. Different phonon energies and
linewidths can be assigned to independent exciton channels; each channel
contains one phonon mode.

## Installation and first calculation

Python 3.10 or later is required. From the repository root:

```sh
python -m pip install .
python -m examples.theory_spectra
```

For a spectrum evaluated directly in Python:

```python
import numpy as np
from phonon_absorption import spectrum
from phonon_absorption.parameters import ANCHORS

energy = np.linspace(0.85, 1.45, 512)
absorption = spectrum(energy, model="II", cutoff=None, **ANCHORS["II"])
```

Figure export uses LaTeX. Install `latex`, `dvipng`, `type1cm`, `type1ec`,
`amsmath` and `textcomp`; on Ubuntu the packages are
`texlive-latex-extra cm-super dvipng ghostscript`.

## Fitting spectra

The CSV input contains `energy_eV` and `absorption`, with an optional `sigma`
column for standard deviations. Exciton energies are fixed in the supplied
configurations. Phonon energies, linewidths, temperature, couplings and mixing
angles are varied within the configured bounds. Nonnegative channel amplitudes
are solved linearly at each nonlinear trial.

To generate a seeded synthetic example and fit it:

```sh
python -m phonon_absorption.cli demo
python -m phonon_absorption.cli fit \
  --data examples/output/synthetic_observations_II.csv \
  --config configs/model_ii_shared.json \
  --output examples/output/my_fit
```

For experimental observations, replace `--data` with the local CSV path.
The output contains fitted parameters, predictions, residuals and fit figures
in PDF, PNG and SVG format.

| Configuration | Phonon and linewidth assumptions |
| --- | --- |
| `configs/model_i_shared.json` | Model I; shared phonon energy and linewidth |
| `configs/model_ii_shared.json` | Model II; shared phonon energy and linewidth |
| `configs/model_ii_independent_frequencies.json` | Model II; independent phonon energies, shared linewidth |
| `configs/model_ii_independent_fixed_modes.json` | Model II; fixed phonon energies 0.033 and 0.0295 eV, independent linewidths |

The public fitting examples use synthetic observations with Gaussian noise
$\sigma=0.02$. Experimental comparison-point arrays are not distributed.
Their literature source and the local-fit summaries are described in
[data provenance](docs/provenance.md) and
[the aggregate fit report](validation/local_refit_summary.json).

## Figures and numerical checks

- [Rounded-parameter spectra](docs/theory_spectra.png)
- [Temperature dependence](docs/temperature_spectra.png)
- [Exciton-channel decomposition](docs/channel_decomposition.png)
- [Model I synthetic fit](docs/synthetic_fit_I.png) and [Model II synthetic fit](docs/synthetic_fit_II.png)

```sh
python -m unittest discover -s tests -v
python -m validation.run
```

A fixed `cutoff=80` sums sidebands from $m=-80$ to $80$.
`cutoff=None` selects adaptive support and reports the omitted thermal weight.
Truncated weights are not renormalized. The checks cover the zero-temperature,
zero-coupling and intraband limits, agreement with the source expressions,
sideband convergence, and synthetic-spectrum recovery. Parameter bounds,
fit diagnostics and truncation estimates are documented in
[numerical methods](docs/numerics.md); measured results are in
[validation/report.json](validation/report.json).

In the residual panels, $\Delta\alpha$ denotes observations minus the fitted spectrum.

The figures use the shared [publication style](docs/figure_style.md).

## Reference

Yifan Liu, Robert Dawson, Nathaniel Gabor and Vivek Aji, “Phonon-assisted
absorption in transition metal dichalcogenide heterostructures,”
*Phys. Rev. B* **113**, 235306 (2026).
[DOI:10.1103/v3sr-lg47](https://doi.org/10.1103/v3sr-lg47) ·
[arXiv:2602.02803](https://arxiv.org/abs/2602.02803).
BibTeX and citation metadata are provided in `citation.bib` and `CITATION.cff`.
