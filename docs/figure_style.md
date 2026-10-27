# Figure style v1: compact square frames, LaTeX

Following the author's updated instructions, all ordinary letters, numbers,
legends and mathematical expressions use real LaTeX (`text.usetex=True`).
The exported PDF embeds Computer Modern Roman/Math Italic fonts (CMR/CMMI),
matching the requested LaTeX appearance. This explicitly supersedes v0's
`text.usetex=False` and LM Roman font-file pinning.

Role-based sizes and the original visual constants are preserved: axis labels
10 pt; ticks, legends and annotations 9 pt; frames 1 pt; PNG exports 300 dpi.
The canvas is fitted to the reference image's font/frame proportions: each
data frame is approximately 2.5 cm square, with `Axes.set_box_aspect(1)`.
This does not set equal x/y data units. Canvas dimensions follow the number
and arrangement of panels, rather than a fixed overall paper width.

Each heatmap has its own thin colorbar: 60% of its frame height and 4.5%
of its frame width, with the bottoms aligned. A compact mathematical label
is placed above the bar. Comparable panels retain the same numerical
normalization even though their colorbars are separate. No long shared
colorbars are used. PDF pages are not cropped; PNG and SVG are also exported.

Numerical arrays, fitted parameters and probability normalization are unchanged.
Flux maps display the oriented phase modulo $2\pi$ on a cyclic scale. The quantum-walk overview is arranged as three
rows and four columns; its detailed two-case plots retain eigenvalue legends.

## LaTeX dependencies

Install a TeX distribution with `latex`, `dvipng`, `type1cm`, `type1ec` and
AMS and textcomp packages. On Ubuntu, the CI uses:

```sh
sudo apt-get install texlive-latex-extra cm-super dvipng ghostscript
```

The legacy `--font-path` option remains accepted for compatibility but does
not override TeX's Computer Modern selection.

The TeX preamble is `\usepackage{amsmath}\usepackage{textcomp}`.
Operator names use math subscripts (`$U_1$`, `$U_2$`), dimensionless ratios
use math mode (`$V/t$`), and phase letters use mathematical italics.
