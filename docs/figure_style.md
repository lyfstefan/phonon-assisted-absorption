# Figure style v1: compact square frames, LaTeX

Following the author's updated instructions, all ordinary letters, numbers,
legends and mathematical expressions use real LaTeX (`text.usetex=True`).
The exported PDF embeds Computer Modern Roman/Math Italic fonts (CMR/CMMI),
matching the requested LaTeX appearance. This explicitly supersedes v0's
`text.usetex=False` and LM Roman font-file pinning.

Role-based sizes and the original visual constants are preserved: axis labels
10 pt; ticks, legends and annotations 9 pt; frames 1 pt; PNG exports 300 dpi.
Simple figures are 7.6 cm wide and multipanel figures 15.2 cm wide.
Every data axes frame uses `Axes.set_box_aspect(1)`: its physical width and
height are equal. This does not require a square overall canvas, nor equal
x/y data units. Canvas height follows the panel arrangement; colorbars
remain slim. PDF pages are not cropped at export. SVG versions
and exact typography/size records accompany every figure.

Numerical arrays, fitted parameters, probability normalization and continuous
color limits are unchanged. The quantum-walk overview is arranged as three
rows and four columns; its detailed two-case plots retain eigenvalue legends.

## LaTeX dependencies

Install a TeX distribution with `latex`, `dvipng`, `type1cm`, `type1ec` and
AMS packages. On Ubuntu, the CI uses:

```sh
sudo apt-get install texlive-latex-extra cm-super dvipng ghostscript
```

The legacy `--font-path` option remains accepted for compatibility but does
not override TeX's Computer Modern selection.
