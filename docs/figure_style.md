# Figure style v0

Figures use the user-supplied `publication_style.py` template without changes
to its style functions. Axis labels are 10 pt; ticks, legends and annotations
are 9 pt; frames are 1 pt. PNG exports are 300 dpi. Standalone PDF pages are
8.6 cm wide for simple single-panel figures and 17.8 cm for multipanel figures,
with no bounding-box crop. SVG versions remain available.

The local release explicitly pins `lmroman10-regular.otf` (LM Roman 10 Regular);
PDF font resources were checked for `LMRoman10-Regular`. Math uses Matplotlib's
Computer Modern (`mathtext.fontset = cm`). Numerical arrays, fitted parameters,
probability normalization and heat-map color limits are unchanged. Captions
are kept in the project documentation; long panels have been rearranged.

For portable regeneration, use `--font-path /path/to/lmroman10-regular.otf`.
Without a supplied font file, the template uses its installed-font fallback
list. Some LM OpenType files report the family name `Latin Modern Roman`, so
the literal `LM Roman 10` name alone does not always find them. Pinning the
file avoids this mismatch. The font file itself is not bundled. Each figure's
`*_typography.json` records the exact file used and printed page width.

```sh
python -m phonon_absorption.cli demo --font-path /path/to/lmroman10-regular.otf
python -m examples.theory_spectra --font-path /path/to/lmroman10-regular.otf
```

Public figures contain synthetic observations or theory only. Real-data
comparison figures were restyled locally from saved fits without refitting.
