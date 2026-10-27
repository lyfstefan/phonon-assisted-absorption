"""Publication style v1: compact figures with square data frames and real LaTeX text.

Retains the user's v0 colors, line styles, frames and role-based type sizes.
User-requested changes: TeX typesetting and smaller 7.6/15.2 cm defaults.
"""

from pathlib import Path
import json
import shutil
import os
import tempfile
from io import BytesIO
import matplotlib as mpl
from matplotlib.text import Text

C = dict(
    blue="#0072B2",
    red="#D55E00",
    green="#009E73",
    orange="#E69F00",
    purple="#CC79A7",
    sky="#56B4E9",
    yellow="#F0E442",
    gray="#666666",
    black="#1a1a1a",
    light="#d9d9d9",
)
DASH = [
    "-",
    (0, (4.0, 1.6)),
    (0, (1.2, 1.3)),
    (0, (5.5, 1.5, 1.2, 1.5)),
    (0, (2.5, 1.2)),
]
AXIS_FONT = 10
DETAIL_FONT = 9


def pin_roman_face(fig):
    """All labels use TeX's Computer Modern; no external Roman font is pinned."""
    for text in fig.findobj(match=Text):
        text.set_usetex(True)


def unify_typography(fig):
    for ax in fig.axes:
        ax.xaxis.label.set_fontsize(AXIS_FONT)
        ax.yaxis.label.set_fontsize(AXIS_FONT)
        ax.tick_params(which="both", labelsize=DETAIL_FONT)
        ax.xaxis.offsetText.set_fontsize(DETAIL_FONT)
        ax.yaxis.offsetText.set_fontsize(DETAIL_FONT)
        for text in [ax.title, ax._left_title, ax._right_title, *ax.texts]:
            text.set_fontsize(DETAIL_FONT)
        legend = ax.get_legend()
        if legend:
            for text in [*legend.get_texts(), legend.get_title()]:
                text.set_fontsize(DETAIL_FONT)
    for text in fig.texts:
        text.set_fontsize(AXIS_FONT)
    for legend in fig.legends:
        for text in [*legend.get_texts(), legend.get_title()]:
            text.set_fontsize(DETAIL_FONT)


def width(columns=1, *, single_cm=7.6, double_cm=15.2):
    if columns not in (1, 2):
        raise ValueError("Use one or two journal columns")
    return (single_cm if columns == 1 else double_cm) / 2.54


def apply_style(*, font_path=None):
    """Use the installed LaTeX toolchain. font_path is a legacy no-op argument.

    TeX's CMR/CMMI fonts replace v0's Matplotlib font-file selection, following
    the explicit request that ordinary letters and mathematics all use LaTeX.
    """
    missing = [name for name in ["latex", "dvipng"] if shutil.which(name) is None]
    if missing:
        raise RuntimeError("LaTeX figure export requires: " + ", ".join(missing))
    mpl.rcParams.update(
        {
            "text.usetex": True,
            "font.family": "serif",
            "font.serif": ["Computer Modern Roman"],
            "text.latex.preamble": r"\usepackage{amsmath,amssymb}",
            "mathtext.fontset": "cm",
            "font.size": DETAIL_FONT,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "figure.dpi": 150,
            "axes.linewidth": 1.0,
            "lines.linewidth": 1.2,
            "xtick.top": False,
            "ytick.right": False,
            "xtick.major.width": 0.9,
            "ytick.major.width": 0.9,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            "axes.titlesize": DETAIL_FONT,
            "axes.labelsize": AXIS_FONT,
            "legend.fontsize": DETAIL_FONT,
            "xtick.labelsize": DETAIL_FONT,
            "ytick.labelsize": DETAIL_FONT,
            "savefig.bbox": None,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save(fig, destination):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    unify_typography(fig)
    for ax in fig.axes:
        if not hasattr(ax, "_colorbar"):
            ax.set_box_aspect(1)
        for spine in ax.spines.values():
            spine.set_linewidth(1.0)
        ax.tick_params(which="both", top=False, right=False, width=0.9)
        ax.tick_params(which="major", length=3.0)
        if hasattr(ax, "_colorbar"):
            ax._colorbar.outline.set_linewidth(0.9)
    pin_roman_face(fig)
    # Render fully before atomically replacing an existing export.
    # This avoids partial PDFs when a synced or open destination cannot be streamed to.
    frames = []
    for extension, kwargs in [("pdf", {}), ("png", {"dpi": 300}), ("svg", {})]:
        buffer = BytesIO()
        fig.savefig(buffer, format=extension, **kwargs)
        if extension == "pdf":
            frames = [
                {
                    "bounds_figure_fraction": list(ax.get_position().extents),
                    "width_cm": float(
                        ax.get_position().width * fig.get_figwidth() * 2.54
                    ),
                    "height_cm": float(
                        ax.get_position().height * fig.get_figheight() * 2.54
                    ),
                }
                for ax in fig.axes
                if not hasattr(ax, "_colorbar")
            ]
        target = destination.parent / (destination.name + "." + extension)
        with tempfile.NamedTemporaryFile(
            dir=destination.parent, prefix=".figure-", delete=False
        ) as stream:
            temporary = Path(stream.name)
            stream.write(buffer.getvalue())
        try:
            os.replace(temporary, target)
        finally:
            temporary.unlink(missing_ok=True)
    (destination.parent / (destination.name + "_typography.json")).write_text(
        json.dumps(
            {
                "template": "v1_compact_square_frames_latex",
                "text_engine": "latex",
                "text_usetex": True,
                "frame_pt": 1.0,
                "axis_labels_pt": AXIS_FONT,
                "ticks_pt": DETAIL_FONT,
                "colorbar_ticks_pt": DETAIL_FONT,
                "annotations_pt": DETAIL_FONT,
                "legends_pt": DETAIL_FONT,
                "page_width_cm": float(fig.get_figwidth() * 2.54),
                "page_height_cm": float(fig.get_figheight() * 2.54),
                "font_family": "Computer Modern Roman",
                "math_font_family": "Computer Modern Math Italic",
                "data_frames": frames,
                "pinned_font_file": None,
            },
            indent=2,
        )
        + "\n"
    )
