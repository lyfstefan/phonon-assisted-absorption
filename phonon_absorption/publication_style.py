"""Publication style v1: compact figures with square data frames and real LaTeX text.

Retains the user's v0 colors, line styles, frames and role-based type sizes.
User-requested changes: TeX, reference-matched 2.5 cm square frames, and
independent short colorbars. Canvas dimensions follow the panel arrangement.
"""

from pathlib import Path
import json
import shutil
import os
import tempfile
from io import BytesIO
from statistics import median
import math
from matplotlib.ticker import MaxNLocator, FuncFormatter
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
FRAME_CM = 2.5


def all_axes(fig):
    """Include per-panel inset colorbar axes in typography and export checks."""
    pending = list(fig.axes)
    seen = set()
    while pending:
        ax = pending.pop(0)
        if id(ax) in seen:
            continue
        seen.add(id(ax))
        yield ax
        pending.extend(ax.child_axes)


def panel_colorbar(fig, ax, mappable, label):
    """One short bar per data frame, bottom-aligned, with its label above."""
    cax = ax.inset_axes([1.045, 0, 0.045, 0.60])
    cb = fig.colorbar(mappable, cax=cax)
    cb.locator = MaxNLocator(nbins=3, steps=[1, 2, 5, 10])
    maximum = max(abs(float(mappable.norm.vmin)), abs(float(mappable.norm.vmax)))
    exponent = math.floor(math.log10(maximum)) if maximum > 0 else 0
    exponent = exponent if exponent <= -2 or exponent >= 2 else 0
    scale = 10.0**exponent
    cb.formatter = FuncFormatter(lambda value, pos: rf"${value/scale:.3g}$")
    cb.update_ticks()
    title = label + (chr(10) + rf"$\times10^{{{exponent}}}$" if exponent else "")
    cax.set_title(title, loc="left", pad=4, fontsize=DETAIL_FONT)
    cax._publication_parent = ax
    cax._publication_label = label
    cax._publication_exponent = exponent
    cax.tick_params(labelsize=DETAIL_FONT, pad=2)
    return cb


def fit_frame_size(fig, target_cm=FRAME_CM):
    """Match the reference's font/frame ratio while retaining the 10/9 pt sizes."""
    data = [ax for ax in all_axes(fig) if not hasattr(ax, "_colorbar")]
    if len(data) == 1:
        # Preserve physical space for external legends and headers. Repeated
        # constrained-layout scaling can collapse a single compact data panel.
        fig.canvas.draw()
        ax = data[0]
        box = ax.get_position()
        width_cm, height_cm = fig.get_size_inches() * 2.54
        left, right = box.x0 * width_cm, (1 - box.x1) * width_cm
        bottom, top = box.y0 * height_cm, (1 - box.y1) * height_cm
        new_width, new_height = left + target_cm + right, bottom + target_cm + top
        renderer = fig.canvas.get_renderer()
        header_width = (
            max(
                [
                    text.get_window_extent(renderer).width / fig.dpi * 2.54
                    for text in fig.texts
                ]
                + [0]
            )
            + 0.24
        )
        extra = max(0, header_width - new_width)
        left += extra / 2
        new_width += extra
        fig.set_layout_engine(None)
        fig.set_size_inches(new_width / 2.54, new_height / 2.54, forward=False)
        ax.set_position(
            [
                left / new_width,
                bottom / new_height,
                target_cm / new_width,
                target_cm / new_height,
            ]
        )
        fig.canvas.draw()
        return
    for _ in range(9):
        fig.canvas.draw()
        data = [ax for ax in all_axes(fig) if not hasattr(ax, "_colorbar")]
        side = median(
            ax.get_position().width * fig.get_figwidth() * 2.54 for ax in data
        )
        if abs(side - target_cm) < 0.035:
            break
        factor = max(0.65, min(1.6, target_cm / side))
        fig.set_size_inches(
            fig.get_figwidth() * factor, fig.get_figheight() * factor, forward=False
        )
    fig.canvas.draw()


def pad_figure(fig):
    """Freeze the finished layout and extend only paper margins for all text."""
    fig.canvas.draw()
    fig.set_layout_engine(None)
    for _ in range(2):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        tight = fig.get_tightbbox(renderer)
        old_w, old_h = fig.get_size_inches()
        left = max(0.0, 0.04 - tight.x0)
        right = max(0.0, tight.x1 + 0.04 - old_w)
        bottom = max(0.0, 0.04 - tight.y0)
        top = max(0.0, tight.y1 + 0.04 - old_h)
        if max(left, right, bottom, top) < 0.001:
            break
        positions = [(ax, ax.get_position().frozen()) for ax in fig.axes]
        texts = [(text, text.get_position()) for text in fig.texts]
        legends = [
            (
                legend,
                fig.transFigure.inverted()
                .transform_bbox(legend.get_bbox_to_anchor())
                .frozen(),
            )
            for legend in fig.legends
        ]
        new_w, new_h = old_w + left + right, old_h + bottom + top
        fig.set_size_inches(new_w, new_h, forward=False)
        for ax, box in positions:
            ax.set_position(
                [
                    (box.x0 * old_w + left) / new_w,
                    (box.y0 * old_h + bottom) / new_h,
                    box.width * old_w / new_w,
                    box.height * old_h / new_h,
                ]
            )
        for text, (x, y) in texts:
            text.set_position(
                ((x * old_w + left) / new_w, (y * old_h + bottom) / new_h)
            )
        for legend, box in legends:
            legend.set_bbox_to_anchor(
                (
                    (box.x0 * old_w + left) / new_w,
                    (box.y0 * old_h + bottom) / new_h,
                    box.width * old_w / new_w,
                    box.height * old_h / new_h,
                ),
                transform=fig.transFigure,
            )
    fig.canvas.draw()


def pin_roman_face(fig):
    """All labels use TeX's Computer Modern; no external Roman font is pinned."""
    for text in fig.findobj(match=Text):
        text.set_usetex(True)


def unify_typography(fig):
    for ax in all_axes(fig):
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


def width(columns=1, *, single_cm=6.8, double_cm=13.6):
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
            "text.latex.preamble": r"\usepackage{amsmath}\usepackage{textcomp}",
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


def save(fig, destination, *, legend_fontsize=None):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    unify_typography(fig)
    if legend_fontsize is not None:
        legends = [ax.get_legend() for ax in all_axes(fig)] + list(fig.legends)
        for legend in legends:
            if legend is not None:
                for text in [*legend.get_texts(), legend.get_title()]:
                    text.set_fontsize(legend_fontsize)
    for ax in all_axes(fig):
        if not hasattr(ax, "_colorbar"):
            ax.set_box_aspect(1)
        for spine in ax.spines.values():
            spine.set_linewidth(1.0)
        ax.tick_params(which="both", top=False, right=False, width=0.9)
        ax.tick_params(which="major", length=3.0)
        if hasattr(ax, "_colorbar"):
            ax._colorbar.outline.set_linewidth(0.9)
    pin_roman_face(fig)
    fit_frame_size(fig)
    pad_figure(fig)
    # Render fully before atomically replacing an existing export.
    # This avoids partial PDFs when a synced or open destination cannot be streamed to.
    frames = []
    bars = []
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
                for ax in all_axes(fig)
                if not hasattr(ax, "_colorbar")
            ]
        if extension == "pdf":
            bars = []
            for ax in all_axes(fig):
                if not hasattr(ax, "_publication_parent"):
                    continue
                parent = ax._publication_parent
                box = ax.get_position()
                frame = parent.get_position()
                bars.append(
                    {
                        "parent_bounds_figure_fraction": list(frame.extents),
                        "bounds_figure_fraction": list(box.extents),
                        "height_fraction_of_frame": float(box.height / frame.height),
                        "width_fraction_of_frame": float(box.width / frame.width),
                        "bottom_offset_cm": float(
                            (box.y0 - frame.y0) * fig.get_figheight() * 2.54
                        ),
                        "label": ax._publication_label,
                        "label_position": "above",
                        "shared": False,
                        "label_bottom_clearance_cm": float(
                            (
                                fig.transFigure.inverted()
                                .transform_bbox(
                                    ax._left_title.get_window_extent(
                                        fig.canvas.get_renderer()
                                    )
                                )
                                .y0
                                - box.y1
                            )
                            * fig.get_figheight()
                            * 2.54
                        ),
                        "scientific_exponent": ax._publication_exponent,
                    }
                )
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
                "legends_pt": (
                    DETAIL_FONT if legend_fontsize is None else legend_fontsize
                ),
                "page_width_cm": float(fig.get_figwidth() * 2.54),
                "page_height_cm": float(fig.get_figheight() * 2.54),
                "font_family": "Computer Modern Roman",
                "math_font_family": "Computer Modern Math Italic",
                "data_frames": frames,
                "colorbars": bars,
                "target_frame_side_cm": FRAME_CM,
                "latex_preamble": str(mpl.rcParams["text.latex.preamble"]),
                "pinned_font_file": None,
            },
            indent=2,
        )
        + "\n"
    )
