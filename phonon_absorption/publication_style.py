"""General publication style v0: portable, with no paper-specific physics."""
from pathlib import Path
import json

import matplotlib as mpl
from matplotlib import font_manager
from matplotlib.text import Text

C = dict(blue='#0072B2', red='#D55E00', green='#009E73', orange='#E69F00',
         purple='#CC79A7', sky='#56B4E9', yellow='#F0E442', gray='#666666',
         black='#1a1a1a', light='#d9d9d9')
DASH = ['-', (0,(4.,1.6)), (0,(1.2,1.3)), (0,(5.5,1.5,1.2,1.5)), (0,(2.5,1.2))]
_ROMAN_FILE = None


def pin_roman_face(fig):
    if _ROMAN_FILE is not None:
        for text in fig.findobj(match=Text):
            prop=text.get_fontproperties().copy();prop.set_file(_ROMAN_FILE);text.set_fontproperties(prop)

AXIS_FONT = 10
DETAIL_FONT = 9


def unify_typography(fig):
    """Enforce role-based sizes before layout/legend checks and at export."""
    for ax in fig.axes:
        ax.xaxis.label.set_fontsize(AXIS_FONT)
        ax.yaxis.label.set_fontsize(AXIS_FONT)
        ax.tick_params(which='both', labelsize=DETAIL_FONT)
        ax.xaxis.offsetText.set_fontsize(DETAIL_FONT)
        ax.yaxis.offsetText.set_fontsize(DETAIL_FONT)
        for text in [ax.title, ax._left_title, ax._right_title, *ax.texts]:
            text.set_fontsize(DETAIL_FONT)
        legend=ax.get_legend()
        if legend:
            for text in [*legend.get_texts(),legend.get_title()]:text.set_fontsize(DETAIL_FONT)
    for text in fig.texts:text.set_fontsize(AXIS_FONT)
    for legend in fig.legends:
        for text in [*legend.get_texts(),legend.get_title()]:text.set_fontsize(DETAIL_FONT)


def width(columns=1, *, single_cm=8.6, double_cm=17.8):
    if columns not in (1, 2):
        raise ValueError('Use one or two journal columns')
    return (single_cm if columns == 1 else double_cm)/2.54


def apply_style(*, font_path=None):
    """Apply defaults; optionally pin a portable font file supplied by the caller.

    No cluster paths or repository data are needed. If LM Roman 10 is not
    installed, Matplotlib uses CMU Serif or DejaVu Serif instead.
    """
    global _ROMAN_FILE
    _ROMAN_FILE=None
    family=['LM Roman 10','CMU Serif','DejaVu Serif']
    if font_path is not None:
        path=Path(font_path).expanduser().resolve()
        if not path.is_file():raise FileNotFoundError(path)
        font_manager.fontManager.addfont(str(path));_ROMAN_FILE=str(path)
        family.insert(0,font_manager.FontProperties(fname=str(path)).get_name())
    mpl.rcParams.update({'text.usetex':False,'font.family':'serif','font.serif':family,
        'mathtext.fontset':'cm','font.size':DETAIL_FONT,
        'xtick.direction':'in','ytick.direction':'in','figure.dpi':150})
    mpl.rcParams.update({'axes.linewidth': 1., 'lines.linewidth': 1.2,
        'xtick.top': False, 'ytick.right': False,
        'xtick.major.width': .9, 'ytick.major.width': .9,
        'xtick.major.size': 3., 'ytick.major.size': 3.,
        'axes.titlesize': DETAIL_FONT, 'axes.labelsize': AXIS_FONT, 'legend.fontsize': DETAIL_FONT,
        'xtick.labelsize': DETAIL_FONT, 'ytick.labelsize': DETAIL_FONT,
        'savefig.bbox': None, 'pdf.fonttype': 42, 'ps.fonttype': 42})


def save(fig, destination):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    unify_typography(fig)
    # Apply at export too: individual legacy plotters sometimes set heavier
    # frames explicitly. Keep all data artists, limits and normalizations.
    for ax in fig.axes:
        for spine in ax.spines.values():
            spine.set_linewidth(1.)
        ax.tick_params(which='both', top=False, right=False, width=.9)
        ax.tick_params(which='major', length=3.)
        if hasattr(ax, '_colorbar'):
            ax._colorbar.outline.set_linewidth(.9)
    pin_roman_face(fig)
    # Do not crop the page: PDF width must equal the requested printed width.
    # Decimal angles in a stem (theta1.05) are not filename extensions.
    fig.savefig(destination.parent/(destination.name+'.pdf'))
    fig.savefig(destination.parent/(destination.name+'.png'), dpi=300)
    (destination.parent/(destination.name+'_typography.json')).write_text(json.dumps({
        'template':'v0','frame_pt':1.,
        'axis_labels_pt':AXIS_FONT,'ticks_pt':DETAIL_FONT,'colorbar_ticks_pt':DETAIL_FONT,
        'annotations_pt':DETAIL_FONT,'legends_pt':DETAIL_FONT,
        'page_width_cm':float(fig.get_figwidth()*2.54),
        'font_family':list(mpl.rcParams['font.serif']), 'math_fontset':'cm',
        'pinned_font_file':_ROMAN_FILE},indent=2)+'\n')
