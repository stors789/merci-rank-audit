"""Shared journal-sized figure styling for MERCI Rank Audit."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

INK = '#243540'
TEAL = '#147D80'
TERRA = '#D57A5E'
SLATE = '#7B8EA2'
GRID = '#E7ECEF'
MUTED = '#657780'
BLUE = '#3D6E88'
OCHRE = '#B68B52'

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 8.1,
    'text.color': INK, 'axes.labelcolor': INK, 'axes.edgecolor': INK,
    'xtick.color': INK, 'ytick.color': INK, 'axes.labelsize': 8.2,
    'xtick.labelsize': 7.7, 'ytick.labelsize': 7.7, 'legend.fontsize': 7.5,
    'axes.titlesize': 8.6, 'axes.titleweight': 'medium', 'axes.linewidth': .6,
    'xtick.major.width': .6, 'ytick.major.width': .6,
    'xtick.major.size': 3, 'ytick.major.size': 3,
    'axes.spines.top': False, 'axes.spines.right': False,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'figure.facecolor': 'white', 'savefig.facecolor': 'white',
})


def panel(ax, label, title):
    ax.text(0, 1.10, label, transform=ax.transAxes, fontsize=9.6,
            fontweight='bold', color=INK, va='bottom')
    ax.text(.105, 1.10, title, transform=ax.transAxes, fontsize=8.5,
            fontweight='semibold', color=INK, va='bottom')


def grid(ax, axis='y'):
    ax.set_axisbelow(True)
    ax.grid(axis=axis, color=GRID, linewidth=.65)
    ax.tick_params(which='minor', length=0)


def save(fig, directory, name):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    for suffix in ('pdf', 'svg', 'png'):
        path = directory / f'{name}.{suffix}'
        fig.savefig(path, dpi=350)
        if suffix == 'svg':
            lines = path.read_text().splitlines()
            path.write_text('\n'.join(line.rstrip(' \t') for line in lines) + '\n')
    plt.close(fig)
