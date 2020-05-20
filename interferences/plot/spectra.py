"""
Visualisation of mass spectra (intensity vs. m/z).
"""
import matplotlib.pyplot as plt
import pyrolite.plot

from adjustText import adjust_text
from ..util.mz import process_window
from ..table import build_table
from ..table.molecules import get_molecule_labels
import logging


logging.getLogger(__name__).addHandler(logging.NullHandler())

try:
    from pyrolite.util.meta import subkwargs

    _have_adjustText = True
except ImportError:
    _have_adjustText = False


def stemplot(
    components=None,
    table=None,
    window=None,
    threshold=10e-8,
    yvar="iso_product",
    adjust_labels=True,
    **kwargs,
):
    """
    Create a stemplot based spectra (zero-width peaks)
    """
    window = process_window(window)
    if table is not None:
        if not "label" in table.columns:
            table["label"] = get_molecule_labels(table)
        if window is not None:
            # filter the table to match the window
            table = table.loc[table.m_z >= window[0] & table.m_z <= window[1], :]
    elif components is not None:
        table = build_table(
            components, window=window, add_labels=True, **subkwargs(kwargs, build_table)
        )
    else:
        raise AssertionError("Either components or a built table must be supplied.")

    ax = table.loc[:, ["m_z", yvar]].pyroplot.stem(**kwargs)
    ax.set_xlim(window)
    ax.set_yscale("log")
    ax.set_ylim(threshold / 10, 10)

    annotations = []
    for row in table.index:
        abund = table.loc[row, yvar]
        if abund > threshold:
            # if abund < 0.5 :
            _an = ax.annotate(
                table.loc[row, "label"],
                xy=(table.loc[row, "m_z"], table.loc[row, yvar]),
                fontsize=12,
                # rotation=90,
            )
            annotations.append(_an)
            # keep x position

    if adjust_labels and _have_adjustText:
        adjust_text(
            annotations,
            table["m_z"].values,
            table[yvar].values,
            force_text=(0.1, 0.8),
            # force_points=(0.2, 0.5),
            # expand_points=(1.05, 1.05),
            # only_move={"text": "y"},
            autoalign=None,
            va="bottom",
            ha="center",
            arrowprops=dict(
                arrowstyle="-",
                connectionstyle="arc3,rad=0.2",
                color="0.5",
                ls="--",
                fc="w",
            ),
        )

    return ax
