"""
Visualisation of mass spectra (intensity vs. m/z).
"""
import matplotlib.pyplot as plt
import pyrolite.plot

from adjustText import adjust_text
from ..util.mz import process_window
from ..table import build_table
from ..table.molecules import get_molecule_labels
from pyrolite.util.plot.helpers import rect_from_centre
from ..util.log import Handle

logger = Handle(__name__)

try:
    from pyrolite.util.meta import subkwargs

    _have_adjustText = True
except ImportError:
    _have_adjustText = False


def stemplot(
    ax=None,
    components=None,
    table=None,
    window=None,
    threshold=10e-8,
    yvar="iso_product",
    adjust_labels=True,
    **kwargs,
):
    """
    Create a stemplot based spectra (zero-width peaks).

    Parameters
    ----------
    ax : :class:`matplotlib.axes.Axes`
        Axes to add the plot to, optional.
    components : :class:`list`
        List of components to include in the table.
    table : :class:`pandas.DataFrame`
        Table of interferences to use for the plot.
    window : :class:`tuple`
        Window in m/z to use for the plot. Can specify (low, high) or (isotope, width).
    threshold : :class:`float`
        Threshold for low-abundance isotopes to ignore.
    yvar : :class:`str`
        Column to use for the y-axis.
    adjust_labels : :class:`bool`
        Whether to adjust the label positions for clarity.

    Returns
    -------
    :class:`matplotlib.axes.Axes`
    """
    window = process_window(window)
    if table is not None:
        if not "label" in table.columns:
            table["label"] = get_molecule_labels(table)
        if window is not None:
            # filter the table to match the window
            table = table.loc[table["m_z"].between(*window), :]
    elif components is not None:
        table = build_table(
            components, window=window, add_labels=True, **subkwargs(kwargs, build_table)
        )
    else:
        raise AssertionError("Either components or a built table must be supplied.")

    ax = table.loc[:, ["m_z", yvar]].pyroplot.stem(ax=ax, **kwargs)

    if window is not None:
        ax.set_xlim(window)

    ax.set_yscale("log")
    ax.set_ylim(threshold / 10.0, 1000.0)

    annotations = []
    for row in table.index:
        abund = table.loc[row, yvar]
        if abund > threshold:
            # if abund < 0.5 :
            # if it's a primary peak (i.e. one elmeent), make it bold
            weights = {ix: weight for ix, weight in enumerate([800, 600, 400])}
            _an = ax.annotate(
                table.loc[row, "label"],
                xy=(table.loc[row, "m_z"], table.loc[row, yvar]),
                fontsize=12,
                fontweight=weights.get(row.count("["), 200)
                # rotation=90,
            )
            annotations.append(_an)
            # keep x position

    if adjust_labels and _have_adjustText:
        # add an empty rectangle over the peaks
        # xm, xw = np.mean(ax.get_xlim()), np.diff(ax.get_xlim())/2
        # patch = rect_from_centre(xm, ym,width fill=False, alpha=0)

        adjust_text(
            annotations,
            table["m_z"].values,
            table[yvar].values,
            force_text=(0.05, 0.6),
            # force_points=(0.2, 0.5),
            # expand_points=(1.05, 1.05),
            # only_move={"text": "y"},
            autoalign=None,
            va="bottom",
            ha="center",
            # add_objects=[patch],
            arrowprops=dict(
                arrowstyle="-",
                connectionstyle="arc3,rad=0.2",
                color="0.5",
                ls="--",
                fc="w",
            ),
        )

    return ax
