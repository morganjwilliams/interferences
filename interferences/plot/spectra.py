"""
Visualisation of mass spectra (intensity vs. m/z).
"""
import matplotlib.pyplot as plt
import pyrolite.plot
import numpy as np
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
    intensity_threshold=0.00001,
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
    intensity_threshold : :class:`float`
        Threshold for low-intensity peaks to ignore.
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
        if window is not None:
            # filter the table to match the window
            table = table.loc[table["m_z"].between(*window), :]

        if not "label" in table.columns:
            logger.debug("Fetching labels.")
            table["label"] = get_molecule_labels(table)
            logger.debug("Labels fetched.")
    elif components is not None:
        table = build_table(
            components,
            window=window,
            add_labels=True,
            **subkwargs(kwargs, build_table),
        )
    else:
        raise AssertionError("Either components or a built table must be supplied.")

    logger.debug("Plotting {} peaks.".format(table.index.size))
    ax = table.loc[:, ["m_z", yvar]].pyroplot.stem(ax=ax, **kwargs)

    if window is not None:
        ax.set_xlim(window)

    ax.set_yscale("log")
    ymin = intensity_threshold / 10.0
    ax.set_ylim(ymin, 1000.0)

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
                xytext=(table.loc[row, "m_z"], table.loc[row, yvar]),
                fontsize=12,
                fontweight=weights.get(row.count("["), 200)
                # rotation=90,
            )
            annotations.append(_an)
            # keep x position

    if adjust_labels and _have_adjustText:
        # add an empty rectangle over the peaks
        if window is None:
            xm, dx = np.mean(ax.get_xlim()), np.diff(ax.get_xlim()) / 2
        else:
            xm, dx = np.mean(window), np.diff(window) / 2
        ymax = 1
        ym, dy = (ymin + ymax) / 2, (ymax - ymin) / 2
        patch = rect_from_centre(xm, ym, dx=dx, dy=dy, fill=False, alpha=0)
        ax.add_patch(patch)
        adjust_text(
            annotations,
            table["m_z"].values,
            table[yvar].values,
            force_text=(0.2, 0.6),
            force_objects=(0, 0.3),
            force_points=(0.5, 0.5),
            expand_text=(1.5, 1.5),
            expand_points=(1.5, 1.5),
            # only_move={"text": "y"},
            # autoalign=False,
            # va="bottom",
            # ha="center",
            add_objects=[patch],
            arrowprops=dict(
                arrowstyle="-",
                connectionstyle="arc3,rad=0.1",
                color="0.5",
                ls="--",
                fc="w",
            ),
        )

    return ax
