"""
Visualisation of mass spectra (intensity vs. m/z).
"""
import matplotlib.pyplot as plt
import pyrolite.plot
import numpy as np
from adjustText import adjust_text
from collections import defaultdict
import matplotlib.lines
from ..util.mz import process_window
from ..table import build_table
from ..table.molecules import get_molecule_labels
from .kernel import peak, peak_kernel
from pyrolite.util.plot.helpers import rect_from_centre
from pyrolite.util.plot.axes import init_axes
from ..util.log import Handle

logger = Handle(__name__)

try:
    from pyrolite.util.meta import subkwargs

    _have_adjustText = True
except ImportError:
    _have_adjustText = False


def _get_table(components=None, table=None, window=None, **kwargs):
    """
    Parameters
    ----------
    components : :class:`list`
        List of components to include in the table.
    table : :class:`pandas.DataFrame`
        Table of interferences to use for the plot.
    window : :class:`tuple`
        Window in m/z to use for the plot. Can specify (low, high) or (isotope, width).
    """
    if table is not None:
        if window is not None:
            # filter the table to match the window
            table = table.loc[table["m_z"].between(*window), :]
        if not "label" in table.columns:
            logger.debug("Fetching labels.")
            table.loc[:, "label"] = get_molecule_labels(table)
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
    return table


def _label_peaks(
    ax,
    table,
    yvar="iso_product",
    window=None,
    intensity_threshold=0.00001,
    adjust_labels=True,
    add_patch=True,
    max_labels=12,
    ymin=0.00001,
    ymax=1,
    iter_lim=100,
):
    """
    Parameters
    ----------
    ax : :class:`matplotlib.axes.Axes`
        Axes to add the labels to.
    table : :class:`pandas.DataFrame`
        Table of interferences to use for the plot.
    intensity_threshold : :class:`float`
        Threshold for low-intensity peaks to ignore for labelling.
    adjust_labels : :class:`bool`
        Whether to adjust the label positions for clarity.
    add_patch : :class:`bool`
        Whether to add a label-deflecting patch over the peak area.
    max_labels : :class:`int`
        Maximum labels to add to the plot.
    """

    annotations = []
    # if it's a primary peak (i.e. one elmeent), make it bold
    weights = defaultdict(lambda: "light")
    weights.update({ix + 1: weight for ix, weight in enumerate(["black", "normal"])})
    for row in table.nlargest(max_labels, yvar).index:
        intensity = table.loc[row, yvar]
        if intensity > intensity_threshold:
            _an = ax.annotate(
                table.loc[row, "label"],
                xy=(table.loc[row, "m_z"], table.loc[row, yvar]),
                xytext=(table.loc[row, "m_z"], table.loc[row, yvar]),
                fontsize=12,
                fontweight=weights[row.count("[")]
                # rotation=90,
            )
            annotations.append(_an)

    if adjust_labels and _have_adjustText:
        add_objs = []
        if add_patch:
            # add an empty rectangle over the peaks
            if window is None:
                xm, dx = np.mean(ax.get_xlim()), np.diff(ax.get_xlim()) / 2
            else:
                xm, dx = np.mean(window), np.diff(window) / 2

            ym, dy = (ymin + ymax) / 2, (ymax - ymin) / 2
            patch = rect_from_centre(xm, ym, dx=dx, dy=dy, fill=False, alpha=0)
            ax.add_patch(patch)
            add_objs.append(patch)
        adjust_text(
            annotations,
            table["m_z"].values,
            table[yvar].values,
            ax=ax,
            force_text=(0.2, 1),
            force_objects=(0, 0.5),
            force_points=(0, 0.5),
            expand_text=(1.5, 1.5),
            expand_points=(4, 4),
            # only_move={"text": "y"},
            # autoalign=False,
            # va="bottom",
            # ha="center",
            add_objects=add_objs,
            arrowprops=dict(
                arrowstyle="-",
                connectionstyle="arc3,rad=0.1",
                color="0.5",
                ls="--",
                fc="w",
            ),
            lim=iter_lim,
            on_basemap=True,
        )


def _format_axes(ax, window=None, ymin=0.00001):
    """
    Parameters
    ----------
    ax : :class:`matplotlib.axes.Axes`
        Axes to format for spectra.
    window : :class:`tuple`
        Window in m/z to use for the plot. Can specify (low, high) or (isotope, width).
    ymin : :class:`float`
        Minimum value for the y-axis.
    """
    ax.set_ylabel("Estimated Relative Intensity")

    if window is not None:
        ax.set_xlim(window)

    ax.set_yscale("log")
    ax.set_ylim(ymin, 1000.0)


def stemplot(
    ax=None,
    components=None,
    table=None,
    window=None,
    yvar="iso_product",
    ymin=0.00001,
    **kwargs
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
    yvar : :class:`str`
        Column to use for the y-axis.
    ymin : :class:`float`
        Minimum value for the y-axis.

    Returns
    -------
    :class:`matplotlib.axes.Axes`
    """

    window = process_window(window)
    table = _get_table(components=components, table=table, window=window, **kwargs)
    logger.debug("Plotting {} peaks.".format(table.index.size))
    ax = table.loc[:, ["m_z", yvar]].pyroplot.stem(ax=ax, **kwargs)

    _format_axes(ax, window=window, ymin=ymin)
    _label_peaks(
        ax,
        table,
        yvar=yvar,
        window=window,
        ymin=ymin,
        **subkwargs(kwargs, _label_peaks),
    )

    return ax


def spectra(
    ax=None,
    components=None,
    table=None,
    window=None,
    mass_resolution=1000,
    image_ratio=0.0,
    yvar="iso_product",
    ymin=0.00001,
    **kwargs
):
    """
    Create a spectrum based on peaks from a table.

    Parameters
    ----------
    ax : :class:`matplotlib.axes.Axes`
        Axes to add the plot to, optional.
    components : :class:`list`
        List of components to include in the table.
    table : :class:`pandas.DataFrame`
        Table of interferences to use for the plot.
    mass_resolution : :class:`float`
        Mass resolution :math:`\Delta M/M` defined at Full Width Half Maximum (FWHM),
        used to scale the width/mass range of the peak.
    image_ratio : :class:`float`
        The ratio of the size of the image to the limiting slit. Values between 0
        (zero width image) and 1 correspond to scenarios where the full image fits
        in a collector slit. Beyond 1, this corresponds to scenarios where the image
        is larger than the slit, with reduced maximum intensities.
    window : :class:`tuple`
        Window in m/z to use for the plot. Can specify (low, high) or (isotope, width).
    yvar : :class:`str`
        Column to use for the y-axis.
    ymin : :class:`float`
        Minimum value for the y-axis.

    Returns
    -------
    :class:`matplotlib.axes.Axes`
    """
    window = process_window(window)
    table = _get_table(components=components, table=table, window=window, **kwargs)
    logger.debug("Plotting {} peaks.".format(table.index.size))

    ax = init_axes(ax=ax)
    # get kernel
    krnl = peak_kernel(
        mass_resolution=mass_resolution,
        image_ratio=image_ratio,
        **subkwargs(kwargs, peak_kernel),
    )
    for (name, p) in table.iterrows():
        idx, signal = peak(p["m_z"], p[yvar], kernel=krnl)
        ax.plot(
            idx,
            signal,
            **subkwargs(kwargs, ax.plot, ax.scatter, matplotlib.lines.Line2D),
        ),

    _format_axes(ax, window=window, ymin=ymin)
    _label_peaks(
        ax,
        table,
        yvar=yvar,
        window=window,
        ymin=ymin,
        **subkwargs(kwargs, _label_peaks),
    )

    return ax
