import pandas as pd
import numpy as np
import periodictable as pt
from .log import Handle

logger = Handle(__name__)


def _fltr_by_n(low=None, high=None):
    _fltr = [lambda x: x]
    if low is not None:
        _fltr += [lambda x: x.number >= low]
    if high is not None:
        _fltr += [lambda x: x.number <= high]

    return [el for el in pt.elements if all([f(el) for f in _fltr])]


def get_periodic_frame():
    """
    Construct a simple periodic table dataframe organised by group and row. Note that
    the lanthanides and actinides are each found in a single cell.

    Returns
    -------
    :class:`pandas.DataFrame`
    """

    table = pd.DataFrame(columns=np.arange(1, 19), index=np.arange(1, 8))

    table.loc[1, [1, 2]] = [pt.H, pt.He]
    table.loc[2, [1, 2]] = [pt.Li, pt.Be]
    table.loc[2, np.arange(13, 19)] = _fltr_by_n(5, 10)
    table.loc[3, [1, 2]] = [pt.Na, pt.Mg]
    table.loc[3, np.arange(13, 19)] = _fltr_by_n(13, 18)
    table.loc[4, np.arange(1, 19)] = _fltr_by_n(19, 36)
    table.loc[5, np.arange(1, 19)] = _fltr_by_n(37, 54)
    table.loc[6, [1, 2]] = [pt.Cs, pt.Ba]
    table.loc[6, 3] = _fltr_by_n(57, 71)
    table.loc[6, np.arange(4, 19)] = _fltr_by_n(72, 86)
    table.loc[7, [1, 2]] = [pt.Fr, pt.Ra]
    table.loc[7, 3] = _fltr_by_n(89, 103)
    table.loc[7, np.arange(4, 19)] = _fltr_by_n(104, 118)
    return table
