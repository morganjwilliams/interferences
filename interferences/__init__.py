import sys
import numpy as np
import periodictable as pt
import pandas as pd
import logging
from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions

logging.getLogger(__name__).addHandler(logging.NullHandler())
logging.captureWarnings(True)

from .table.build import build_table
from .plot.spectra import stemplot
from .util.mz import process_window

# note that only some of these methods will be valid for series
@pd.api.extensions.register_series_accessor("mz")
@pd.api.extensions.register_dataframe_accessor("mz")
class mz(object):
    def __init__(self, obj):
        """
        Custom dataframe accessor for interferences.
        """
        self._validate(obj)
        self._obj = obj

    @staticmethod
    def _validate(obj):
        pass

    def get_window(self, window):
        """
        Get a m/z window from a table.

        Parameters
        ----------
        window : :class:`tuple`
            Window specification. Either low and high m/z, or isotope and width.

        Returns
        -------
        :class:`pandas.DataFrame`
            Filtered dataframe.
        """
        return self._obj.loc[self._obj.m_z.between(*process_window(window))]

    def stemplot(self, *args, **kwargs):
        """
        Get a m/z window from a table.

        Parameters
        ----------
        window : :class:`tuple`
            Window specification. Either low and high m/z, or isotope and width.

        Returns
        -------
        :class:`pandas.DataFrame`
            Filtered dataframe.
        """
        return stemplot(table=self._obj, *args, **kwargs)
