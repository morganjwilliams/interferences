"""
Functions for working with tables of molecular ions.
"""
import pandas as pd
import numpy as np
import periodictable as pt
from .build import build_table
import logging

logging.getLogger(__name__).addHandler(logging.NullHandler())
logger = logging.getLogger(__name__)

# note that only some of these methods will be valid for series
@pd.api.extensions.register_series_accessor("interf")
@pd.api.extensions.register_dataframe_accessor("interf")
class interf(object):
    def __init__(self, obj):
        """
        Custom dataframe accessor for pyrolite compositional transforms.
        """
        self._validate(obj)
        self._obj = obj

    @staticmethod
    def _validate(obj):
        pass

    def get_window(self, low, high):
        """
        Get a m/z window from a table.

        Parameters
        ----------
        components : :class:`list`
            Option subcompositon to renormalise to 100. Useful for the use case
            where compostional data and non-compositional data are stored in the
            same dataframe.
        scale : :class:`float`, :code:`100.`
            Closure parameter. Typically either 100 or 1.

        Returns
        -------
        :class:`pandas.DataFrame`
            Renormalized dataframe.

        Notes
        ------
        This won't modify the dataframe in place, you'll need to assign it to something.
        If you specify components, those components will be summed to 100%,
        and others remain unchanged.
        """
        obj = self._obj
        return renormalise(obj, components=components, scale=scale)
