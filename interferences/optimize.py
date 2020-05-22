"""
Functions for looking into ways for minimising interferences.
"""
import pandas as pd
import periodictable as pt
from .util.log import Handle

logger = Handle(__name__)


def sum_of_interferences(ion, composition=None, window=0.1):
    """
    Calculate the sum of all interferences for a given isotope.
    """
    result = 0.0

    if composition is None:
        if ion.isotope in _ELEMENTS[ion.element]:
            _interf = [v for (k, v) in _ELEMENTS[ion.element][ion.isotope].items()]
            if _interf:
                result += np.array(_interf).sum()
    else:
        assert isinstance(composition, pd.Series)
        int_list = _ELEMENTS[ion.element]
        if ion.isotope in int_list:
            table = _ELEMENTS[ion.element][ion.isotope].items()
            _interf = [
                v * constrained_abundance_estimate(composition, k) for (k, v) in table
            ]
            if _interf:
                result += np.array(_interf).sum() / composition[str(ion.element)]

    return result


def minimum_interference_isotope(element, composition=None):
    """
    Find the minimum total interference isotope for a given element.
    """

    isos = list(_ELEMENTS[element].keys())
    _interfs = [
        sum_of_interferences(element.add_isotope(iso), composition=composition)
        for iso in isos
    ]

    return isos[_interfs.index(np.min(_interfs))]
