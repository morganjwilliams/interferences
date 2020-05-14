from itertools import groupby, product, combinations_with_replacement
import numpy as np
import pandas as pd
import periodictable as pt
from periodictable.core import isisotope, iselement, isatom
import matplotlib.pyplot as plt
from collections import defaultdict
from copy import deepcopy
from pyrolite.geochem.norm import get_reference_composition

_MASS = defaultdict(dict)
_ELEMENTS = defaultdict(dict)


def constrained_abundance_estimate(composition, formula):
    """
    Get an abundance estimate for a specific molecule constrained by a
    starting composition.

    Parameters
    ----------
    composition : :class:`pandas.Series`
        Composition of the target.
    formula : :class:`~periodictable.formulas.Formula`
        Formula of the molecule.

    Returns
    -------
    :class:`float`
    """

    if isinstance(formula, pt.core.Element):
        return getattr(composition, str(formula.element), 0.0)
    elif isinstance(formula, pt.formulas.Formula):
        prod = 1.0
        for a in formula.atoms.keys():
            prod *= getattr(composition, str(a.element), 0.0)
        return prod
    elif isinstance(formula, str):
        return constrained_abundance_estimate(composition, pt.formula(compound=formula))
    else:

        print(formula, type(formula))
        raise AssertionError


def get_reference_abundance(molecule, reference="Chondrite_PON", unknown_val=1000):
    """
    Get a reference abundance for a molecule based on a specified reference composition.

    Parameters
    ----------
    molecule : :class:`~periodictable.formulas.Formula`

    reference : :class:`str`
        Reference composition to calculate molecular abundance for.
    """
    norm = get_reference_composition(reference).comp  # in ppm
    abund = 10 ** 6  # 100%
    for iso in molecule.atoms.keys():
        el = get_element_name(iso)
        abund *= getattr(norm, el)
    if not np.isfinite(abund):
        abund = unknown_val  # unknown abundance%
    return abund


def sum_of_interferences(ion, composition: pd.Series = None):
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
