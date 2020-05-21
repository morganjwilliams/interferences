import periodictable as pt
import numpy as np
import pandas as pd
from interferences.util.pt import get_periodic_frame
from .log import Handle

logger = Handle(__name__)


def _build_relative_electronegativities(reverse=True):
    """
    Construct a dictionary of ordering of electronegativity across the elements,
    according to IUPAC nomenclature.

    Parameters
    ----------
    reverse : :class:`bool`
        Whether to reverse the ordering of the list (i.e such that low numbers, which
        would occur first, correspond to less electronegative elements) allowing
        electronegative elements to occur last (or, on the RHS of formulae).

    Returns
    -------
    :class:`dict`

    References
    ----------
    https://en.wikipedia.org/wiki/IUPAC_nomenclature_of_inorganic_chemistry_2005
    """
    pt_df = get_periodic_frame()
    ordering = []
    ordering += pt_df[17].to_list()
    ordering += pt_df[16].to_list()
    ordering += [pt.H]
    for grp in np.arange(4, 16)[::-1]:
        ordering += pt_df[grp].to_list()
    ordering += pt_df.loc[:5, 3].to_list()  # Sc, Y
    ordering += pt_df.loc[6, 3]  # lanthanoids
    ordering += pt_df.loc[7, 3]  # actinoids
    ordering += pt_df[2].to_list()
    ordering += pt_df.loc[1:, 1].to_list()  # grp one excluding H
    ordering += pt_df[18].to_list()  # noble gases
    ordering = [el for el in ordering if not pd.isnull(el)]
    if reverse:
        ordering = ordering[::-1]
    return {el.number: ix for ix, el in enumerate(ordering)}


def get_first_atom(molecule):
    """
    Get the first atom in a molecular formula.

    Parameters
    ----------
    molecule : :class:`~periodictable.core.Element` | :class:`~periodictable.formulas.Formula`
        Molecule to check.

    Returns
    -------
    :class:`~periodictable.core.Element`
        Element or isotope.
    """
    if isinstance(molecule, pt.core.Element):
        return molecule
    else:
        return list(pt.formula(molecule).atoms.keys())[0]


def get_relative_electronegativity(element, reverse=True):
    """
    Get an index of the relative electronegativity of an element, for use in
    sorting elements (e.g. for chemical formulae). If a list of elements is supplied,
    a list will be returned.

    Parameters
    ----------
    element : :class:`str` | :class:`periodictable.core.Element` | :class:`list`

    Returns
    -------
    :class:`int` | :class:`list`

    Note
    -----
    Electronegativity check uses numbers as these are provided by both
    Element and Isotope objects.
    """
    en = _build_relative_electronegativities(reverse=reverse)
    if isinstance(element, list):
        return [en[get_first_atom(pt.formula(e)).number] for e in element]
    else:
        return en[get_first_atom(pt.formula(element)).number]


get_relative_electronegativity("Ca[40]")
