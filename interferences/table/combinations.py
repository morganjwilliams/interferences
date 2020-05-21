"""
Functions for calculating combinations (in the combinatorics sense) of elements and
isotopes into isotope-specified molecular ions.
"""
import pandas as pd
import numpy as np
import periodictable as pt
from collections import Counter
from itertools import product, combinations_with_replacement
from ..util.sorting import get_relative_electronegativity
from .intensity import isotope_abundance_threshold, get_isotopic_abundance_product
from .molecules import molecule_from_components, repr_formula
from ..util.log import Handle

logger = Handle(__name__)


def get_elemental_combinations(elements, max_atoms=3):
    """
    Combine a list of elements into lists of molecular combinations up to a maximum
    number of atoms per molecule. Successively adds smaller molecules until down to
    single atoms.

    Parameters
    ----------
    elements : :class:`list`
        Elements or isotopes to combine into molecules.
    max_atoms : :class:`int`
        Maximum number of atoms per molecule. This limits the number of molecules
        returned to the generally most relevant simple molecules.

    Todo
    ----
    Check that isotopes supplied to this function are propogated
    """
    poss_mol_parts = []
    n = max_atoms
    # sorting here should ensure sorted collections later
    elements = sorted(elements, key=get_relative_electronegativity)
    while n:
        components = combinations_with_replacement(elements, n)
        poss_mol_parts += list(components)
        n -= 1
    return poss_mol_parts[::-1]  # backwards so small ones come first


def get_isotopic_combinations(element_comb, threshold=10e-8):
    """
    Take a combination of elements and expand it to generate the potential combinations
    of elements.

    Parameters
    ----------
    element_comb : :class:`list`
        List of elements for which to combine lists of isotopes.
    threshold : :class:`float`
        Threshold below which to ignore low-abundance isotopes.

    Returns
    -------
    :class:`list`
    """
    iso_components = [
        [el.add_isotope(i) for i in el.isotopes]
        if not isinstance(el, pt.core.Isotope)
        else [el]
        for el in element_comb
    ]
    iso_components = [
        isotope_abundance_threshold(lst, threshold=threshold) for lst in iso_components
    ]
    # Counters used for unorderd comparison of lists,
    # otherwise could use list(product(*(isotope_components)))
    iso_counters = [Counter(comb) for comb in product(*(iso_components))]
    # check for duplicates O(n^2) ~ n isn't likely very large, so this might be ok
    iso_combinations = [
        list(c.elements())
        for n, c in enumerate(iso_counters)
        if c not in iso_counters[:n]
    ]
    return iso_combinations


def component_subtable(components, charges=[1, 2], threshold=10e-8):
    """
    Build a sub-table from a set of elemental components.

    Parameters
    ----------
    components : :class:`list`
        List of elements to combine in the subtable.
    charges : :class:`list` ( :class:`int` )
        Ionic charges to include in the model.
    threshold : :class:`float`
        Threshold for isotopic abundance for inclusion of low-abudance/non-stable
        isotopes.

    Returns
    -------
    :class:`pandas.DataFrame`
    """
    df = pd.DataFrame(
        columns=["m_z", "molecule", "components", "mass", "charge", "iso_product",]
    )
    isocombs = get_isotopic_combinations(components, threshold=threshold)
    df["charge"] = np.repeat(charges, len(isocombs))
    df["components"] = isocombs * len(charges)  # multiplied by number of charges
    # calculate the iso_abund_product ##################################################
    df["iso_product"] = df["components"].apply(get_isotopic_abundance_product)
    # build molecules for each set of components #######################################
    # note we could add charge here, but it would have to be assigned to atoms,
    # not formulae! # break up the following process?
    df["molecule"] = df["components"].apply(molecule_from_components)
    df["mass"] = df["molecule"].apply(lambda x: x.mass)
    df["m_z"] = df["mass"] / df["charge"]
    # get a string-based index #########################################################
    df.index = df["molecule"].apply(repr_formula)
    df.index += df["charge"].apply(lambda c: "+" * c)
    # for consistency, we could string-convert object columns here
    return df
