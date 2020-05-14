import pandas as pd
import numpy as np
import periodictable as pt
from itertools import product, combinations_with_replacement
from collections import Counter
from .molecules import (
    molecule_from_components,
    get_molecule_labels,
    repr_formula,
    _get_isotope,
)
from .store import load_store, dump_table, get_table
from ..util.sorting import get_relative_electronegativity
import logging

logging.getLogger(__name__).addHandler(logging.NullHandler())


def isotope_abundance_threshold(isotopes, threshold=10e-8):
    """
    Remove isotopes from a list which have no or zero abundance.

    Parameters
    ----------
    isotopes : :class:`list`

    Returns
    -------
    :class:`list`
    """
    constraint_functions = [
        lambda x: hasattr(x, "abundance"),
        lambda x: x.abundance > threshold,
    ]
    for valid in constraint_functions:
        isotopes = [i for i in isotopes if valid(i)]
    return isotopes


def get_isotopic_abund_product(components):
    """
    Estimates the abundance of a molecule based on the abundance of the isotopic
    components.

    Returns
    -------
    :class:`float`

    Notes
    ------
    This is essentially a simplistic activity model.
    Isotopic abundances from periodictable are in %, and are hence divded by 100 here.
    """
    abund = 1.0
    for iso in components:
        abund *= iso.abundance / 100.0
    return abund


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
    return poss_mol_parts


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
        columns=[
            "m_z",
            "molecule",
            "components",
            "mass",
            "charge",
            "iso_abund_product",
        ]
    )
    # take these elemental combinations and exapand to a list of isotopes ##########
    isotope_components = [
        [el.add_isotope(i) for i in el.isotopes]
        if not isinstance(el, pt.core.Isotope)
        else [el]
        for el in components
    ]
    isotope_components = [
        isotope_abundance_threshold(lst, threshold=threshold)
        for lst in isotope_components
    ]
    # get isotopic combinations ####################################################
    # Counters used for unorderd comparison of lists,
    # otherwise could use list(product(*(isotope_components)))
    isocombs = [Counter(comb) for comb in product(*(isotope_components))]
    # check for duplicates O(n^2) ~ n isn't likely very large, so this might be ok
    isocombs = [
        list(c.elements()) for n, c in enumerate(isocombs) if c not in isocombs[:n]
    ]
    df["charge"] = np.repeat(charges, len(isocombs))
    df["components"] = isocombs * len(charges)  # multiplied by number of charges
    # calculate the iso_abund_product ##################################################
    df["iso_abund_product"] = df["components"].apply(get_isotopic_abund_product)
    # build molecules for each set of components #######################################
    # note we could add charge here, but it would have to be assigned to atoms,
    # not formulae! # break up the following process?
    df["molecule"] = df["components"].apply(molecule_from_components)
    df["mass"] = df["molecule"].apply(lambda x: x.mass)
    df["m_z"] = df["mass"] / df["charge"]
    # remove duplicate m/z #############################################################
    df.drop(df.index[df["m_z"].duplicated()], axis="index", inplace=True)
    # get a string-based index #########################################################
    df.index = df["molecule"].apply(repr_formula)
    df.index += df["charge"].apply(lambda c: "+" * c)
    # for consistency, we could string-convert object columns here
    return df


def build_table(
    elements=None,
    max_atoms=3,
    sortby=["m_z", "charge", "mass"],
    charges=[1, 2],
    add_labels=False,
    threshold=10e-8,
    complevel=3,
):
    """
    Build the interferences table.

    Parameters
    ----------
    elements : :class:`list`
        List of elements to include in the table.
    max_atoms : :class:`int`
        Largest size of molecule to build, in atoms.
    sortby : :class:`str` | :class:`list`
        Column or list of columns to sort the final table by.
    charges : :class:`list` ( :class:`int` )
        Ionic charges to include in the model.
    add_labels : :class:`bool`
        Whether to produce molecule names which are nicely formatted. This takes
        additional computation time.
    threshold : :class:`float`
        Threshold for isotopic abundance for inclusion of low-abudance/non-stable
        isotopes.

    Todo
    -----
    Consider options for parellizing this to reduce build time. This would allow
    larger molecules to be included.

    Invalid molecules (e.g. `H{2+}`) will currently be present, but will ideally be
    filtered out

    In some cases, mass peaks will be duplicated, and we want to keep the simplest
    version (e.g. `Ar[40]+` and `Ar[40]2{2+}`). We here remove duplicate mass peaks
    before sorting (i.e. take the first one, as higher charges would be penalised),
    but we could potentially add a check that both contain the same isotopic
    components for verificaiton (this would be slow..).

    While "m/z" would be an appropriate column name, it can't be used in HDF indexes.
    """

    table = pd.DataFrame(
        columns=[
            "m_z",
            "molecule",
            "components",
            "mass",
            "charge",
            "iso_abund_product",
        ]
    )
    # set numeric datatypes
    table = table.astype(
        {
            "mass": "float",
            "charge": "int8",
            "iso_abund_product": "float",
            "m_z": "float",
        }
    )
    # build up combinations of elements, forming the components column
    molecular_components = get_elemental_combinations(
        elements, max_atoms=max_atoms
    )  # this can't be split easily
    store = load_store()
    for components in molecular_components:
        # convert elements to pt.core.Element
        components = [
            getattr(pt, el) if isinstance(el, str) else el for el in components
        ]
        identifier = "-".join([repr(c) for c in components])
        try:  # check if these are in the database
            df = get_table(store, identifier)
        except KeyError:  # if not, build it
            df = component_subtable(components, charges=charges)
            dump_table(df, identifier)  # append to the HDF store for later use

        table = pd.concat([table, df], axis=0, ignore_index=False)
    store.close()
    # filter out invalid entries, eg. H{2+} ############################################
    # TODO
    # sort table #######################################################################
    table.sort_values(sortby, axis="index", inplace=True)
    # additional columns ###############################################################
    if add_labels:  # this step is string-operation intensive, and hence very slow
        table["label"] = get_molecule_labels(df)
    return table
