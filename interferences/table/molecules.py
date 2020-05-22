"""
Functions for creating, formatting and serialising representaitons of molecules.
"""
import re
import pandas as pd
import numpy as np
import periodictable as pt
from pyrolite.mineral.transform import merge_formulae
from ..util.sorting import get_relative_electronegativity
from ..util.meta import interferences_datafolder
from ..util.log import Handle

logger = Handle(__name__)

_COMPLEVEL = 4
_COMPLIB = "lzo"
_ITEMSIZES = {"label": 50, "index": 40}


def components_from_index_value(idx):
    return re.findall(r"\w+\[\d+\]", idx)


def _find_duplicate_multiples(df, charges=None):
    """
    Remove multiples of moleclues which have the same m/z (e.g. OH+, H2O2++).

    Parameters
    ----------
    df : :class:`pandas.DataFrame`
        Dataframe to check the index of.
    charges : :class:`list`
        List of valid charges for the frame.

    Returns
    -------
    :class:`list:
    """
    counts = df.index.map(lambda s: s.count("["))
    target_charges = [c for c in np.arange(np.max(charges)) + 1 if c // 2 == c / 2]
    source_n_atoms = [c for c in np.arange(counts.max()) + 1 if c <= (counts.max() / 2)]

    drop_mols = []
    for n_atoms in source_n_atoms:
        src = df.index[counts == n_atoms]  # get e.g. 1-atom molecules
        for m in src.str.strip("+"):
            potential_multiples = [
                repr_formula(merge_formulae([m] * c)) + "+" * c for c in target_charges
            ]

            drop_mols += df.index.intersection(potential_multiples).to_list()
    return drop_mols


def deduplicate(df, charges=None, multiples=True):
    """
    De-duplicate a dataframe index based on index values and and molecule-multiples.

    Parameters
    ----------
    df : :class:`pandas.DataFrame`
        Dataframe to check the index of.
    charges : :class:`list`
        List of valid charges for the frame.
    multiples : :class:`bool`
        Whether to remove molecule-multiples.

    Returns
    -------
    :class:`pandas.DataFrame`
    """
    # remove duplicate m/z #############################################################
    idx = df.index
    if idx.duplicated().any():
        duplicates = df.index[df.index.duplicated(keep="first")]
        logger.debug("Dropping duplicate indexes: {}".format(", ".join(duplicates)))
        df.drop_duplicates(
            subset="index", keep="first", inplace=True
        )  # drop any duplicate indexes

    if multiples:
        dup_multiples = _find_duplicate_multiples(df, charges=charges)
        if dup_multiples:
            logger.debug(
                "Dropping multiples (duplicate m_z): {}".format(
                    ", ".join(dup_multiples)
                )
            )
            df.drop(dup_multiples, axis=0, inplace=True)  # drop any duplicate m_z
    return df


def _get_isotope(element):
    """
    Parameters
    ----------
    element : :class:`periodictable.core.Element`
        Element or isotope.
    Returns
    -------
    :class:`int`
    """
    try:
        return element.isotope
    except AttributeError:
        return 0


def repr_formula(molecule):
    """
    Get a string representation of a formula which preserves element and isotope
    information.
    """
    parts = [
        "{}".format(repr(el)) if cnt == 1 else "{}".format(repr(el)) * cnt
        for el, cnt in molecule.atoms.items()
    ]
    return "".join(parts)


def get_formatted_formula(molecule, sorted=False):
    """
    Construct a formatted name for a molecule.

    Parameters
    -----------
    molecule : :class:`~periodictable.formulas.Formula`
        Molecule to name.
    sorted : :class:`bool`
        Whether a molecular formula is already sorted, so sorting can
        be skipped.

    Returns
    -------
    :class:`str`
    """
    molecule = pt.formula(molecule)
    components = list(molecule.atoms.keys())
    if not sorted:
        components = sorted(
            components,
            key=lambda x: (get_relative_electronegativity(x), _get_isotope(x)),
        )
    name = r"$\mathrm{"  # remove italicized text effect
    for c in components:
        part = ""
        if hasattr(c, "isotope"):
            part += "^{" + "{}".format(c.isotope) + "}"  # superscript isotope
        part += str(c.element)
        count = molecule.atoms[c]
        if count > 1:
            part += "_{" + "{:d}".format(molecule.atoms[c]) + "}"
        name += part
    name += "}$"  # finish TeX formatting
    return name


def get_molecule_labels(df, **kwargs):
    """
    Get labels for molecules based on their composition and charge.

    Parameters
    -----------
    df : :class:`pandas.DataFrame`

    Returns
    -------
    :class:`pandas.Series`
    """
    # look up index values which are pre-computed
    label_src = interferences_datafolder(subfolder="table") / "labels.h5"
    labels = pd.DataFrame(index=df.index, columns=["label"])
    try:
        with pd.HDFStore(
            label_src, complevel=_COMPLEVEL, complib=_COMPLIB, **kwargs
        ) as store:
            label_store = store.select("/table")

        known = label_store.index.intersection(df.index)
        unknown = df.index.difference(known)
        if known.size:
            labels.loc[known, "label"] = label_store["label"]

    except (KeyError, FileNotFoundError):
        label_store = pd.DataFrame(columns=["label"])
        unknown = df.index  # assume they're all unknown

    if unknown.size:
        logger.debug("Buiding {} labels.".format(unknown.size))
        # fill in the gaps

        mols = unknown.map(lambda x: get_formatted_formula(x.strip('+'), sorted=True))
        charges = df.loc[unknown, "charge"].apply(
            lambda c: r"$\mathrm{^{" + "+" * c + "}}$"
        )
        labels.loc[unknown, "label"] = mols + charges
        # append new index values to the datafile
        logger.debug("Dumping {} labels to file.".format(unknown.size))

        if label_src.exists():
            labels.loc[unknown].to_hdf(
                label_src,
                key="table",
                mode="a",
                append=True,
                format="table",
                min_itemsize=_ITEMSIZES,
                complevel=_COMPLEVEL,
                complib=_COMPLIB,
            )
        else:  # write and create the file with headers
            labels.loc[unknown].to_hdf(
                label_src,
                key="table",
                mode="w",
                append=True,
                format="fixed",
                min_itemsize=_ITEMSIZES,
                complevel=_COMPLEVEL,
                complib=_COMPLIB,
            )
    return labels


def molecule_from_components(components):
    """
    Builds a :class:`~periodictable.formulas.Formula` from a list of atom or
    isotope components.

    Parameters
    ----------
    components : :class:`list`
        Atomic, isotope or molecular components to construct an ionic molecule from.

    Returns
    -------
    :class:`~periodictable.formulas.Formula`

    Todo
    -----
    * Modify to accept consumption of molecular components (e.g. Fe2O3+)

    See Also
    ---------
    :func:`pyrolite.mineral.transform.merge_formulae`
    """
    return merge_formulae(components)
