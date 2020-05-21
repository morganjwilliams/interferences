"""
Functions for creating, formatting and serialising representaitons of molecules.
"""
import pandas as pd
import periodictable as pt
from pyrolite.mineral.transform import merge_formulae
from ..util.sorting import get_relative_electronegativity
from ..util.meta import interferences_datafolder
import logging

logging.getLogger(__name__).addHandler(logging.NullHandler())
logger = logging.getLogger(__name__)


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
    mols = [pt.formula(i) for i in df.index.str.strip("+")]
    drop = []
    _charges = [i for i in charges if i > 1]
    for m in mols:
        for i in _charges:
            mult = merge_formulae([m] * i)  # a multiple of this molecule
            if mult in mols:
                drop.append(repr_formula(mult) + "+" * i)
    return drop


def _find_duplicate_indexes(df):
    return df.index[df.index.duplicated(keep="first")]


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
        duplicates = _find_duplicate_indexes(df)
        logger.debug("Dropping duplicate indexes: {}".format(", ".join(duplicates)))
        df.drop_duplicates(
            subset="index", keep="first", inplace=True
        )  # drop any duplicate indexes

    if multiples:
        idx = df.index
        dup_multiples = _find_duplicate_multiples(df, charges=charges)
        if duplicated_mz:
            logger.debug(
                "Dropping multiples (duplicate m_z): {}".format(
                    ", ".join(dup_multiples)
                )
            )
            df = df.loc[idx.difference(dup_multiples), :]  # drop any duplicate m_z
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


def get_molecule_labels(df):
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
    label_src = interferences_datafolder(subfolder="table") / "molecular_labels.csv"
    labels = pd.Series(index=df.index, name="label")
    try:
        _label_index = pd.read_csv(label_src, index_col=0)
        known_labels = _label_index.reindex(index=df.index).dropna().index
        unknown_labels = df.index.difference(known_labels)

        if known_labels.size:
            labels[known_labels] = _label_index["label"]
    except (ValueError, FileNotFoundError):
        _label_index = pd.DataFrame(columns=["label"])
        unknown_labels = df.index  # assume they're all unknown

    # fill in the gaps
    unkwn = df.loc[unknown_labels, "molecule"].apply(
        get_formatted_formula, sorted=True
    ) + df.loc[unknown_labels, "charge"].apply(
        lambda c: r"$\mathrm{^{" + "+" * c + "}}$"
    )
    labels.loc[unknown_labels] = unkwn
    # append new index values to the datafile
    if unknown_labels.size:
        if label_src.exists():
            labels[unknown_labels].sort_index().to_csv(
                label_src, mode="a", header=False, index=True
            )
        else:  # write and create the file with headers
            labels[unknown_labels].sort_index().to_csv(
                label_src, mode="w", header=True, index=True
            )
    return labels


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
