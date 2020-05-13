import pandas as pd
import numpy as np
import periodictable as pt
from itertools import product, combinations_with_replacement
from .molecules import molecule_from_components, get_molecular_formula
from ..util.meta import interferences_datafolder
from ..util.sorting import get_relative_electronegativity


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
    while n:
        components = combinations_with_replacement(elements, n)
        poss_mol_parts += list(components)
        n -= 1
    return poss_mol_parts


def build_table(
    elements=None,
    max_atoms=3,
    sortby="m/z",
    abundance_threshold=10e-8,
    charges=[1, 2],
    add_labels=False,
):
    """
    Build the interferences table from scratch.

    Parameters
    ----------
    elements : :class:`list`
        List of elements to include in the table.
    max_atoms : :class:`int`
        Largest size of molecule to build, in atoms.
    sortby : :class:`str`
        Column to sort the final table by.
    abundance_threshold : :class:`float`
        Threshold for isotopic abundance for inclusion of low-abudance/non-stable
        isotopes.
    charges : :class:`list` ( :class:`int` )
        Ionic charges to include in the model.
    add_labels : :class:`bool`
        Whether to produce molecule names which are nicely formatted. This takes
        additional computation time.

    Todo
    -----
    Consider options for parellizing this to reduce build time. This would allow
    larger molecules to be included.

    This could be more efficient if the isotopes are expanded first.

    In some cases, mass peaks will be duplicated, and we want to keep the simplest
    version (e.g. Ar[40]+ and Ar[40]2{2+}). We here remove duplicate mass peaks
    before sorting (i.e. take the first one, as higher charges would be penalised),
    potentially with a check that both contain the same isotopic components.
    """

    table = pd.DataFrame(
        columns=[
            "m/z",
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
            "mass": "float64",
            "charge": "int",
            "iso_abund_product": "float32",
            "m/z": "float64",
        }
    )
    # convert elements to pd.Element
    elements = [getattr(pt, el) for el in elements]
    # build up combinations of elements, forming the components column
    molecular_components = get_elemental_combinations(
        elements, max_atoms=max_atoms
    )  # this can't be split easily

    for components in molecular_components:
        # check if these are in the database
        df = pd.DataFrame(columns=table.columns)
        # take these elemental combinations and add the isotopic axis ##################
        isotope_components = [
            isotope_abundance_threshold(
                [el.add_isotope(i) for i in el.isotopes], threshold=abundance_threshold
            )
            for el in components
        ]
        # get isotopic combinations ####################################################
        # NOTE -  itertools will return a generator for product
        #  should this be combinations, not product??
        isocombs = [comb for comb in product(*(isotope_components))]
        # need to do complex list comparison here to avoid duplicates ...

        df["components"] = isocombs * len(charges)  # multiplied by number of charges
        df["charge"] = np.repeat(charges, len(isocombs))

        # build molecules for each set of components ###################################
        # note we could add charge here, but it would have to be assigned to atoms,
        # not formulae! # break up the following process?
        df["molecule"] = df["components"].apply(molecule_from_components)
        # get a string-based index #####################################################
        # adding the charge to the index repr string means that we would need to
        # remove this if it were to be again parsed by periodictable

        # append the table #############################################################
        table = pd.concat([table, df], axis=0, ignore_index=True)

    # filter out invalid entries, eg. H{2+}

    # calculate the iso_abund_product
    table["iso_abund_product"] = table["components"].apply(get_isotopic_abund_product)
    table["mass"] = table["molecule"].apply(lambda x: x.mass)
    table["m/z"] = table["mass"] / table["charge"]
    # calculate abundances
    # sort df by m/z >
    table.sort_values(["m/z", "charge", "mass"], axis="index", inplace=True)
    # remove duplicate m/z
    # dup_rows = table.index[table["m/z"].value_counts() > 1]

    table.drop(table.index[table["m/z"].duplicated()], axis="index", inplace=True)

    # clean up index names
    table.index = table["molecule"].apply(str)
    table.index += table["charge"].apply(lambda c: "+" * c)
    if add_labels:  # this step is string-operation intensive, and hence very slow
        # look up index values which are pre-computed
        src = interferences_datafolder(subfolder="table") / "molecular_labels.json"
        try:
            _label_index = pd.read_json(src)
            known_labels = _label_index.reindex(index=table.index).dropna().index
            unknown_labels = table.index.difference(known_labels)

            if known_labels.size:
                table.loc[known_labels, "label"] = _label_index["label"]
        except (ValueError, FileNotFoundError):
            _label_index = pd.DataFrame(columns=["label"])
            unknown_labels = table.index  # assume they're all unknown
        # fill in the gaps
        table.loc[unknown_labels, "label"] = table.loc[
            unknown_labels, "molecule"
        ].apply(get_molecular_formula) + table["charge"].apply(
            lambda c: r"$^{" + "+" * c + "}$"
        )

        # write out new index values to the datafile
        pd.concat([_label_index, table.loc[unknown_labels, ["label"]]]).to_json(src)
    return table


def dump_table(df, dirpath=None):
    """
    Dump the interferences table to file.

    Parameters
    ----------
    df : :class:`pandas.DataFrame`
        Dataframe to dump.
    dirpath : :class:`str` | :class:`pathlib.Path`
        Path to the directory to save the file in.
    """
    dirpath = dirpath or interferences_datafolder(subfolder="table")
    filepath = dirpath / "interferences.csv"
    df.to_csv(filepath)


def load_table(dirpath=None):
    """
    Dump the interferences table to file.

    Parameters
    ----------
    dirpath : :class:`str` | :class:`pathlib.Path`
        Path to the directory to save the file in.

    Returns
    -------
    :class:`pandas.DataFrame`
    """
    dirpath = dirpath or interferences_datafolder(subfolder="table")
    filepath = dirpath / "interferences.csv"
    df = pd.read_csv(filepath, index=False)
    #
