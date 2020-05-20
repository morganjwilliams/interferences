import pandas as pd
import numpy as np
import periodictable as pt
from .molecules import (
    molecule_from_components,
    get_molecule_labels,
    repr_formula,
    _get_isotope,
)
from .store import (
    load_groups_store,
    dump_element_group,
    lookup_component_subtable,
    consoliate_groups,
)
from .intensity import isotope_abundance_threshold, get_isotopic_abundance_product
from .combinations import get_elemental_combinations, component_subtable
from ..util.sorting import get_first_atom
from ..util.mz import process_window
import logging

logging.getLogger(__name__).addHandler(logging.NullHandler())


def build_table(
    elements=None,
    max_atoms=3,
    sortby=["m_z", "charge", "mass"],
    charges=[1, 2],
    add_labels=False,
    threshold=10e-8,
    window=None,
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
    mass_window : :class:`tuple`
        Window of interest to filter out irrelevant examples (here a mass window,
        which directly translates to m/z window with z=1).

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
        columns=["m_z", "molecule", "components", "mass", "charge", "iso_product",]
    )
    # set numeric datatypes
    table = table.astype(
        {"mass": "float", "charge": "int8", "iso_product": "float", "m_z": "float",}
    )
    window = process_window(window)
    # build up combinations of elements, forming the components column
    molecular_components = get_elemental_combinations(
        elements, max_atoms=max_atoms
    )  # this can't be split easily
    groupstore = load_groups_store()
    updated_groups = False
    for components in molecular_components:
        # convert elements to pt.core.Element
        components = [
            get_first_atom(el) if isinstance(el, str) else el for el in components
        ]
        identifier = "-".join([repr(c) for c in components])
        try:  # check if these are in the database
            df = lookup_component_subtable(groupstore, identifier, window=window)
        except (KeyError, IndexError):  # if not, build it
            # create the whole table, ignoring window, to dump into refernce.
            df = component_subtable(components, charges=charges)
            dump_element_group(
                df, identifier, path=groupstore
            )  # append to the HDF store for later use
            updated_groups = True

            # filter for window here

        table = pd.concat([table, df], axis=0, ignore_index=False)
    groupstore.close()
    if updated_groups:
        consoliate_groups()  # update the interferences.h5 file with new data
    # filter out invalid entries, eg. H{2+} ############################################
    # TODO
    # sort table #######################################################################
    table.sort_values(sortby, axis="index", inplace=True)
    # additional columns ###############################################################
    if add_labels:  # this step is string-operation intensive, and hence very slow
        table["label"] = get_molecule_labels(table)
    # for consistency with prevsiouly serialized data:
    return table.astype({"molecule": str, "components": str})
