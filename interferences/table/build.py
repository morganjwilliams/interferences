import pandas as pd
import numpy as np
import periodictable as pt
from tqdm import tqdm
from pyrolite.util.meta import ToLogger
from .molecules import (
    molecule_from_components,
    get_molecule_labels,
    repr_formula,
    _get_isotope,
    deduplicate,
)
from .store import (
    load_store,
    process_subtables,
    lookup_components,
)
from .intensity import isotope_abundance_threshold, get_isotopic_abundance_product
from .combinations import get_elemental_combinations, component_subtable
from ..util.sorting import get_first_atom
from ..util.mz import process_window
from ..util.log import Handle

logger = Handle(__name__)


def build_table(
    elements=None,
    max_atoms=3,
    sortby=["m_z", "charge", "mass"],
    charges=[1, 2],
    add_labels=False,
    threshold=None,
    window=None,
    cache_results=True,
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
    cache_results : :class:`bool`
        Whether to store the results on disk for

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
    table = pd.DataFrame(columns=["m_z", "mass", "charge", "iso_product",])
    # set numeric datatypes
    table = table.astype(
        {"mass": "float", "charge": "int8", "iso_product": "float", "m_z": "float",}
    )
    window = process_window(window)
    # build up combinations of elements, forming the components column
    # this can't be split easily
    combinations = get_elemental_combinations(elements, max_atoms=max_atoms)
    logger.info("Building {:d} component combinations.".format(len(combinations)))
    # convert elements to pt.core.Element
    combinations = [
        [get_first_atom(el) if isinstance(el, str) else el for el in components]
        for components in combinations
    ]
    identifiers = [
        "-".join([repr(c) for c in components]) for components in combinations
    ]
    # do a loookup here for all the identifiers, then go build the unkonwn ones
    cached_combinations = []
    try:
        lookup = lookup_components(identifiers)  # ignore window here
        cached_combinations = list(pd.unique(lookup.index.get_level_values("elements")))
        lookup = lookup.droplevel("elements")
        if window is not None:  # process_window for lookup
            lookup = lookup.loc[lookup.m_z.between(*window)]
        table = pd.concat([table, lookup], axis=0, ignore_index=False)
    except KeyError as e:
        pytables_expect = "No object named /table in the file"
        if pytables_expect in str(e):
            logger.debug("Table not built - Defaulting to build-all.")
            build = identifiers
        else:
            raise e
    except IndexError as e:
        logger.debug("No identifiers in table - Defaulting to build-all.")
        build = identifiers

    if threshold is not None:
        # build them from scratch, as it's more difficult to filter after the fact
        cached_combinations = []

    beyond_bounds = []
    if window is not None:  # check potential m_z relevance
        # check whether mz is within margin of target
        margin = 0.10  # 10%
        in_mass_bounds = lambda ID: any(
            [
                window[0] * (1 - margin)
                < pt.formula(ID.replace("-", "")).mass / c
                < window[1] * (1 + margin)
                for c in charges
            ]
        )
        beyond_bounds = []
        beyond_bounds = [ID for ID in identifiers if not in_mass_bounds(ID)]
        if beyond_bounds:
            logger.debug(
                "Skipping tables outside m/z bounds {}.".format(",".join(beyond_bounds))
            )

    need_to_build = [
        (ID, components)
        for (ID, components) in zip(identifiers, combinations)
        if (ID not in cached_combinations + beyond_bounds)
    ]

    new_tables = []
    if need_to_build:
        progressbar = tqdm(need_to_build)  # file=ToLogger(logger)
        barwidth = 16 + 3 * max_atoms
        progressbar.set_description(" " * barwidth)
        for ID, components in progressbar:
            # create the whole table, ignoring window, to dump into refernce.
            df = component_subtable(components, charges=charges, threshold=threshold)
            df.name = ID
            msg = "{} @ {:d} rows".format(ID, df.index.size)
            msg += " " * (barwidth - len(msg))
            progressbar.set_description(msg)
            logger.debug("Building table for {} @ {:d} rows".format(ID, df.index.size))

            new_tables.append(df)

    if new_tables:  # append new dfs to the HDF store for later use
        # if we use threshold, we'll put an incomplete table into the reference store
        additions = process_subtables(
            new_tables, charges=charges, dump=(cache_results and (threshold is None))
        )
        # should de-duplicate the new_tables in this table
        # could rearrange and return deduped tables from dump_subtables
        if window is not None:
            additions = additions.loc[additions["m_z"].between(*window), :]
        table = pd.concat([table, additions], axis=0, ignore_index=False,)

    # filter out invalid entries, eg. H{2+} ############################################
    # TODO
    # sort table #######################################################################
    table.sort_values(sortby, axis="index", inplace=True)
    # additional columns ###############################################################
    if add_labels:  # this step is string-operation intensive, and hence very slow
        logger.info("Adding labels to the table.")
        table["label"] = get_molecule_labels(table)
    # for consistency with prevsiouly serialized data:
    return table
