import os
import pandas as pd
import pathlib
from ..util.meta import interferences_datafolder
from ..util.mz import process_window
from .molecules import deduplicate, _find_duplicate_multiples
from ..util.log import Handle

logger = Handle(__name__)

_COMPLEVEL = 4
_COMPLIB = "lzo"
_ITEMSIZES = {"elements": 30, "parts": 40}


def load_store(path=None, complevel=_COMPLEVEL, complib=_COMPLIB, **kwargs):
    """
    Load the interferences HDF store.

    Parameters
    ----------
    path : :class:`str` | :class:`pathlib.Path`
        Path to the store.
    complevel : :class:`int`
        Compression level option for the HDF store. Uncompressed tables can easily
        reach a few hundred MB - this isn't an issue on a local disk, but can be
        limiting for web transfer.
    complib : :class:`str`
        Which compression library to use.

    Returns
    -------
    :class:`pandas.HDFStore`
    """
    path = path or interferences_datafolder(subfolder="table") / "interferences.h5"
    if not path.exists():
        reset_table(
            path=path, complevel=complevel, complib=complib, remove=False
        )  # init table
    store = pd.HDFStore(path, complevel=complevel, complib=complib, **kwargs)
    return store


def lookup_components(identifier, path=None, key="table", window=None, **kwargs):
    """
    Look up a a list of components from the store based on their identifiers.

    Parameters
    ----------
    identifiers : :class:`str`
        Identifiers for the components to look up.
    path : :class:`str` | :class:`pathlib.Path`
        Path to store to search.
    key : :class:`str`
        Key for the table within the store.
    window : :class:`tuple`
        Window for indexing along m/z to return a subset of results.
    drop_first_level : :class:`bool`
        Whether to drop the first level of the index for simplicity.

    Returns
    -------
    :class:`pandas.DataFrame`
    """
    logger.debug("Attempting identifier lookup.")
    window = process_window(window)
    name = "/" + key

    multi = ""
    if isinstance(identifier, str):
        multi_lookup = False
    elif isinstance(identifier, (list, pd.Index)) and len(identifier) == 1:
        multi_lookup = False
        identifier = identifier[0]
    else:
        multi_lookup = True
        multi = "multi-"
    # try:
    with load_store(path, **kwargs) as store:
        where = []
        empty = False
        if not multi_lookup:
            where += ["elements == '{}'".format(identifier)]
            empty = store.select(name, where=" & ".join(where)).empty
        if window:  # add the m_z window information
            where += ["m_z >= {:5f} & m_z <= {:5f}".format(*window)]

        msg = "Performing {}lookup".format(multi)
        if where:
            msg += " & ".join(where)
        logger.debug(msg)

        if not empty:
            df = store.select(name, where=" & ".join(where))
        else:
            raise IndexError("Identifer(s) not in table.")

        if multi_lookup:
            tbl_idents = pd.unique(df.index.droplevel("parts"))
            df = df.loc[[i for i in identifier if i in tbl_idents], :]

    return df
    # except KeyError:
    #    raise KeyError("Key not in HDFStore.")


def _get_default_multiindex():
    """
    Build an empty multi-index for the table.

    Returns
    -------
    :class:`pandas.MultiIndex`
    """
    return pd.MultiIndex.from_product([[], []], names=["elements", "parts"])


def get_store_index(path, drop_first_level=True, **kwargs):
    """
    """
    with pd.HDFStore(path, **kwargs) as store:
        if "/table" in store.keys():
            index = store.select("/table", columns=["elements", "parts"]).index
        else:
            index = _get_default_multiindex()  # empty index
    if drop_first_level:
        index = index.droplevel("elements")
    return index


def process_subtables(
    dfs,
    charges=None,
    dump=True,
    path=None,
    mode="a",
    data_columns=["elements", "m_z", "iso_abund_product"],
    complevel=_COMPLEVEL,
    complib=_COMPLIB,
    **kwargs
):
    """
    Process and optionally dump a set of subtables to file,
    appending to the hierarchically-indexed table.

    Parameters
    ----------
    dfs : :class:`list`(:class:`pandas.DataFrame`)
        Dataframes to dump.
    charges : :class:`list`
        Charges used to create for the table.
    path : :class:`str` | :class:`pathlib.Path`
        Path to the file to add the table to.
    mode : :class:`str`
        Mode for accessing the HDF file.
    data_columns : :class:`list`
        List of columns to create an indexes for to allow query-by-data.
    complevel : :class:`int`
        Compression level option for the HDF store. Uncompressed tables can easily
        reach a few hundred MB - this isn't an issue on a local disk, but can be
        limiting for web transfer.
    complib : :class:`str`
        Which compression library to use.

    Returns
    -------
    :class:`pandas.DataFrame`
        De-duplicated concatenated version of new tables.
    """
    path = path or interferences_datafolder(subfolder="table") / "interferences.h5"
    logger.debug("Checking Store")
    current_index = get_store_index(path).to_list()
    logger.debug("Combining DataFrames")
    df = pd.concat(dfs, axis=0, ignore_index=False)
    df.index.rename("parts", inplace=True)
    df["elements"] = [id for d in dfs for id in [d.name] * d.index.size]
    ####################################################################################
    logger.debug("Deduplicating")
    output = df.loc[~df.index.duplicated(keep="first"), :]  # remove duplicated indexes
    # take the index from df, and the index from the store and combine them to dedupe
    new_duplicates = _find_duplicate_multiples(
        pd.DataFrame(index=output.index.to_list() + current_index), charges=charges
    )
    if len(new_duplicates):
        # all of these should be in the output.index, so we can just drop them
        logger.debug("Removing duplicates: {}".format(", ".join(new_duplicates)))
        output.drop(index=new_duplicates, inplace=True)

    if dump:
        logger.debug("Reindexing")
        # create hierarchical indexes for a copy of the table to dump into the store
        to_store = output.set_index("elements", append=True)
        to_store = to_store.reorder_levels(["elements", "parts"], axis=0)
        # convert non-string. non-numerical objects to string
        # append to the existing dataframe
        # somehow S[34]S[34]++ sneaks past
        logger.debug(
            "Dumping {} tables to HDF store.".format(
                ",".join(pd.unique(to_store.index.get_level_values("elements")))
            )
        )
        to_store.to_hdf(
            path,
            key="table",
            mode="a",
            append=True,
            format="table",
            data_columns=data_columns,
            min_itemsize=_ITEMSIZES,
            complevel=_COMPLEVEL,
            complib=_COMPLIB,
        )
    return output


def reset_table(
    path=None,
    remove=True,
    key="table",
    format="table",
    complevel=_COMPLEVEL,
    complib=_COMPLIB,
    **kwargs
):
    """
    Reset or remove a HDF store.

    Parameters
    ----------
    path : :class:`str` | :class:`pathlib.Path`
        Path to store.
    remove : :class:`bool`
        Whether to remove the table from disk, if possible.
    format : :class:`str`
        Format to set for the new tables.
    complevel : :class:`int`
        Compression level option for the HDF store. Uncompressed tables can easily
        reach a few hundred MB - this isn't an issue on a local disk, but can be
        limiting for web transfer.
    complib : :class:`str`
        Which compression library to use.
    """
    path = path or interferences_datafolder(subfolder="table") / "interferences.h5"
    if not path.parent.exists():
        logger.debug("Creating folder for store.")
        path.parent.mkdir(parents=True)  # ensure directory exists
    if remove:
        logger.debug("Removing store.")
        try:
            os.remove(path)  # remove the file
        except FileNotFoundError:
            logger.debug("Store already removed or not present.")
    else:  # keep table keys, set them to empty frames
        logger.debug("Resetting store table: {}/{}".format(path.name, key))
        df = pd.DataFrame(
            index=_get_default_multiindex(),
            columns=["m_z", "mass", "charge", "iso_product",],
        )
        # note - this will not work until a pytables bug is fixed,
        # where the table doesnt' generate from an empty frame.
        df.to_hdf(
            path,
            key=key,
            format=format,
            mode="w",
            append=True,
            complevel=complevel,
            complib=complib,
            min_itemsize=_ITEMSIZES,
            **kwargs
        )
