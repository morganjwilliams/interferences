import os
import pandas as pd
import pathlib
from ..util.meta import interferences_datafolder
from ..util.mz import process_window
from .molecules import deduplicate, _find_duplicate_multiples
import logging

logging.getLogger(__name__).addHandler(logging.NullHandler())
logger = logging.getLogger(__name__)

_COMPLEVEL = 4
_COMPLIB = "lzo"
_ITEMSIZES = {"elements": 30, "parts": 40, "components": 40, "molecule": 30}


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


def lookup_component_subtable(
    store, identifier, key="table", window=None, drop_first_level=True
):
    """
    Look up a component-subtable from the store based on an identifier.

    Parameters
    ----------
    store : :class:`pandas.HDFStore`
        Store to search.
    identifier : :class:`str`
        Identifier for the subtable.
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
    logger.debug("Attempting lookup for Identifer: {}".format(identifier))
    window = process_window(window)
    name = "/" + key
    if name in store.keys():
        where = "elements == '{}'".format(identifier)
        if not store.select(name, where=where).empty:
            if window:  # add the m_z window information
                where += " & m_z >= {:5f} & m_z <= {:5f}".format(*window)
            logger.debug("Performing lookup where: " + where)
            # get the sub-table, and drop the extra index level for simplicity
            if drop_first_level:
                return store.select(name, where=where).droplevel("elements")
            else:
                return store.select(name, where=where)
        else:
            raise IndexError("Identifer not in table.")
    raise KeyError("Key not in HDFStore.")


def dump_subtable(
    df,
    identifier,
    charges=None,
    path=None,
    mode="a",
    data_columns=["elements", "m_z", "iso_abund_product"],
    complevel=_COMPLEVEL,
    complib=_COMPLIB,
    **kwargs
):
    """
    Dump the interferences group to file, appending to the heirarchical-indexed table.

    Parameters
    ----------
    df : :class:`pandas.DataFrame`
        Dataframe to dump.
    identifier : :class:`str`
        Identifier for the group.
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
    """
    store = path or interferences_datafolder(subfolder="table") / "interferences.h5"

    if isinstance(store, (str, pathlib.Path)):
        store = load_store(
            path=path,
            complevel=complevel,
            complib=complib,
            min_itemsize=_ITEMSIZES,
            **kwargs
        )
    # try to get current index
    if "/table" in store.keys():  #
        current_index = store.select("/table", columns=[]).droplevel("elements").index
    else:
        current_index = pd.DataFrame().index  # empty index
    output = deduplicate(df, charges=charges, multiples=False)
    # take the index from df, and the index from the store and combine them
    dup_multiples = output.index.intersection(
        _find_duplicate_multiples(
            pd.DataFrame(index=output.index.to_list() + current_index.to_list()),
            charges=charges,
        )
    )
    if dup_multiples.size:
        logger.debug(
            "Removing duplicates before dump: {}".format(", ".join(dup_multiples))
        )
        output.drop(dup_multiples, axis="index", inplace=True)
    # create hierarchical indexes
    output = output.set_index(
        pd.MultiIndex.from_product(
            [[identifier], output.index.to_list()], names=["elements", "parts"]
        )
    )
    # convert non-string. non-numerical objects to string
    output = output.astype({"molecule": "str", "components": "str"})
    # append to the existing dataframe
    output.to_hdf(
        store,
        key="table",
        mode="a",
        append=True,
        format="table",
        data_columns=data_columns,
        min_itemsize=_ITEMSIZES,
    )


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
        os.remove(path)  # remove the file
    else:  # keep table keys, set them to empty frames
        logger.debug("Resetting store table: {}/{}".format(path.name, key))
        df = pd.DataFrame(
            index=pd.MultiIndex.from_product([[], []], names=["elements", "parts"]),
            columns=["m_z", "molecule", "components", "mass", "charge", "iso_product",],
        )
        df.to_hdf(
            path,
            key=key,
            format=format,
            mode="w",
            complevel=complevel,
            complib=complib,
            min_itemsize=_ITEMSIZES,
            **kwargs
        )
