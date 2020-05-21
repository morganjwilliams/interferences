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

    try:
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
                tbl_idents = df.index.droplevel("parts")
                df = df.loc[[i for i in identifier if i in tbl_idents]]

            return df
    except KeyError:
        raise KeyError("Key not in HDFStore.")


'''
def lookup_component_subtable(
    identifier, path=None, key="table", window=None, drop_first_level=True, **kwargs
):
    """
    Look up a component-subtable from the store based on an identifier.

    Parameters
    ----------
    identifier : :class:`str`
        Identifier for the subtable.
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
    logger.debug("Attempting lookup for Identifer: {}".format(identifier))
    window = process_window(window)
    name = "/" + key
    with load_store(path, **kwargs) as store:
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
'''


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
    with load_store(path=path, **kwargs) as store:
        store.keys()
        if "/table" in store.keys():
            index = store.select("/table", columns=[]).index
        else:
            index = _get_default_multiindex()  # empty index
    if drop_first_level:
        index = index.droplevel("elements")
    return index


def dump_subtable(
    df,
    ID,
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
    ID : :class:`str`
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
    path = path or interferences_datafolder(subfolder="table") / "interferences.h5"
    current_index = get_store_index(path)

    output = deduplicate(df, charges=charges, multiples=True)
    # take the index from df, and the index from the store and combine them
    tabledelta_duplicates = _find_duplicate_multiples(
        pd.DataFrame(index=output.index.to_list() + current_index.to_list()),
        charges=charges,
    )
    new_duplicates = output.index.intersection(tabledelta_duplicates)

    if new_duplicates.size:
        logger.debug(
            "Removing duplicates before dump: {}".format(", ".join(new_duplicates))
        )
        output.drop(new_duplicates, axis="index", inplace=True)
    # create hierarchical indexes
    output = output.set_index(
        pd.MultiIndex.from_product(
            [[ID], output.index.to_list()], names=["elements", "parts"]
        )
    )
    # convert non-string. non-numerical objects to string
    # append to the existing dataframe
    logger.debug("Dumping {} table to HDF store.".format(ID))
    output.astype({"molecule": "str", "components": "str"}).to_hdf(
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
            columns=["m_z", "molecule", "components", "mass", "charge", "iso_product",],
        )
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
