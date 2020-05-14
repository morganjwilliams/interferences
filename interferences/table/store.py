import os
import pandas as pd
from ..util.meta import interferences_datafolder
import logging

logging.getLogger(__name__).addHandler(logging.NullHandler())

_COMPLEVEL = 4
_COMPLIB = "lzo"


# consider using heirarchical dataframes instead of multiple tables?
def dump_element_group(
    df,
    identifier,
    path=None,
    mode="a",
    data_columns=["m_z", "iso_abund_product"],
    complevel=_COMPLEVEL,
    complib=_COMPLIB,
    **kwargs
):
    """
    Dump the interferences group to file, appending to the heirarchical-indexed
    table.

    Parameters
    ----------
    df : :class:`pandas.DataFrame`
        Dataframe to dump.
    identifier : :class:`str`
        Identifier for the group.
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
    store = path or interferences_datafolder(subfolder="table") / "groups.h5"

    # convert non-string. non-numerical objects to string
    # create hierarchical indexed
    output = df.set_index(
        pd.MultiIndex.from_product(
            [[identifier], df.index.to_list()], names=["elements", "parts"]
        )
    )
    # append to the existing dataframe
    if not isinstance(path, pd.HDFStore):
        store = load_groups_store(
            path=path, complevel=complevel, complib=complib, **kwargs
        )
    output.astype({"molecule": "str", "components": "str"}).to_hdf(
        store,
        key="table",
        mode="a",
        append=True,
        format="table",
        data_columns=data_columns,
        min_itemsize={"elements": 15, "parts": 25},
    )


def load_groups_store(path=None, complevel=_COMPLEVEL, complib=_COMPLIB, **kwargs):
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
    path = path or interferences_datafolder(subfolder="table") / "groups.h5"
    store = pd.HDFStore(path, complevel=complevel, complib=complib)
    return store


def get_table(store, key):
    """
    Load a specific table from the HDF store.

    Parameters
    ----------
    store : :class:`pandas.HDFStore`
        Store to access.
    key : :class:`str`
        Key to query.

    Returns
    -------
    :class:`pandas.DataFrame`
    """
    return store.get("/" + key)


def get_element_group(
    store, identifier, key="table",
):
    name = "/" + key
    if name in store.keys():
        df = store.get(name)
        if identifier in df.index.get_level_values("elements"):
            return df.loc[identifier]
        else:
            raise KeyError("Identifier not in index.")
    raise KeyError("Key not in table.")


def reset_group_tables(path=None, remove=True, format="table", **kwargs):
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
    """
    path = path or interferences_datafolder(subfolder="table") / "groups.h5"
    if remove:
        os.remove(path)  # remove the file
    else:  # keep table keys, set them to empty frames
        with load_groups_store(path=path) as store:
            df = pd.DataFrame(
                index=pd.MultiIndex.from_product([[], []], names=["elements", "parts"]),
                columns=[
                    "m/z",
                    "molecule",
                    "components",
                    "mass",
                    "charge",
                    "iso_abund_product",
                ],
            )
            df.to_hdf(
                store,
                key="table",
                format=format,
                mode="w",
                min_itemsize={"elements": 15, "parts": 25},
                **kwargs
            )


def consoliate_groups(
    path=None,
    mode="a",
    data_columns=["m_z", "mass", "iso_abund_product"],
    complevel=_COMPLEVEL,
    complib=_COMPLIB,
):
    """
    Consolidate the group store into a singluar interferences table.

    Parameters
    ----------
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
    """
    output_path = (
        path or interferences_datafolder(subfolder="table") / "interferences.h5"
    )
    store = load_groups_store()
    df = get_table(store, "table").droplevel("elements")
    df.sort_values("m_z", inplace=True)
    # pass off to separate h5 file which could be distributed and easily used
    # for queries etc
    df.to_hdf(
        output_path,
        key="table",
        mode=mode,
        data_columns=data_columns,
        complevel=complevel,
        complib=complib,
    )
    return df
