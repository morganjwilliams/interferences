import os
import pandas as pd
from ..util.meta import interferences_datafolder
import logging

logging.getLogger(__name__).addHandler(logging.NullHandler())

_COMPLEVEL = 4
_COMPLIB = "lzo"


def dump_table(
    df,
    key,
    path=None,
    mode="a",
    data_columns=["m_z", "mass", "iso_abund_product"],
    complevel=_COMPLEVEL,
    complib=_COMPLIB,
    **kwargs
):
    """
    Dump the interferences table to file.

    Parameters
    ----------
    df : :class:`pandas.DataFrame`
        Dataframe to dump.
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
    # convert non-string. non-numerical objects to string
    df.astype({"molecule": "str", "components": "str"}).to_hdf(
        path,
        key=key,
        mode=mode,
        format="table",
        data_columns=data_columns,
        complevel=complevel,
        complib=complib,
        **kwargs
    )


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


def reset_tables(path=None, remove=True, format="table", **kwargs):
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
    path = path or interferences_datafolder(subfolder="table") / "interferences.h5"
    if remove:
        os.remove(path)  # remove the file
    else:  # keep table keys, set them to empty frames
        with load_store() as store:
            for key in store.keys():
                df = pd.DataFrame(
                    index=pd.MultiIndex.from_product(
                        [[], []], names=["elements", "components"]
                    ),
                    columns=[
                        "m/z",
                        "molecule",
                        "components",
                        "mass",
                        "charge",
                        "iso_abund_product",
                    ],
                )
                df.to_hdf(store, key=key, format=format, mode="w", **kwargs)
