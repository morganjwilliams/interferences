from pathlib import Path
from pyrolite.util.meta import get_module_datafolder
import logging

logging.getLogger(__name__).addHandler(logging.NullHandler())
logger = logging.getLogger(__name__)


def interferences_datafolder(subfolder=None):
    """
    Returns the path of the interferences data folder.

    Parameters
    -----------
    subfolder : :class:`str`
        Subfolder within the interferences data folder.

    Returns
    -------
    :class:`pathlib.Path`
    """
    return get_module_datafolder(module="interferences", subfolder=subfolder)
