"""
Functions to threshold, combine and estimate intensities of elements and isotopes
based on their abundances.
"""
from ..util.log import Handle

logger = Handle(__name__)


def isotope_abundance_threshold(isotopes, threshold=None):
    """
    Remove isotopes from a list which have no or zero abundance.

    Parameters
    ----------
    isotopes : :class:`list`
        List of isotopes to filter.
    threshold : :class:`float`
        Minimum isotope abundance for inclusion.

    Returns
    -------
    :class:`list`
    """

    threshold = threshold or 10e-8
    constraint_functions = [
        lambda x: hasattr(x, "abundance"),
        lambda x: x.abundance >= threshold,
    ]
    for valid in constraint_functions:
        isotopes = [i for i in isotopes if valid(i)]
    return isotopes


def get_isotopic_abundance_product(components):
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
