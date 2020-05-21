"""
Placing constraints on mass spectra using compositional information.
"""
import pandas as pd
import periodictable as pt
from pyrolite.geochem.norm import get_reference_composition
from ..util.log import Handle

logger = Handle(__name__)


def constrained_abundance_estimate(composition, formula):
    """
    Get an abundance estimate for a specific molecule constrained by a
    starting composition.

    Parameters
    ----------
    composition : :class:`pandas.Series`
        Composition of the target.
    formula : :class:`~periodictable.formulas.Formula`
        Formula of the molecule.

    Returns
    -------
    :class:`float`
    """

    if isinstance(formula, pt.core.Element):
        return getattr(composition, str(formula.element), 0.0)
    elif isinstance(formula, pt.formulas.Formula):
        prod = 1.0
        for a in formula.atoms.keys():
            prod *= getattr(composition, str(a.element), 0.0)
        return prod
    elif isinstance(formula, str):
        return constrained_abundance_estimate(composition, pt.formula(compound=formula))
    else:
        # print(formula, type(formula))
        raise AssertionError


def get_reference_abundance(molecule, reference="Chondrite_PON", unknown_val=1000):
    """
    Get a reference abundance for a molecule based on a specified reference composition.

    Parameters
    ----------
    molecule : :class:`~periodictable.formulas.Formula`

    reference : :class:`str`
        Reference composition to calculate molecular abundance for.
    """
    norm = get_reference_composition(reference).comp  # in ppm
    abund = 10 ** 6  # 100%
    for iso in molecule.atoms.keys():
        el = get_element_name(iso)
        abund *= getattr(norm, el)
    if not np.isfinite(abund):
        abund = unknown_val  # unknown abundance%
    return abund
