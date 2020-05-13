from collections import defaultdict
import periodictable as pt
from pyrolite.mineral.transform import merge_formulae
from ..util.sorting import get_relative_electronegativity


def _get_isotope(element):
    """
    Parameters
    ----------
    element : :class:`periodictable.core.Element`
        Element or isotope.
    Returns
    -------
    :class:`int`
    """
    try:
        return element.isotope
    except AttributeError:
        return 0


def get_molecular_formula(molecule):
    """
    Construct a name for a molecule.

    Parameters
    -----------
    molecule : :class:`~periodictable.formulas.Formula`

    Returns
    -------
    :class:`str`
    """
    components = list(molecule.atoms.keys())
    components = sorted(
        components, key=lambda x: (get_relative_electronegativity(x), _get_isotope(x))
    )
    name = r"$\mathrm{"  # remove italicized text effect
    for c in components:
        part = ""
        if hasattr(c, "isotope"):
            part += "^{" + "{}".format(c.isotope) + "}"  # superscript isotope
        part += str(c.element)
        count = molecule.atoms[c]
        if count > 1:
            part += "_{" + "{:d}".format(molecule.atoms[c]) + "}"
        name += part
    name += "}$"  # finish TeX formatting
    return name


# see pyrolite.mineral.transform.merge_formulae
def molecule_from_components(components):
    """
    Builds a :class:`~periodictable.formulas.Formula` from a list of atom or
    isotope components.

    Parameters
    ----------
    components : :class:`list`
        Atomic, isotope or molecular components to construct an ionic molecule from.

    Returns
    -------
    :class:`~periodictable.formulas.Formula`

    Todo
    -----
    * Modify to accept consumption of molecular components (e.g. Fe2O3+)
    """
    return merge_formulae(components)
