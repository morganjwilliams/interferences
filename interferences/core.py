from itertools import groupby, product, combinations_with_replacement
import numpy as np
import pandas as pd
import periodictable as pt
from periodictable.core import isisotope, iselement, isatom
import matplotlib.pyplot as plt
from collections import defaultdict
from copy import deepcopy
from pyrolite.geochem.norm import get_reference_composition

_MASS = defaultdict(dict)
_ELEMENTS = defaultdict(dict)


def isotope_abundance_threshold(isotopes, threshold=10e-8):
    """
    Remove isotopes from a list which have no or zero abundance.

    Parameters
    ----------
    isotopes : :class:`list`

    Returns
    -------
    :class:`list`
    """
    constraint_functions = [
        lambda x: hasattr(x, "abundance"),
        lambda x: x.abundance > threshold,
    ]
    for valid in constraint_functions:
        isotopes = [i for i in isotopes if valid(i)]
    return isotopes


def molecule_from_components(components, charge=0):
    """
    Builds a :class:`~periodictable.formulas.Formula` from a list of atom or
    isotope components.

    Parameters
    ----------
    components : :class:`list`
        Atomic, isotope or molecular components to construct an ionic molecule from.
    charge : :class:`int`
        Molecular ionic charge.

    Returns
    -------
    :class:`~periodictable.formulas.Formula`

    Todo
    -----
    * Modify to accept consumption of molecular components (e.g. Fe2O3+)
    """

    def zero():
        return 0

    atoms = defaultdict(zero)
    for iso in components:
        if isinstance(iso, pt.formulas.Formula):
            try:
                iso = get_first_atom(iso)
            except:
                raise AssertionError
        atoms[iso] += 1

    molecule = [(c, n) for n, c in atoms.items()]
    return pt.formula(molecule)


def get_element_name(el):
    """
    Construct a name for an element based on its symbol and isotope.

    Parameters
    -----------
    el : :class:`~periodictable.core.Element`
    """
    name = ""
    if hasattr(el, "element"):
        name += "{}".format(el.element)
        if hasattr(el, "isotope"):
            name += "[{}]".format(el.isotope)
    else:
        name += "{}".format(el.symbol)
    return name


def get_first_atom(molecule):
    """
    Get the first atom in a molecular formula.

    Parameters
    ----------
    molecule : :class:`~periodictable.core.Element` | :class:`~periodictable.formulas.Formula`
        Molecule to check.

    Returns
    -------
    :class:`~periodictable.core.Element`
        Element or isotope.
    """
    if isinstance(molecule, pt.core.Element):
        return molecule
    else:
        return list(molecule.atoms.keys())[0]


def get_molecular_abundance(components):
    """
    Estimates the abundnace of a molecule based on the abundance of the isotopic
    components.

    Returns
    -------
    :class:`float`

    Notes
    ------
    This is essentially a simplistic activity model.
    """
    abund = 1.0
    for iso in components:
        abund *= iso.abundance
    return abund


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

        print(formula, type(formula))
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


def build_mz_ratios(formula, charges={1: 1.0, 2: 1.0}):
    """
    Calculate the m/z ratios and relative abundance of different isotopic
    compositions from a single molecule. Optionally, higher charges can be
    penalized.
    """
    mz_ratios = {}
    if isinstance(formula, pt.core.Element):
        formula = pt.formula(formula)
    atoms = formula.atoms
    if len(atoms.keys()) < 2:
        el = get_first_atom(formula)
        isotopes = [el.add_isotope(i) for i in el.isotopes]
        molecules = isotope_abundance_threshold(isotopes)
        abundances = [i.abundance for i in molecules]
    else:
        components = []
        for a in atoms.keys():
            components += [a for _ in range(atoms[a])]

        components = [[el.add_isotope(i) for i in el.isotopes] for el in components]
        components = [isotope_abundance_threshold(_i) for _i in components]
        poss_mol_parts = list(product(*components))
        molecules = [
            molecule_from_components(set_parts) for set_parts in poss_mol_parts
        ]
        abundances = [molecular_abundance(set_parts) for set_parts in poss_mol_parts]

    masses = [m.mass for m in molecules]

    for mol, abund, mass in zip(molecules, abundances, masses):
        for z in charges.keys():
            m_on_z = mass / z
            zmol = pt.formula(deepcopy(mol))
            # setattr(list(zmol.atoms.keys())[-1], 'charge', z)
            mz_ratios[m_on_z] = (zmol, abund * charges[z])
    return mz_ratios


def sum_of_interferences(ion, composition: pd.Series = None):
    """
    Calculate the sum of all interferences for a given isotope.
    """
    result = 0.0

    if composition is None:
        if ion.isotope in _ELEMENTS[ion.element]:
            _interf = [v for (k, v) in _ELEMENTS[ion.element][ion.isotope].items()]
            if _interf:
                result += np.array(_interf).sum()
    else:
        assert isinstance(composition, pd.Series)
        int_list = _ELEMENTS[ion.element]
        if ion.isotope in int_list:
            table = _ELEMENTS[ion.element][ion.isotope].items()
            _interf = [
                v * constrained_abundance_estimate(composition, k) for (k, v) in table
            ]
            if _interf:
                result += np.array(_interf).sum() / composition[str(ion.element)]

    return result


def minimum_interference_isotope(element, composition=None):
    """
    Find the minimum total interference isotope for a given element.
    """

    isos = list(_ELEMENTS[element].keys())
    _interfs = [
        sum_of_interferences(element.add_isotope(iso), composition=composition)
        for iso in isos
    ]

    return isos[_interfs.index(np.min(_interfs))]


def get_molecular_combinations(elements, max_atoms=3):
    """
    Combine a set of elements into molecular combinations up to a maximum
    number of atoms per molecule. Successively adds smaller molecules until
    down to single atoms.

    Parameters
    ----------
    elements : :class:`list`
        Elements or isotopes to combine into molecules.
    max_atoms : :class:`int`
        Maximum number of atoms per molecule. This limits the number of molecules
        returned to the generally most relevant simple molecules.

    Todo
    ----
    Check that isotopes supplied to this function are propogated
    """
    poss_mol_parts = []
    n = max_atoms
    while n:
        components = combinations_with_replacement(elements, n)
        poss_mol_parts += list(components)
        n -= 1
    molecules = [molecule_from_components(set_parts) for set_parts in poss_mol_parts]

    return molecules
