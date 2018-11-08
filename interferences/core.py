import itertools
import numpy as np
import pandas as pd
import periodictable as pt
import matplotlib.pyplot as plt
from collections import defaultdict

_MASS = defaultdict(dict)
_ELEMENTS = defaultdict(dict)


def cleanup_isotope_list(isotopes):
    """Remove isotopes from a list which have no or zero abundance."""
    for const in [lambda x: hasattr(x, 'abundance'),
                  lambda x: x.abundance > 0]:
        isotopes = [i for i in isotopes if const(i)]
    return isotopes


def build_molecule(components, charge=0):
    """Builds a pt.Formula from a list of atom or isotope components."""
    molecule = ''
    for iso in components:
        molecule += '{}[{}]'.format(iso.element, iso.isotope)
    return pt.formula(molecule)


def molecular_abundance(components):
    """Builds a pt.Formula from a list of atom or isotope components."""
    abund = 1.
    for iso in components:
        abund *= iso.abundance
    return abund


def get_first_atom(molecule):
    return list(molecule.atoms.keys())[0]


def build_mz_ratios(formula, charges={1:1., 2:1.}):
    """
    Calculate the m/z ratios and relative abundance of different isotopic
    compositions from a single molecule. Optionally, higher charges can be
    penalized.
    """
    mz_ratios = {}
    atoms = formula.atoms
    if len(atoms.keys()) < 2:
        el = list(formula.atoms.keys())[0]
        isotopes = [el.add_isotope(i) for i in el.isotopes]
        molecules = cleanup_isotope_list(isotopes)
        abundances = [i.abundance for i in molecules]
    else:
        components = []
        for a in atoms.keys():
            components += [a for _ in range(atoms[a])]

        components = [[el.add_isotope(i) for i in el.isotopes]
                      for el in components]
        components = [cleanup_isotope_list(_i) for _i in components]
        poss_mol_parts = list(itertools.product(*components))
        molecules = [build_molecule(parts) for parts in poss_mol_parts]
        abundances = [molecular_abundance(parts) for parts in poss_mol_parts]

    masses = [m.mass for m in molecules]

    for mol, abund, mass in zip(molecules, abundances, masses):
        for z in charges.keys():
            m_on_z = mass/z
            mz_ratios[m_on_z] = (mol, abund * charges[z])
    return mz_ratios.copy()


def sum_of_interferences(ion, composition:pd.Series=None):
    """Calculate the sum of all interferences for a given isotope."""
    result = 0.

    if composition is None:
        if isotope.isotope in _ELEMENTS[isotope.element]:
            _interf = [v for (k, v) in
                       _ELEMENTS[isotope.element][isotope.isotope].items()]
            if _interf:
                result += np.array(_interf).sum()
    else:
        assert isinstance(composition, pd.Series)
        int_list = _ELEMENTS[isotope.element]
        if isotope.isotope in int_list:
            table = _ELEMENTS[isotope.element][isotope.isotope].items()
            _interf = [v * abundance_estimiate(composition, k)
                       for (k, v) in table]
            if _interf:
                result += np.array(_interf).sum() / \
                          composition[str(isotope.element)]

    return result


def minimum_interference_isotope(element, composition=None):
    """Find the minimum total interference isotope for a given element."""

    isos = list(_ELEMENTS[element].keys())
    _interfs = [sum_of_interferences(element.add_isotope(iso),
                                     composition=composition)
                for iso in isos]

    return isos[_interfs.index(np.min(_interfs))]


def abundance_estimiate(composition, formula):

    if isinstance(formula, pt.core.Element):
        return getattr(composition, str(formula.element), 0.)
    elif isinstance(formula, pt.formulas.Formula):
        prod = 1.
        for a in formula.atoms.keys():
            prod *= getattr(composition, str(a.element), 0.)
        return prod
    elif isinstance(formula, str):
        return abundance_estimiate(composition, pt.formula(compound=formula))
    else:

        print(formula, type(formula))
        raise AssertionError




# populate dict by mass - elements then molecular ions
_els = [i for i in pt.elements if not i.symbol == 'n']
for el1 in _els:

    for z1 in el1.isotopes:
        iso1 = el1.add_isotope(z1)
        abund1 = iso1.abundance
        if abund1:
            _MASS[z1].update({iso1: abund1})

            for el2 in _els:
                for z2 in el2.isotopes:
                    iso2 = el2.add_isotope(z2)
                    abund2 = iso2.abundance
                    if abund2:
                        molecule = '{}[{}]{}[{}]'.format(iso1.element,
                                                         iso1.isotope,
                                                         iso2.element,
                                                         iso2.isotope)
                        molecule = pt.formula(molecule)
                        molecule_abund = abund1*abund2
                        _MASS[z1+z2].update({repr(molecule): molecule_abund})

for el in pt.elements:
    for z in el.isotopes:
        iso = el.add_isotope(z)
        abund = iso.abundance

        if abund:
            # could divide by isotopic composition here /abund
            __interferences = {k: v for (k,v) in _MASS[z].items()
                              if k!=el}
            if __interferences:
                _ELEMENTS[el][z] = __interferences




if __name__ == '__main__':

    from pyrolite.util.pd import test_serf

    from pyrolite.normalisation import ReferenceCompositions
    norm = ReferenceCompositions()['Chondrite_PON']
    isos = [pt.formula('NiS'),
            pt.formula('Pt'), pt.formula('Au'),
            pt.formula('Ru'), pt.formula('Pd'),
            pt.formula('Rh'), pt.formula('Ir')]
    rel_abunds =  [getattr(norm, str(get_first_atom(iso))) for iso in isos]
    mzdicts = [build_mz_ratios(iso, charges={1:1., 2:1.}) for iso in isos]
    fig, ax = plt.subplots(1)
    ax.set_yscale('log')
    for iso, mz, abund in zip(isos, mzdicts, rel_abunds):
        c = ax.hist(list(mz.keys()),
                    weights=[mz[m][1]*abund for m in list(mz.keys())],
                    bins=100,
                    label=str(iso))
    ax.legend(frameon=False, facecolor=None)
    # populate dict by element

    c = test_ser(index=['Si', 'Ca', 'Mg', 'Ar'])
    c

    isotp = pt.Ca.add_isotope(40)
    isotp = pt.formula(compound='Ca[40]')
    isinstance(isotp, pt.formulas.Formula)
    #print(isotp.atoms)

    sum_of_interferences(isotp, composition=c)
    pt.Zr.isotopes
