from itertools import groupby, product, combinations_with_replacement
import numpy as np
import pandas as pd
import periodictable as pt
import matplotlib.pyplot as plt
from collections import defaultdict
from copy import deepcopy
from pyrolite.normalisation import ReferenceCompositions

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

    def zero():
        return 0

    atoms = defaultdict(zero)
    for iso in components:
        if isinstance(iso, pt.formulas.Formula):
            try:
                iso = get_first_atom(iso)
            except:
                raise AssertionError
        name = element_name(iso)
        atoms[name] += 1

    molecule = ''
    for n, count in atoms.items():
        if count == 1:
            molecule += '{}'.format(n)
        else:
            molecule += '{}{}'.format(n, count)
    return pt.formula(molecule)


def element_name(el):
    name = ''
    if hasattr(el, 'element'):
        name += '{}'.format(el.element)
        if hasattr(el, 'isotope'):
            name += '[{}]'.format(el.isotope)
    else:
        name += '{}'.format(el.symbol)
    return name


def molecular_abundance(components):
    """Builds a pt.Formula from a list of atom or isotope components."""
    abund = 1.
    for iso in components:
        abund *= iso.abundance
    return abund


def get_first_atom(molecule):
    if isinstance(molecule, pt.core.Element):
        return molecule
    else:
        return list(molecule.atoms.keys())[0]


def build_mz_ratios(formula, charges={1:1., 2:1.}):
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
        molecules = cleanup_isotope_list(isotopes)
        abundances = [i.abundance for i in molecules]
    else:
        components = []
        for a in atoms.keys():
            components += [a for _ in range(atoms[a])]

        components = [[el.add_isotope(i) for i in el.isotopes]
                      for el in components]
        components = [cleanup_isotope_list(_i) for _i in components]
        poss_mol_parts = list(product(*components))
        molecules = [build_molecule(parts) for parts in poss_mol_parts]
        abundances = [molecular_abundance(parts) for parts in poss_mol_parts]

    masses = [m.mass for m in molecules]

    for mol, abund, mass in zip(molecules, abundances, masses):
        for z in charges.keys():
            m_on_z = mass/z
            zmol = pt.formula(deepcopy(mol))
            #setattr(list(zmol.atoms.keys())[-1], 'charge', z)
            mz_ratios[m_on_z] = (zmol, abund * charges[z])
    return mz_ratios


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

def get_molecular_combinations(elements, max_atoms=3):
    """
    Combine a set of elements into molecular combinations up to a maximum
    number of atoms per molecule. Successively adds smaller molecules until
    down to single atoms.
    """
    poss_mol_parts = []
    n = max_atoms
    while n:
        components = combinations_with_replacement(elements,n)
        poss_mol_parts += list(components)
        n -= 1
    molecules = [build_molecule(parts) for parts in poss_mol_parts]

    return molecules


def get_reference_abundance(molecule,
                            norm='Chondrite_PON',
                            unknown_val=1000):
    norm = ReferenceCompositions()[norm] # in ppm
    abund = 10**6 # 100%
    for iso in molecule.atoms.keys():
        el = element_name(iso)
        abund *= getattr(norm, el)
    if not np.isfinite(abund):
        abund = unknown_val # unknown abundance%
    return abund




elements = ['Fe', 'Ti', 'Ni', 'S',  'Co', 'Au', 'Pd', 'Pt', 'Ir', 'Rh', 'Ru',
            'As', 'Se', 'Zn', 'Cu', 'Mn',
            'O', 'He', 'H', 'Ar', 'N']


els = sorted([pt.formula(el) for el in elements],
             key=lambda x: x.mass,
             reverse=True)


molecules = get_molecular_combinations(els, max_atoms=2)
gs = set(map(get_first_atom, molecules))
_mols = [list(filter(lambda x: get_first_atom(x)==g, molecules)) for g in gs]

for g in _mols:
    # estimate some natural relative abundances for the molecules
    rel_abunds =  np.array([get_reference_abundance(mol)
                            for mol in g]).astype(np.float)
    rel_abunds /= np.sum(rel_abunds)
    # get the m/z ratios and relative isotopic abundances for the molecules
    mzdicts = [build_mz_ratios(mol, charges={1:1., 2:0.05}) for mol in g]
    fig, ax = plt.subplots(1, figsize=(10,5))

    _mzs = np.array([mz for mzlist in mzdicts for mz in list(mzlist.keys())])
    min, max = np.round(_mzs.min(),0)-0.5, np.round(_mzs.max(),0)+1.5
    bins = np.arange(min, max, 1)
    #ax.set_yscale('log')
    #ax.set_ylim((np.float(10**10) ,np.float(10**20)))
    series = 0
    bottom = np.zeros(bins.size-1)
    for mol, mz, abund in zip(g, mzdicts, rel_abunds):
        mzs = np.array(list(mz.keys()))
        if np.logical_and(mzs>=np.min(bins), mzs<=np.max(bins)).any():
            n, bins, patches = ax.hist(list(mz.keys()),
                        weights=[mz[m][1]*abund for m in list(mz.keys())],
                        label=str(mol),
                        bins=bins,
                        bottom=bottom,
                        alpha=0.5,)

            bottom += n
            series += 1

    ax.set_ylabel('Est. Relative Intensity')
    ax.set_xlabel('m/z')
    ax.legend(frameon=False, facecolor=None,
              bbox_to_anchor=(1,1),
              bbox_transform=ax.transAxes, #ncol=series//17
              )

build_mz_ratios(pt.formula('CoAr'), charges={1:1., 2:0.05})
build_mz_ratios(pt.formula('Ru'), charges={1:1., 2:0.05})

if __name__ == '__main__':

    from pyrolite.util.pd import test_ser

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

    # populate dict by element

    c = test_ser(index=['Si', 'Ca', 'Mg', 'Ar'])
    c

    isotp = pt.Ca.add_isotope(40)
    isotp = pt.formula(compound='Ca[40]')
    isinstance(isotp, pt.formulas.Formula)
    #print(isotp.atoms)

    sum_of_interferences(isotp, composition=c)
    pt.Zr.isotopes
