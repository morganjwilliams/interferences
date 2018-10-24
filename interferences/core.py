import numpy as np
import pandas as pd
import periodictable as pt
from collections import defaultdict

_MASS = defaultdict(dict)
_ELEMENTS = defaultdict(dict)

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


#pt.formula('{}{}'.format(str(iso), str(iso2)))

iso = pt.formula('1H')

# populate dict by element
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


def abundance_estimiate(composition, formula):

    if isinstance(formula, pt.core.Element):
        return getattr(composition, str(formula.element), 0.)
    else isinstance(formula, pt.core.Formula):
        prod = 1.
        for a in formula.atoms:
            prod *= getattr(composition, str(a.element), 0.)

        return prod

def sum_of_interferences(isotope, composition:pd.Series=None):
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


if __name__ == '__main__':

    from pyrolite.util.pd import test_ser

    c = test_ser(index=['Si', 'Ca', 'Mg', 'Ar'])
    c
    isotp = pt.Ca.add_isotope(40)

    sum_of_interferences(isotp, composition=c)
    pt.Zr.isotopes
