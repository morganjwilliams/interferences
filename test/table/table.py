import unittest
import pandas as pd
import periodictable as pt
from interferences.table import (
    build_table,
    component_subtable,
    get_elemental_combinations,
    get_isotopic_abund_product,
    isotope_abundance_threshold,
)
from interferences.table.molecules import _get_isotope
from interferences.util.sorting import get_relative_electronegativity
from interferences.util.meta import interferences_datafolder


class TestTable(unittest.TestCase):
    def setUp(self):
        self.dirpath = interferences_datafolder(subfolder="table")
        self.filepath = self.dirpath / "interferences.h5"
        self.elements = ["H", "O"]

    def test_default_build(self):
        df = build_table(self.elements)

        # check columns
        for column in [
            "m_z",
            "molecule",
            "components",
            "mass",
            "charge",
            "iso_abund_product",
        ]:
            with self.subTest(column=column):
                self.assertIn(column, df.columns)

        # check index
        self.assertIsInstance(df.index, pd.Index)  # not multi-index
        self.assertIsInstance(df.index[0], str)  # index should all be string values
        # check values

    def test_add_labels(self):
        for add_labels in [True, False]:
            with self.subTest(add_labels=add_labels):
                df = build_table(self.elements, add_labels=True)

    def test_element_object_components(self):
        df = build_table([getattr(pt, el) for el in self.elements])

    def test_isotope_string_components(self):
        isotopes = [repr(pt.O.add_isotope(18)), repr(pt.H.add_isotope(1))]
        df = build_table(isotopes)
        for row in df.index:
            components = [
                el.strip() for el in df.loc[row, "components"][1:-1].split(",")
            ]
            set_final = set(components)
            set_start = set(isotopes)
            self.assertTrue(set_final.issubset(set_start))

    def test_isotope_object_components(self):
        isotopes = [pt.O.add_isotope(18), pt.H.add_isotope(1)]
        df = build_table(isotopes)
        for row in df.index:
            components = [
                el.strip() for el in df.loc[row, "components"][1:-1].split(",")
            ]
            set_final = set(components)
            set_start = set([repr(i) for i in isotopes])
            self.assertTrue(set_final.issubset(set_start))

    def test_components_sorted(self):
        df = build_table(self.elements)
        # These are stringified lists
        components0 = [
            el.strip()
            for el in df["components"].values[df.index.size // 2][1:-1].split(",")
        ]
        components0_sorted = sorted(
            components0,
            key=lambda x: (get_relative_electronegativity(x), _get_isotope(x)),
        )
        self.assertEqual(components0_sorted, components0)


if __name__ == "__main__":
    unittest.main()
