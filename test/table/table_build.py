import unittest
import pandas as pd
import periodictable as pt
from interferences.table import build_table
from interferences.table.molecules import _get_isotope, components_from_index_value
from interferences.util.sorting import get_relative_electronegativity
from interferences.util.meta import interferences_datafolder


class TestBuildTable(unittest.TestCase):
    def setUp(self):
        self.dirpath = interferences_datafolder(subfolder="table")
        self.filepath = self.dirpath / "interferences.h5"
        self.elements = ["H", "O"]

    def test_default_build(self):
        df = build_table(self.elements)

        # check columns
        for column in [
            "m_z",
            "mass",
            "charge",
            "iso_product",
        ]:
            with self.subTest(column=column):
                self.assertIn(column, df.columns)

        # check index
        self.assertIsInstance(df.index, pd.Index)  # not multi-index
        self.assertIsInstance(df.index[0], str)  # index should all be string values
        # check values

    def test_window(self):
        window = (10, 14)
        df = build_table(self.elements, window=window)
        self.assertTrue((df["m_z"] >= window[0]).all())
        self.assertTrue((df["m_z"] <= window[1]).all())

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
            set_final = set(components_from_index_value(row))
            set_start = set(isotopes)
            self.assertTrue(set_final.issubset(set_start))

    def test_isotope_object_components(self):
        # this fails on new builds
        isotopes = [pt.O.add_isotope(18), pt.H.add_isotope(1)]
        df = build_table(isotopes)
        for row in df.index:
            set_final = set(components_from_index_value(row))
            set_start = set([repr(i) for i in isotopes])
            self.assertTrue(set_final.issubset(set_start))

    def test_components_sorted(self):
        df = build_table(self.elements)
        # These are stringified lists
        components0 = components_from_index_value(df.index.values[df.index.size // 2])
        components0_sorted = sorted(
            components0,
            key=lambda x: (get_relative_electronegativity(x), _get_isotope(x)),
        )
        self.assertEqual(components0_sorted, components0)


if __name__ == "__main__":
    unittest.main()
