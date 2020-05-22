import unittest
import pandas as pd
import periodictable as pt
from interferences.table import build_table
from interferences.util.meta import interferences_datafolder


class TestMZAccessor(unittest.TestCase):
    def setUp(self):
        self.dirpath = interferences_datafolder(subfolder="table")
        self.filepath = self.dirpath / "interferences.h5"
        self.elements = ["H", "O"]

    def test_default_build(self):
        df = build_table(self.elements)
        # check the interf interface is present
        self.assertTrue(hasattr(df, "mz"))


if __name__ == "__main__":
    unittest.main()
