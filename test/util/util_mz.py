import unittest
import numpy as np
import periodictable as pt
from interferences.util.mz import process_window


class TestProcessWindow(unittest.TestCase):
    def test_floats(self):
        window = (0.1, 10.0)
        pw = process_window(window)
        self.assertTrue(all([isinstance(v, (float, int)) for v in pw]))

    def test_isotope_string(self):
        window = ("Ca[40]", 0.05)
        pw = process_window(window)
        self.assertTrue(all([isinstance(v, (float, int)) for v in pw]))

    def test_element_string(self):
        window = ("Ca", 0.05)
        pw = process_window(window)
        self.assertTrue(all([isinstance(v, (float, int)) for v in pw]))

    def test_isotope_object(self):
        window = (pt.Ca.add_isotope(40), 0.05)
        pw = process_window(window)
        self.assertTrue(all([isinstance(v, (float, int)) for v in pw]))

    def test_element_object(self):
        window = (pt.Ca, 0.05)
        pw = process_window(window)
        self.assertTrue(all([isinstance(v, (float, int)) for v in pw]))


if __name__ == "__main__":
    unittest.main()
