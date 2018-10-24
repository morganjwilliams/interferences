import unittest

from interferences import *


class TestSumInferferences(unittest.TestCase):

    def setUp(self):
        pass

    def test_no_interferences(self):
        sum_of_interferences(pt.H.add_isotope(2))

    def test_positive_interference_sum(self):
        pass

    def test_no_composition(self):
        pass

    def test_composition_provided(self):
        pass


class TestOptimalIsotope(unittest.TestCase):

    def setUp(self):
        pass

    def test_one_isotope(self):
        pass

    def test_multiple_one_optimal(self):
        minimum_interference_isotope(pt.Zr)

    def test_mutiple_optimal(self):
        pass

    def test_no_composition(self):
        pass

    def test_composition_provided(self):
        pass


if __name__ == '__main__':
    unittest.main()
