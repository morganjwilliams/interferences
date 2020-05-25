import unittest
import periodictable as pt
import matplotlib.axes
from interferences.plot.spectra import stemplot, spectra
from interferences.table.build import build_table


class TestSpectraStemPlot(unittest.TestCase):
    def setUp(self):
        self.components = ["H", "He"]
        self.tbl = build_table(self.components)
        self.window = (0.5, 1.5)

    def test_nottable_no_components_raises(self):
        with self.assertRaises(AssertionError):
            ax = stemplot()

    def test_default_table(self):
        ax = stemplot(table=self.tbl)
        self.assertIsInstance(ax, matplotlib.axes.Axes)

    def test_default_components(self):
        ax = stemplot(components=self.components)
        self.assertIsInstance(ax, matplotlib.axes.Axes)

    def test_window(self):
        ax = stemplot(table=self.tbl, window=self.window)
        self.assertIsInstance(ax, matplotlib.axes.Axes)

    def test_no_text_adjust(self):
        ax = stemplot(table=self.tbl, adjust_text=False)
        self.assertIsInstance(ax, matplotlib.axes.Axes)


class TestSpectraSpectraPlot(unittest.TestCase):
    def setUp(self):
        self.components = ["H", "He"]
        self.tbl = build_table(self.components)
        self.window = (0.5, 1.5)

    def test_nottable_no_components_raises(self):
        with self.assertRaises(AssertionError):
            ax = spectra()

    def test_default_table(self):
        ax = spectra(table=self.tbl)
        self.assertIsInstance(ax, matplotlib.axes.Axes)

    def test_default_components(self):
        ax = spectra(components=self.components)
        self.assertIsInstance(ax, matplotlib.axes.Axes)

    def test_window(self):
        ax = spectra(table=self.tbl, window=self.window)
        self.assertIsInstance(ax, matplotlib.axes.Axes)

    def test_no_text_adjust(self):
        ax = spectra(table=self.tbl, adjust_text=False)
        self.assertIsInstance(ax, matplotlib.axes.Axes)

    def test_abb(self):
        ax = spectra(table=self.tbl, abb=1.0)
        self.assertIsInstance(ax, matplotlib.axes.Axes)

    def test_mass_resolution(self):
        ax = spectra(table=self.tbl, mass_resolution=500)
        self.assertIsInstance(ax, matplotlib.axes.Axes)



if __name__ == "__main__":
    unittest.main()
