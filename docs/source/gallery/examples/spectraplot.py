"""
Spectra : spectra
====================

If you would like something that looks similar to your peak scans, you can use
:func:`~interferences.plot.spectra.spectra`, which you can also access from your
dataframe using the :class:`~interferences.mz` accessor:
"""
import matplotlib.pyplot as plt
from interferences import build_table
from pyrolite.geochem.ind import REE

# sphinx_gallery_thumbnail_number = 2

########################################################################################
# Here build a table based on some low-mass isotopes, and focus in on the BO+ ion:
#
window = ("B[10]O[16]", 0.05)
df = build_table(["C", "B", "N", "O"], window=window, max_atoms=2)
########################################################################################
# From this table, we can create our plot, limiting the labelling to the
# five peaks with highest estimated intensity. Note we should specify the mass
# resolution for the simulated peaks:
#
ax = df.mz.spectra(window=window, mass_resolution=3000, max_labels=5, figsize=(8, 4))
plt.show()
########################################################################################
# Thesse peaks better show the 'interference' aspect of these ions at relatively low
# mass resolution, but are notably unnaturally square. To simulate some shoulders for
# your peaks (i.e. a non-zero-width image of your source, and a signal limited by a
# collector slit or similar ) you can specify a ratio for the `image_ratio` keyword
# argument. Here we explore the effect this parameter with a few different values
# (0, <1, 1, and >1):
#
fig, ax = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(8, 8))
for a, ratio in zip(ax.flat, [0, 0.2, 1, 1.5]):
    df.mz.spectra(
        ax=a, window=window, mass_resolution=3000, image_ratio=ratio, max_labels=5
    )
    a.annotate(
        "ratio={:.1f}".format(ratio),
        xy=(0.9, 0.9),
        xycoords=a.transAxes,
        ha="right",
        fontsize=12,
    )
plt.show()
