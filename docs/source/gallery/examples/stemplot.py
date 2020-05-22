"""
Spectra : stemplot
====================

Once you have a table, one way to visualise how the ion peaks are distributed is to use
:func:`~interferences.plot.spectra.stemplot`, which you can also access from your
dataframe using the :class:`~interferences.mz` accessor:
"""
import matplotlib.pyplot as plt
from interferences import build_table
from pyrolite.geochem.ind import REE

########################################################################################
# Let's build a table to play with first, focusing on potential interferences for
# thulium (Tm), which has only one stable isotope (:math:`\mathrm{^{169}Tm}`),
# and is susceptible to interferences, especially for quadrupole ICP-MS:
#
window = ("Tm[169]", 0.1)
df = build_table(REE() + ["O", "N", "H"], window=window, max_atoms=2)
########################################################################################
# From this table, we can create our plot, here limiting the labelling to the
# five peaks with highest estimated intensity:
#
ax = df.mz.stemplot(window=window, max_labels=5, figsize=(8, 4))
plt.show()
########################################################################################
# While the production of the doubly-charged double-REE ions is likely less significant
# than shown here (no penalisation for higher charges/larger molecules is included
# in generating these spectra), we can see that :math:`\mathrm{^{153}Eu^{16}O}` could
# be a potential interference issue if the conditions are relatively oxidised,
# and if there's sufficient hydrogen, :math:`\mathrm{^{168}Er^{1}H}` may similarly
# contribute to problems.
#
# Notably, there's a number of other potential ions in vicinity of
# :math:`\mathrm{^{169}Tm}`. However, most of these are doubly-charged double-REE ions.
# Given the highly-correlated nature of the REE, these may not pose as significant
# issues for standardisation as the hydride and oxide ions.
#
