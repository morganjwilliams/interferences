"""
Spectra : stemplot
====================

Once you have a table, one way to visualise how the ion peaks are distributed is to use
:func:`~interferences.plot.spectra.stemplot`, which you can also access from your
dataframe using the :class:`~interferences.mz` accessor:
"""
import matplotlib.pyplot as plt
from interferences import build_table

df = build_table(["Ca", "O", "Ar", "H"], window=("Ar[40]", 0.01))
ax = df.mz.stemplot()
plt.show()
