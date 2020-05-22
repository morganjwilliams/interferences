"""
Building Ion Tables
============================

:mod:`interferences` core funciton is to generate, filter and visualise tables
of isotope-specified small inorganic ionic molecules. This example demonstrates how
to build a small table of ions, and some of the options available. Note that as
:mod:`interferences` is largely built around :mod:`pandas`, you can expect to be
working with :class:`~pandas.DataFrame` objects most of the time.
"""
import pandas as pd
from interferences.table import build_table

########################################################################################
# In the simplest case, where you have a list of elments you want to find viable set of
# ions which they may produce, you can simply specify these elements in a call to
# :func:`~interferences.table.build_table`:
df = build_table(["Ca", "O", "Ar", "H"])
df.info()
########################################################################################
# Note that the table is indexed by the ions themselves, with the respective mass/charge
# ratio found under the `m_z` column. Even with a small set of elements, with the
# combination of isotopes, we've generated a fair number of potential ions:
#
df.index.size
########################################################################################
# As this is probably more ions than you want to consider, let's narrow the foucs
# using a mass window (here for :math:`39 \ge m/z \le 41`):
df = build_table(["Ca", "O", "Ar", "H"], window=(39, 41))
df.index.size
########################################################################################
# While you can use m/z ratios if you know them specifically, it's likely more useful
# to specify an ion mass and a mass-window either side, which you can do as follows:
#
df = build_table(["Ca", "O", "Ar", "H"], window=("Ca[40]", 0.01))
df.index.size
########################################################################################
# By default, :mod:`interferences` builds ions for molecules with up to three atoms
# with ionic charges of either `+1` or `+2` (note the sign of the charge is largely
# irrelevant, given a mass spectrometer will be set up for either positive or negative
# ions). If you wanted to use different parameters to generate a table, you can use the
# `max_atoms` and `charges` keyword arguemnts:
df = build_table(["Ca", "O", "Ar", "H"], charges=[1], max_atoms=2)
df.index.size
########################################################################################
# Also, to save time for futher computation, :mod:`interferences` uses a local
# HDFStore to cache results. You can disable this behaviour and gain some speed
# if you're generating one-off large tables by using :code:`cache_results=False`:
df = build_table(["Ca", "O", "Ar", "H"], cache_results=False)
########################################################################################
# If you're finding that you end up with a table which includes minor ions which
# likely have too low an isotopic abundance to influence results, you can
# use the `threshold` keyword argument to adjust the isotopic abundance threshold for
# isotopes used to build the table. Note that these won't be added to the cached
# reference:
df = build_table(["N", "K"], threshold=0.5)
df.index.size
########################################################################################
# Finally, if you're likely to do some plotting with the ion data, you can specify
# this using the `add_labels` keyword argument, which will add nicely formatted labels
# comapatible with :mod:`matplotlib`:
df = build_table(["Ca", "O", "Ar", "H"], add_labels=True)
