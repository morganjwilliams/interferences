# interferences

[![interferences Documentation](https://readthedocs.org/projects/interferences/badge/?version=develop)](https://interferences.readthedocs.io/)
[![License: CSIRO Modified BSD/MIT License](https://img.shields.io/badge/License-CSIRO_BSD/MIT_License-blue.svg?style=flat)](https://github.com/morganjwilliams/interferences/blob/master/LICENSE)

Tools for inorganic mass spectra and interference patterns.

## What does it do?

This is an under-development package to facilitate debugging and identification of
issues during method development and analysis of novel samples for geologically-focused
mass spectrometry.

Currently, this package generates simulated combinations of atomic and molecular ions
based on an input set of elements (e.g. `Si, O, Al, Ca, Na, K`), or a specific
composition (e.g. an approximation of the target composition) and some user-controlled
constraints. The relative abundances of these ions are approximated with functions which
modulate estimated ionic abundances using i) the *natural* isotopic abundances of their
constituent ions, ii) user-input constraints on the target compositions, and iii) some
simplified cost functions penalizing higher molecular charges.

The package **does not** account for beam/matrix/source/plasma effects, and likely will
never be able to with any accuracy - it's intended as a practical tool which might help
disentangle issues when analysing difficult samples.

## How should I use it?

The tool is well placed for use when developing methods for new types of samples to
identify potential issues *a-priori* (in an ideal world), but it might also be useful
dig into where issues may have arisen after the fact (to understand *why* the data is
bad).

Where you're finding unknown interferences or mass peaks, the tool can help you identify
them based on relative intensities and mass positions, or give you a set of candidate
interference ions which you can systematically test for.

## Install

Note that this package is in relatively early stages of development, but is planned to
be released via PyPI with v0.1.

If you want the most up to date *development* version, you can install directly from the
GitHub repo. Note that breaking changes occur on this branch, and is not guaranteed to
remain stable (check the [Development and Build Status](#development--build-status)
below). If you still want to try out the most recent bugfixes and yet-to-be-released
features, you can install this version with:

```bash
pip install git+git://github.com/morganjwilliams/interferences.git@develop#egg=interferences
```

## Planned Functionality

* Periodic table of potential issues based on a (default or specified) composition
* Histogram-based spectra (where mass resolution is known)
* Stem-plot based spectra (where mass resolution is not known)
* Simulated spectra given specified resolution (likely with log-scales by default)

## Development & Build Status

[![Formatted with Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/python/black)
[![Code Quality](https://api.codacy.com/project/badge/Grade/fd9912a3faae43bf84a47e3da685d84c)](https://www.codacy.com/manual/morganjwilliams/interferences?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=morganjwilliams/interferences&amp;utm_campaign=Badge_Grade)

|                                                                                  **master**                                                                                  |                                                                                  **develop**                                                                                   |
|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
|                     [![Build Status](https://travis-ci.org/morganjwilliams/interferences.svg?branch=master)](https://travis-ci.org/morganjwilliams/interferences)                      |                      [![Build Status](https://travis-ci.org/morganjwilliams/interferences.svg?branch=develop)](https://travis-ci.org/morganjwilliams/interferences)                      |
| [![Coverage Status](https://coveralls.io/repos/github/morganjwilliams/interferences/badge.svg?branch=master)](https://coveralls.io/github/morganjwilliams/interferences?branch=master) | [![Coverage Status](https://coveralls.io/repos/github/morganjwilliams/interferences/badge.svg?branch=develop)](https://coveralls.io/github/morganjwilliams/interferences?branch=develop) |

**Maintainer**: Morgan Williams (morgan.williams _at_ csiro.au)
