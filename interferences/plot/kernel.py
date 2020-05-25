"""
Kernel functions for plotting mass spectra.
"""
import scipy.signal
import numpy as np
from copy import deepcopy
from ..util.log import Handle

logger = Handle(__name__)


def peak_kernel(mass_resolution=1000, image_ratio=0.0, sig_res=1000):
    """
    Kernel function for a peak given a mass resolutiona and level of abbheration.

    Parameters
    ----------
    mass_resolution : :class:`float`
        Mass resolution :math:`\Delta M/M` defined at Full Width Half Maximum (FWHM),
        used to scale the width/mass range of the peak.
    image_ratio : :class:`float`
        The ratio of the size of the image to the limiting slit. Values between 0
        (zero width image) and 1 correspond to scenarios where the full image fits
        in a collector slit. Beyond 1, this corresponds to scenarios where the image
        is larger than the slit, with reduced maximum intensities.
    sig_res : :class:`int`
        Resolution for the signal. Returned signal will have a length equal to three
        times the signal resolution.

    Returns
    -------
    index : :class:`numpy.ndarray`
        Array of indexes corresponding to fractional mass values (:math:`M / M_{peak}`).
    signal : :class:`numpy.ndarray`
        Signal corresponding to fractional intensity values (:math:`I / I_{peak}`).
    ratio : :class:`float`
        Integrated intensity within the FWHM.
    """
    # smoothing from 0 to 1 - full beam fits in slit; above 1 max < 100%
    # assume xvars of -range, +range relative to mass M
    index = np.linspace(
        1 - 3 / (2 * mass_resolution), 1 + 3 / (2 * mass_resolution), 3 * sig_res
    )
    _sig = np.repeat([0.0, 1.0, 0.0], sig_res)  # length of signal is 3x sigres
    if image_ratio:
        win = scipy.signal.triang(int(sig_res * image_ratio))
        sig = scipy.signal.convolve(_sig, win, mode="same") / sum(win)
        ratio = sig[sig_res : 2 * sig_res + 1].sum() / sig.sum()
        return index, sig, ratio
    else:
        sig = _sig
        return index, sig, 1.0


def peak(mz, intensity, kernel=None, mass_resolution=1000, image_ratio=0.0, **kwargs):
    """
    Get arrays corresponding to the intensity vs m/z for a specific peak,
    at a specified mass_resolution.

    Parameters
    ----------
    mass_resolution : :class:`float`
        Mass resolution :math:`\Delta M/M` defined at Full Width Half Maximum (FWHM),
        used to scale the width/mass range of the peak.
    image_ratio : :class:`float`
        The ratio of the size of the image to the limiting slit. Values between 0
        (zero width image) and 1 correspond to scenarios where the full image fits
        in a collector slit. Beyond 1, this corresponds to scenarios where the image
        is larger than the slit, with reduced maximum intensities.
    """
    if kernel is None:
        idx, signal, perc = peak_kernel(
            mass_resolution=mass_resolution, image_ratio=image_ratio, **kwargs
        )
    else:
        idx, signal, perc = deepcopy(kernel)

    idx *= mz
    signal *= intensity
    return idx, signal
