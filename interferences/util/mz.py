import periodictable as pt


def process_window(window):
    """
    Process the two allowable verions of a mass window (element/isotope and width, or
    a low- and high-mass).

    Parameters
    -----------
    window : :class:`tuple`
        Window parameters to process.

    Returns
    --------
    :class:`tuple`
    """
    if window is None:
        return window
    elif isinstance(window[0], (str, pt.core.Element, pt.core.Isotope)):
        peak, width = window
        m_z = pt.formula(peak).mass
        return (m_z - width / 2, m_z + width / 2)
    else:
        return tuple(sorted(list(window)))
