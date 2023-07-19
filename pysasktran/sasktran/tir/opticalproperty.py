from sasktran.opticalproperty import HITRANChemical
from sasktran.exceptions import wrap_skif_functionfail


class HITRANChemicalTIR(HITRANChemical):
    """
    Calculates the optical absorption and extinction of various atmospheric molecules using the Voigt line-shape and
    the HITRAN spectral line database. The object supports all of the HITRAN species specified in the HITRAN database
    file molparam.txt.

    This is a modified version of the HITRANChemical class for use with the TIR engine. The use_cache property is
    forced to False because the TIR engine uses the internal array-based cross section calculation. This property
    may cause issues if set to True, particularly when calculate_radiance is more than once on the same engine object.

    Examples
    --------
    >>> from sasktran.tir.opticalproperty import HITRANChemicalTIR
    >>> from sasktran import MSIS90
    >>> hitran_co2 = HITRANChemicalTIR('CO2', micro_window_margin=50)
    >>> atmospheric_state = MSIS90()
    >>> hitran_co2.calculate_cross_sections(atmospheric_state, latitude=0, longitude=0, altitude=20000, mjd=54372,\
                                            wavelengths=[10000, 10001])
    CrossSections(wavelengths=array([10000, 10001]), absorption=array([8.23591895e-30, 2.03913420e-28]), scattering=array([0., 0.]), total=array([8.23591895e-30, 2.03913420e-28]))

    Parameters
    ----------
    chemical_name : str
        Chemical abbreviation of the molecule of interest.
    isotope_filter : int, optional
        Allows the HITRAN object to load in just one isotope of the requested molecule. The value set must match one of
        the isotope labels used for the given molecule in the HITRAN database file, molparam.txt. Note that the code
        does not adjust the line strength but uses the line strength value as written in the HITRAN database. This
        means you may have to account for and/or remove the abundance automatically built into the HITRAN database line
        strength values.
    line_tolerance : float, optional
        Allows the user to set the tolerance used to reject weak lines from the current micro-window. The default value
        is 1.0E-09. A larger value will speed up calculation of spectra but may result in choppy spectra at the smaller
        intensities, especially in extinction/absorption spectra which typically follow the log of the cross-section. A
        smaller value will reduce choppiness but increase computational speed. A similar result can be achieved by
        changing property max_line_strength; the choice is really down to the users preference. Only values greater
        than zero are acceptable.
    max_line_strength : float, optional
        Allows the user to manually set the maximum line strength within a micro-window. By default the object will use
        the strongest line in the micro-window. The value is used with the line tolerance to reject weak lines from
        spectral calculations. Reducing the value of the maximum line strength can reduce choppiness in the spectra.
        Its use is similar to property line_tolerance. Only values greater than zero are acceptable. A negative value
        is an error. A value of zero will disable the manual setting and reinstate usage of the default.
    micro_window_margin: double, optional
        Sets the margin of the micro-window in wavenumbers. This margin extends the upper and lower bounds of the
        micro-window once it is loaded into memory. The margin value is used to ensure that cross-section calculations
        near the edge of the micro-window are accurate and have contributions from lines outside the micro-window. The
        default value is 10 wavenumbers. It is the users responsibility to choose a value for the margin that provides
        the necessary accuracy for their application.
    """
    @wrap_skif_functionfail
    def __init__(self, chemical_name: str, line_tolerance=None, max_line_strength=None, isotope_filter=None,
                 micro_window_margin=None):
        self._micro_window_margin = micro_window_margin
        super().__init__(chemical_name=chemical_name, line_tolerance=line_tolerance,
                         max_line_strength=max_line_strength, isotope_filter=isotope_filter, use_cache=False)

    def _update_opticalproperty(self):
        super()._update_opticalproperty()
        if self._micro_window_margin is not None:
            self._iskopticalproperty.SetProperty('setmicrowindowmargin', self._micro_window_margin)
