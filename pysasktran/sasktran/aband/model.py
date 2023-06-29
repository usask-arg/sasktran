import sasktran as sk
import sasktran.aband.util as util
import numpy as np


class ABandEmission(sk.Emission):
    def __init__(self, altArray_km):
        """
        Photochemical emissions in the oxygen A Band

        Parameters
        ----------
        altArray_km : np.array
            Altitudes to run the model at in km
        """
        self._altArray_km = altArray_km

    def run_model(self, atmosphere, refPoint, sunPos):
        """

        Parameters
        ----------
        atmosphere : sk.Atmosphere
            Atmospheric state to run the emission model at.  Note that the atmosphere MUST contain the keys
            'ozone', 'air', and 'O2'
        refPoint : List
            4 element list of [latitude, longitude, altitude, mjd] which defines the location to run the model at.
            The altitude parameter is generally not used.
        sunPos : np.array
            3 element unit vector of (x, y, z) in geocentric coordinates pointing towards the sun.

        Returns
        -------
        Jemi_nm, wavlenA
            Source function for emission as a two dimensional array (wavelength, altitude), and corresponding
            wavelength array (wavlenA)
        """
        lat = refPoint[0]
        lon = refPoint[1]
        geo = sk.Geodetic()
        geo.from_lat_lon_alt(lat, lon, 0.0)
        upVec = geo.local_up

        # get source function in terms of wavenumber
        Jemi, wavnumA = util.runEmissionModel(refPoint, sunPos, upVec, atmosphere, self._altArray_km)

        # Convert from wavenumber to wavelength
        wavlenA = 1e7 / wavnumA

        Jemi_nm = np.zeros(Jemi.shape)
        for i in range(Jemi.shape[1]):
            Jemi_nm[:, i] = Jemi[:, i] * wavnumA ** 2 / 1e7

        # transpose and rearrange
        wavlenA = wavlenA[::-1]
        Jemi_nm = np.transpose(Jemi_nm)
        Jemi_nm = Jemi_nm[:, ::-1]

        return Jemi_nm, wavlenA

    def skif_object(self, **kwargs):
        altArray_km = self._altArray_km

        engine = kwargs['engine']
        atmosphere = kwargs['atmosphere']

        refPoint = engine.model_parameters['referencepoint']  # [lat, lon, alt, mjd]
        sunPos = engine.model_parameters['sun']  # [x, y, z] geocentric coordinates

        Jemi_nm, wavlenA = self.run_model(atmosphere, refPoint, sunPos)
        table = sk.EmissionTable(altArray_km * 1000, wavlenA, Jemi_nm * 100)  # SASKTRAN requires a cm-m conversion

        return table.skif_object()
