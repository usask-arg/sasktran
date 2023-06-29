import numpy as np
import numpy.matlib as ml
import math
import time
import logging
import sasktran.aband.load_data as load_data


def mjdToYYYYMMDD(mjd):
    """
    Converts mjd style date to "yyyymmdd" format.

    Parameters
    ----------
    mjd: float
        In MJD format.

    Returns
    -------
    date: float
        In "yyyymmdd" format.
    """
    mjd = int(mjd)
    jd = mjd + 2400001

    F, I = math.modf(jd)
    I = int(I)
    A = math.trunc((I - 1867216.25) / 36524.25)

    if I > 2299160:
        B = I + 1 + A - math.trunc(A / 4.)
    else:
        B = I

    C = B + 1524
    D = math.trunc((C - 122.1) / 365.25)
    E = math.trunc(365.25 * D)
    G = math.trunc((C - E) / 30.6001)
    day = C - E + F - math.trunc(30.6001 * G)

    if G < 13.5:
        month = G - 1
    else:
        month = G - 13

    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715

    year = year * 1e4
    month = month * 1e2
    date = year + month + day

    return date


def getHITRANdata(T):
    """
    Loads and parses HITRAN data for O2 to provide transition wavenumbers and line strengths for standard temperature Ts
    in the A and B bands, and calculates the temperature-dependent line strengths based on input temperature.

    Parameters
    ----------
    T: np.ndarray
        Temperature profile (K), same length as altitude profile

    Returns
    -------
    nuTransA: np.ndarray
        Transition wavenumbers (cm-1) for O2 in the A-band

    nuTransB: np.ndarray
        Transition wavenumbers (cm-1) for O2 in the B-band

    Sj_TA: np.ndarray
        Temperature-dependent line strengths (cm-1/molec/cm-2) corresponding to transition wavenumbers in the A-band

    Sj_TB: np.ndarray
        Temperature-dependent line strengths (cm-1/molec/cm-2) corresponding to transition wavenumbers in the B-band
    """
    # -------------------------------------------------------------------------------------------
    # Define some constants
    # -------------------------------------------------------------------------------------------
    h = 6.6256e-27  # Planck's constant [cm2*kg]
    c = 2.9979e10  # speed of light [cm/s]
    k = 1.3806e-16  # Boltzmann constant [cm2*kg/s2/K]
    hck = h * c / k  # [cm*K]
    Ts = 296  # standard temperature [K]

    # -------------------------------------------------------------------------------------------
    # Parse some data
    # -------------------------------------------------------------------------------------------
    hitran_file = load_data.hitran_o2()
    O2dat = np.genfromtxt(hitran_file,
                          delimiter=[2, 1, 12, 10, 10, 5, 5, 10, 4, 8, 8, 7, 8, 7, 15, 2, 3, 1, 3, 6, 6, 12, 1, 7, 7],
                          dtype=np.str)

    O2datAB = O2dat[5672:6341, :]  # HARDCODE = rows containing A-B bands

    nuTrans = O2datAB[:, 2]  # transition wavenumbers (line positions)   [cm-1]
    Sj_Ts = O2datAB[:, 3]  # line intensity @ Ts = 296K                [cm-1/molec/cm-2]
    E0j = O2datAB[:, 7]  # lower state energy                        [cm-1]

    # convert to floats or strip whitespace from strings
    nuTrans = nuTrans.astype(float)
    Sj_Ts = Sj_Ts.astype(float)
    E0j = E0j.astype(float)

    # -------------------------------------------------------------------------------------------
    # Calculate line strengths (Sj_T)
    # -------------------------------------------------------------------------------------------
    rows = nuTrans.shape[0]
    cols = T.shape[0]
    Sj_T = np.zeros((rows, cols))
    for i in range(cols):
        Sj_T[:, i] = Sj_Ts * Ts / T[i] * np.exp(hck * E0j * (T[i] - Ts) / (T[i] * Ts))

    # HARDCODE - separate out the A and B portions
    nuTransA = nuTrans[0:430]
    nuTransB = nuTrans[431:-1]
    Sj_TA = Sj_T[0:430, :]
    Sj_TB = Sj_T[431:-1, :]

    return nuTransA, nuTransB, Sj_TA, Sj_TB


def getSolFluxData(date):
    """
    Uses a date to select relevant incident solar flux values for the desired bandwidths (Lyman-a, Schumann-Runge,
    Hartley, A-band, B-band). Parses and converts the raw data into appropriate flux funits and pairs up the output with
    the related wavelength array.

    Parameters
    ----------
    date: float
        In yyyymmdd format.

    Returns
    -------
    datLya: np.ndarray
        Column 1 = wavelength for Lyman-a line (nm)
        Column 2 = incident solar flux data (photons/s/cm2/nm)

    datSRC: np.ndarray
        Column 1 = wavelengths for Schumann-Runge continuum (nm)
        Column 2 = incident solar flux data (photons/s/cm2/nm)

    datHar: np.ndarray
        Column 1 = wavelengths for Hartley bands (nm)
        Column 2 = incident solar flux data (photons/s/cm2/nm)

    datA: np.ndarray
        Column 1 = wavelength for A-band (nm)
        Column 2 = incident solar flux data (photons/s/cm2/nm)

    datB: np.ndarray
        Column 1 = wavelength for B-band (nm)
        Column 2 = incident solar flux data (photons/s/cm2/nm)
    """

    # -------------------------------------------------------------------------------------------
    # Load data from txt files -- SORCE data starts from date May 15, 2003 (mjd=52775)
    # -------------------------------------------------------------------------------------------

    rawFUV, rawMUV, rawSIM = load_data.incoming_flux(date)

    # -------------------------------------------------------------------------------------------
    # Set up constants
    # -------------------------------------------------------------------------------------------
    h = 6.6256e-34  # Planck's constant [J*s]
    c = 2.9979e8  # speed of light [m/s]

    # -------------------------------------------------------------------------------------------
    # Convert flux data from [W/m2/nm] to [phot/s/cm2/nm]
    # -------------------------------------------------------------------------------------------
    datFUV = np.zeros((rawFUV.shape[0], 2))
    datFUV[:, 0] = rawFUV[:, 2]  # wavelength [nm]
    datFUV[:, 1] = (1e-13) / h / c * rawFUV[:, 2] * rawFUV[:, 6]  # flux [phot/s/cm2/nm]

    datMUV = np.zeros((rawMUV.shape[0], 2))
    datMUV[:, 0] = rawMUV[:, 2]
    datMUV[:, 1] = (1e-13) / h / c * rawMUV[:, 2] * rawMUV[:, 6]

    datSIM = np.zeros((rawSIM.shape[0], 2))
    datSIM[:, 0] = rawSIM[:, 2]
    datSIM[:, 1] = (1e-13) / h / c * rawSIM[:, 2] * rawSIM[:, 6]

    # -------------------------------------------------------------------------------------------
    # Sort raw data into spectral regions of interest
    # -------------------------------------------------------------------------------------------

    # Lyman-alpha hydrogen line (121.6 nm, 82237 cm-1)
    ix = np.where(datFUV[:, 0] == 121)
    datLya = np.zeros((ix[0].shape[0], 2))
    datLya[:, 0] = 121
    datLya[:, 1] = datFUV[ix[0], 1]

    # Schumann-Runge continuum (130-175 nm, 57143-76923 cm-1)
    ix = np.where(datFUV[:, 0] >= 130)
    datSRC_ = np.zeros((ix[0].shape[0], 2))
    datSRC_[:, 0] = datFUV[ix[0], 0]
    datSRC_[:, 1] = datFUV[ix[0], 1]
    ix = np.where(datSRC_[:, 0] <= 175)
    datSRC = np.zeros((ix[0].shape[0], 2))
    datSRC[:, 0] = datSRC_[ix[0], 0]
    datSRC[:, 1] = datSRC_[ix[0], 1]

    # Hartley bands (198-309 nm, 32362-50505 cm-1)
    ix = np.where(datMUV[:, 0] >= 198)
    datHar_ = np.zeros((ix[0].shape[0], 2))
    datHar_[:, 0] = datMUV[ix[0], 0]
    datHar_[:, 1] = datMUV[ix[0], 1]
    ix = np.where(datHar_[:, 0] <= 309)
    datHar = np.zeros((ix[0].shape[0], 2))
    datHar[:, 0] = datHar_[ix[0], 0]
    datHar[:, 1] = datHar_[ix[0], 1]

    # A-band central wavelength (762 nm, 13123 cm-1)
    ix = np.where(datSIM[:, 0] >= 761)
    datA_ = np.zeros((ix[0].shape[0], 2))
    datA_[:, 0] = datSIM[ix[0], 0]
    datA_[:, 1] = datSIM[ix[0], 1]
    ix = np.where(datA_[:, 0] <= 763)
    datA = np.zeros((ix[0].shape[0], 2))
    datA[:, 0] = 762  # A band wavelength
    datA[:, 1] = datA_[ix[0], 1]  # A band flux

    # B-band central wavelength (688 nm, 14535 cm-1)
    ix = np.where(datSIM[:, 0] >= 687)
    datB_ = np.zeros((ix[0].shape[0], 2))
    datB_[:, 0] = datSIM[ix[0], 0]
    datB_[:, 1] = datSIM[ix[0], 1]
    ix = np.where(datB_[:, 0] <= 689)
    datB = np.zeros((ix[0].shape[0], 2))
    datB[:, 0] = 688  # B band wavelength
    datB[:, 1] = datB_[ix[0], 1]  # B band flux

    return datLya, datSRC, datHar, datA, datB


def DopplerBroadening(T, nuTrans, array_j, nuArr):
    """
    Performs Doppler broadening on an input array, originally dependent on transition wavenumbers. After being Doppler
    broadened, the new spectrum exhibits Gaussian spikes and is dependent on an evenly-spaced array of chosen
    wavenumbers.

    Parameters
    ----------
    T: float
        Single temperature value (K)

    nuTrans: np.ndarray
        Transition wavenumbers (cm-1) for the spectrum

    array_j: np.ndarray
        Original array, one value per line position (xx)

    nuArr: np.ndarray
        Evenly-spaced high resolution wavenumber array (cm-1)

    Returns
    -------
    array_nu: np.ndarray
        Final array, dependent on smooth nuArr (xx/cm-1)
    """
    # -------------------------------------------------------------------------------------------
    # Define some constants
    # -------------------------------------------------------------------------------------------
    c = 2.9979e10  # speed of light [cm/s]
    k = 1.3806e-16  # Boltzmann's constant [cm2*kg/s2*K]
    M_o2 = 31.9988  # mass of O2 [g/mol]
    m_o2 = M_o2 / 6.0221413e23  # mass of O2 [g/molec]

    nu_resolution = nuArr[1] - nuArr[0]

    central_indicies = np.searchsorted(nuArr, nuTrans)

    array_nu = np.zeros_like(nuArr)
    for i in range(nuTrans.shape[0]):
        if array_j[i] / np.max(array_j) < 0:
            continue
        alphaj = nuTrans[i] / c * np.sqrt(2 * k * T / m_o2)
        num_include = alphaj / nu_resolution * 5

        min_index = int(np.max([0, central_indicies[i] - num_include]))
        max_index = int(np.min([len(nuArr), central_indicies[i] + num_include]))

        array_nu[min_index:max_index] += array_j[i] * 1 / alphaj / np.sqrt(np.pi) * np.exp(-((nuArr[min_index:max_index] - nuTrans[i]) / alphaj) ** 2)

    return array_nu, None


def calcAttSolFlux(sza, alt, dens, absXsec, incSolFlx):
    """
    Determines the solar flux spectrum due to attenuation through the atmopshere for each value in the altitude array.
    Uses the Lambert-Beer law and numerical integration for one species across a predefined bandwidth. The incident
    solar flux input must match the unit dependence of absorption x-section (ie. wavenumber vs wavelength). Also
    rearranges the absXsec variable to be in a matrix with the right proportions for subsequent calculations.

    Parameters
    ----------
    sza: float
        Solar zenith angle (degrees)

    alt: np.ndarray
        Altitude profile (km)

    dens: np.ndarray
        Density profile of attenuating species (cm-3)

    absXsec: np.ndarray
        Absorption cross-section of attenuating species in bandwidth of interest
        Depends on wavenumber (cm2/molec/cm-1) or wavelength (cm2/molec/nm)

    incSolFlx: np.ndarray
        Incident solar flux in bandwidth of interest
        Depends on wavenumber (photons/s/cm2/cm-1) or wavelength (photons/s/cm2/cm-1)
        Spectral dependence must match absXsec

    Returns
    -------
    attSolFlux: np.ndarray
        Attenuated solar flux in bandwidth of interest
        Depends on wavenumber (photons/s/cm2/cm-1) or wavelength (photons/s/cm2/nm)

    absXsec_: np.ndarray
        Absorption cross-section matrix properly configured for subsequent calculations.
    """
    # -------------------------------------------------------------------------------------------
    # Restructure vector dimensions
    # -------------------------------------------------------------------------------------------
    # Put density into row vector (1 x alt)
    dens_ = np.zeros((1, alt.shape[0]))
    dens_[0, :] = dens

    # Put absorption X-section into column vectors and repeat into 2d array (alt x spectrum)
    if absXsec.ndim == 1:  # only for non-alt dependent absXsec (Lya, SRC, Har)
        absXsec_ = np.zeros((1, absXsec.shape[0]))  # make row array
        absXsec_[0, :] = absXsec
        absXsec_ = ml.repmat(absXsec_, alt.shape[0], 1)  # make 2d array with alt=rows and spec=cols
    else:
        absXsec_ = np.transpose(absXsec)

    # -------------------------------------------------------------------------------------------
    # Calculate attenuated solar flux
    # -------------------------------------------------------------------------------------------
    # Find the average step size of altitude array & convert units for dz & sza
    dz = np.mean(np.diff(alt)) * 1e5  # [km] to [cm]
    sza = np.deg2rad(sza)  # [degrees] to [radians]

    # Perform the integration
    attSolFlx = np.ones(absXsec_.shape) * incSolFlx

    if absXsec_.shape[1] < 2:  # for Lya (only one flux value/one column)
        for i in range(alt.shape[0]):
            int = np.zeros(absXsec_.shape)  # integrate from current altitude to top
            tau = np.zeros(absXsec_.shape)
            int[i, 0] = np.dot(dens_[0, i:alt.shape[0]], absXsec_[i:alt.shape[0], 0]) * dz
            tau[i, 0] = np.exp(- (1 / np.cos(sza)) * int[i, 0])
            attSolFlx[i, 0] *= tau[i, 0]

    else:  # for SRC, Har, A, B (multiple wavelengths/multiple columns)
        nonzero = np.nonzero(absXsec_[0, :])[0]
        for i in range(alt.shape[0]):
            int = np.dot(dens_[0, i:alt.shape[0]], absXsec_[i:alt.shape[0], nonzero]) * dz
            attSolFlx[i, nonzero] *= np.exp(- (1 / np.cos(sza)) * int)

    return attSolFlx, absXsec_


def calc_gAndJfactors(wvlnSRC, wvlnHar, wavnumA, wavnumB, altArray_km,
                      attSolFlx_Lya, attSolFlx_SRC, attSolFlx_Har, attSolFlx_A, attSolFlx_B,
                      absXsecLya, absXsecSRC, absXsecHar, absXsecA, absXsecB):
    """
    Calculates the photochemical reaction rates for photon absorption in the A and B-bands (g factors) and the rate
    coefficients for photolysis of O2 in the Lyman-a line and Schumann-Runge continuum, and of O3 in the Hartley bands.
    Absorption reaction rates gA and gB are wavenumber dependent (cm-1), and photolysis coefficients JL, JS, and JH are
    wavelength dependent.

    Parameters
    ----------
    wvlnSRC: np.ndarray
        Wavelength array (nm) in Schumann-Runge continuum

    wvlnHar: np.ndarray
        Wavelength array (nm) in Hartley bands

    wavnumA: np.ndarray
        Wavenumber array (cm-1) in A-band

    wavnumB: np.ndarray
        Wavenumber array (cm-1) in B-band

    altArray_km: np.ndarray
        Altitude profile (km)

    attSolFlx_Lya: np.ndarray
        Attenuated solar flux in the Lyman-a line (photons/s/cm2/nm)

    attSolFlx_SRC: np.ndarray
        Attenuated solar flux in the Schumann-Runge continuum (photons/s/cm2/nm)

    attSolFlx_Har: np.ndarray
        Attenuated solar flux in the Hartley bands (photons/s/cm2/nm)

    attSolFlx_A: np.ndarray
        Attenuated solar flux in the A-band (photons/s/cm2/cm-1)

    attSolFlx_B: np.ndarray
        Attenuated solar flux in the B-band (photons/s/cm2/cm-1)

    absXsec_Lya: np.ndarray
        Absorption cross-section of O2 in Lyman-a line (cm2/molec/nm)

    absXsec_SRC: np.ndarray
        Absorption cross-section of O2 in Schumann-Runge continuum (cm2/molec/nm)

    absXsec_Har: np.ndarray
        Absorption cross-section of O3 in Hartley bands (cm2/molec/nm)

    absXsec_A: np.ndarray
        Absorption cross-section of O2 in A-band (cm2/molec/cm-1)

    absXsec_B: np.ndarray
        Absorption cross-section of O2 in B-band (cm2/molec/cm-1)

    Returns
    -------
    gA: np.ndarray
        Photochemical reaction rate profile for O2 photon absorption across the A-band (s-1)

    gB: np.ndarray
        Photochemical reaction rate profile for O2 photon absorption across the B-band (s-1/km)

    JL: np.ndarray
        Photolysis coefficient profile for O2 in the Lyman-a line (s-1)

    JS: np.ndarray
        Photolysis coefficient profile for O2 across the Schumann-Runge continuum (s-1)

    JH: np.ndarray
        Photochemical reaction rate for O3 photon absorption across the Hartley bands (s-1)
    """
    # -------------------------------------------------------------------------------------------
    # Find average step size of arrays
    # -------------------------------------------------------------------------------------------
    dwvlnSRC = np.mean(np.diff(wvlnSRC))  # should be 1 nm
    dwvlnHar = np.mean(np.diff(wvlnHar))
    dnuA = np.mean(np.diff(wavnumA))  # should be 0.0001 cm-1
    dnuB = np.mean(np.diff(wavnumB))

    # -------------------------------------------------------------------------------------------
    # Define constants and coefficients
    # -------------------------------------------------------------------------------------------
    # Quantum yield data
    qyLya = 0.48
    qySRC = 1
    x = np.genfromtxt(load_data.xs_har())
    qyHar = np.zeros((x.shape[0], 1))
    qyHar[:, 0] = x[:, 2]  # make column vector of quantum yield
    qyHar = np.transpose(ml.repmat(qyHar, 1, altArray_km.shape[0]))  # make 2d array to match attSolFlx_Har shape

    # Products of quantum yield & solar flux for ease of J calculations
    qyFluxLya = qyLya * attSolFlx_Lya
    qyFluxSRC = qySRC * attSolFlx_SRC
    qyFluxHar = qyHar * attSolFlx_Har

    # -------------------------------------------------------------------------------------------
    # Calculate g and J factor profiles
    # -------------------------------------------------------------------------------------------
    # Initialize arrays
    gA = np.zeros(altArray_km.shape[0])
    gB = np.zeros(altArray_km.shape[0])
    JL = np.zeros(altArray_km.shape[0])
    JS = np.zeros(altArray_km.shape[0])
    JH = np.zeros(altArray_km.shape[0])

    # Perform integration
    for i in range(altArray_km.shape[0]):
        # Inner product array multiplcation (dot product) multiplied by step size approximates the integral
        gA[i] = np.dot(attSolFlx_A[i, :], absXsecA[i, :]) * dnuA  # g factor data in wavenumber units
        gB[i] = np.dot(attSolFlx_B[i, :], absXsecB[i, :]) * dnuB

        JL[i] = qyFluxLya[i, 0] * absXsecLya[i, 0]  # J factor data in wavelength units
        JS[i] = np.dot(qyFluxSRC[i, :], absXsecSRC[i, :]) * dwvlnSRC
        JH[i] = np.dot(qyFluxHar[i, :], absXsecHar[i, :]) * dwvlnHar

    return gA, gB, JL, JS, JH


def runPhotochemModel(T, altArray_km, gA, gB, J2, J3, densO2, densO3, densN2, densO):
    """
    Runs the full photochemical model O2 excited to its second electronic state. Takes in density profiles, temperature
    profile, and reaction rates, and produces the density profile for the excited O2A due the four contributory
    production mechanisms, separated by O2 and O3 sourced photolysis.

    Parameters
    ----------
    T: np.ndarray
        Temperature profile (K)

    altArray_km: np.ndarray
        Altitude profile (km)

    gA: np.ndarray
        Photochemical reaction rate profile for O2 photon absorption across the A-band (s-1)

    gB: np.ndarray
        Photochemical reaction rate profile for O2 photon absorption across the B-band (s-1/km)

    J2: np.ndarray
        Photolysis coefficient profile for O2 in the Lyman-a line and Schumann-Runge continuum (s-1)

    J3: np.ndarray
        Photochemical reaction rate for O3 photon absorption across the Hartley bands (s-1)

    densO2: np.ndarray
        Density profile of O2 (cm-3)

    densO3: np.ndarray
        Density profile of O3 (cm-3)

    densN2: np.ndarray
        Density profile of N2 (cm-3)

    densO: np.ndarray
        Density profile of O (cm-3)

    Returns
    -------
    densO2Apros: np.ndarray
        Density profiles of excited O2A (cm-3) separated by production mechanisms
        Column 1 = Absorption in the A-band
        Column 2 = Absorption in he B-band and electronic quenching
        Column 3 = Excitation by collision with O1D (produced by photolysis of O2)
        Column 4 = Excitation by collision with O1D (produced by photolysis of O3)
        Column 5 = Two-step Barth process
    """
    # -------------------------------------------------------------------------------------------
    # Define reaction rate coefficients (see Table 3.3)
    # -------------------------------------------------------------------------------------------
    # (g and J factors are calculated in the parent emission model and passed into this function)

    #  Photon absorption in B-band
    A771 = 0.070  # [s-1]
    k3B = 3.0e-10  # [cm3/s]
    k0B = 4.5e-12  # [cm3/s]
    k1B = 4.2e-11 * np.exp(-312 / T)  # [cm3/s]
    k2B = 5.0e-13  # [cm3/s]

    # Collisional excitation due to O1D
    A1D = 6.81e-3  # [s-1]
    k1 = 3.3e-11 * np.exp(55 / T)  # [cm3/s]
    k2 = 2.15e-11 * np.exp(110 / T)  # [cm3/s]
    phi = 0.95  # []

    # Two-step Barth process
    k1star = 4.7e-33 * (300 / T) ** 2  # [cm6/s]
    Co2 = 6.6  # []
    Co = 19  # []

    # Loss mechanisms of O2A
    A1Sig = 0.085  # [s-1]
    k1A = 1.8e-15 * np.exp(45 / T)  # [cm3/s]
    k2A = 3.5e-11 * np.exp(-135 / T)  # [cm3/s]
    k3A = 3.9e-17  # [cm3/s]
    k4A = 8.0e-14  # [cm3/s]

    # -------------------------------------------------------------------------------------------
    # Production mechanisms of A-band O2
    # -------------------------------------------------------------------------------------------
    # 1) Absorption in the A-band
    PO2A1 = gA * densO2

    # 2) Absorption in the B-band, then electronic quenching
    PO2B = gB * densO2

    Lemi = A771
    LO3 = k3B * densO3
    Lqu0 = k0B * densO
    Lqu1 = k1B * densO2
    Lqu2 = k2B * densN2
    LO2B = Lemi + LO3 + Lqu0 + Lqu1 + Lqu2

    PO2A2 = (Lqu0 + Lqu1 + Lqu2) * PO2B / LO2B

    # 3) Collisions with O1D
    PO1D_O2 = J2 * densO2
    PO1D_O3 = J3 * densO3

    Lemi = A1D
    Lqu1 = k1 * densO2
    Lqu2 = k2 * densN2
    LO1D = Lemi + Lqu1 + Lqu2

    PO2A3_O2 = phi * k1 * densO2 * PO1D_O2 / LO1D
    PO2A3_O3 = phi * k1 * densO2 * PO1D_O3 / LO1D

    # 4) Barth process
    densM = densN2 + densO2
    PBar = k1star * densO ** 2 * densM * densO2
    LBar = Co2 * densO2 + Co * densO
    PO2A4 = PBar / LBar

    # -------------------------------------------------------------------------------------------
    # Loss mechanisms of A-band O2
    # -------------------------------------------------------------------------------------------
    Lemi = A1Sig
    Lqu = k1A * densN2 + k2A * densO3 + k3A * densO2 + k4A * densO

    LO2A = Lemi + Lqu

    # -------------------------------------------------------------------------------------------
    # Density profiles of A-band O2
    # -------------------------------------------------------------------------------------------

    densO2Apros = np.zeros((len(T), 5))
    densO2Apros[:, 0] = PO2A1 / LO2A  # due to absorption in A-band
    densO2Apros[:, 1] = PO2A2 / LO2A  # due to absorption in B-band
    densO2Apros[:, 2] = PO2A3_O2 / LO2A  # due to collisions with O1D (from O2 photolysis)
    densO2Apros[:, 3] = PO2A3_O3 / LO2A  # due to collisions with O1D (from O3 photolysis)
    densO2Apros[:, 4] = PO2A4 / LO2A  # due to the Barth process

    return densO2Apros


def runEmissionModel(refPoint, sunPos, upVec, atmosphere, altArray_km):
    """
    Runs the full photochemical emisson model of O2 photon emission in the A-band. Takes in geometry and atmospheric
    parameters and produces the photochemical emission source function for use in the SASKTRAN radiative transfer model.

    Parameters
    ----------
    refPoint: np.ndarray
        Four-element array containing latitude, longitude, altitude, and MJD

    sunPos: np.ndarray
        Three-element array describing sun position in geodetic coordinates

    upVec: np.ndarray
        Three-element array describing direction of local up vector at reference point location

    atmosphere: object
        SASKTRAN object containing temperature and density profiles for O2 nd O3

    altArray_km: np.ndarray
        User-defined altitude array (km) upon which all profiles depend.
        *Must reach at least 200 km altitude for proper solar flux attenuation.

    Returns
    -------
    Jemi: np.ndarray
        Photochemical emission source function (photons/s/cm3/cm-1/ster)
        One spectrum per altitude across the A-band

    wavnumA: np.ndarray
        Wavenumber (cm-1) array used for calculations across the A-band
    """

    # -------------------------------------------------------------------------------------------
    # Set up geometry
    # -------------------------------------------------------------------------------------------
    lat = refPoint[0]
    lon = refPoint[1]
    mjd = refPoint[3]

    dotProd = np.dot(sunPos, upVec)
    sza = np.rad2deg(np.arccos(dotProd))

    # -------------------------------------------------------------------------------------------
    # Set up atmosphere
    # -------------------------------------------------------------------------------------------
    # from SASKTRAN:
    alt_m = altArray_km * 1000
    temp = atmosphere.atmospheric_state.get_parameter('SKCLIMATOLOGY_TEMPERATURE_K', lat, lon, alt_m, mjd)
    densO2 = atmosphere['O2'].climatology.get_parameter('SKCLIMATOLOGY_O2_CM3', lat, lon, alt_m, mjd)
    densO3 = atmosphere['ozone'].climatology.get_parameter('SKCLIMATOLOGY_O3_CM3', lat, lon, alt_m, mjd)

    # from elsewhere:
    MSIS2 = np.genfromtxt(load_data.msis())
    densO = MSIS2[:, 1]
    densN2 = MSIS2[:, 2]

    densO = np.interp(altArray_km, MSIS2[:, 0], densO)
    densN2 = np.interp(altArray_km, MSIS2[:, 0], densN2)

    # -------------------------------------------------------------------------------------------
    # Spectroscopic data and calculations
    # -------------------------------------------------------------------------------------------
    # Get line positions and line strengths for A and B bands from HITRAN
    [nuTransA, nuTransB, Sj_TA, Sj_TB] = getHITRANdata(temp)

    # -------------------------------------------------------------------------------------------
    # Absorption cross-sections
    # -------------------------------------------------------------------------------------------
    logging.info('    Calculating absorption cross-sections')
    t0 = time.time()

    # Get X-section data for Lya, SRC, and Har from data sources
    specDataLya = np.genfromtxt(load_data.xs_lya())
    absXsecLya = np.array([specDataLya[1]])

    specDataSRC = np.genfromtxt(load_data.xs_src())
    absXsecSRC = specDataSRC[:, 1]

    specDataHar = np.genfromtxt(load_data.xs_har())
    absXsecHar = specDataHar[:, 1]

    # Calculate X-sections of O2 in the A and B bands from line strengths
    wavnumA = np.arange(12844, 13165, 0.001)  # must be at least 0.001 nm
    wavnumB = np.arange(14301, 14557, 0.001)  # must be at least 0.001 nm

    absXsecA = np.zeros((wavnumA.shape[0], temp.shape[0]))
    absXsecB = np.zeros((wavnumB.shape[0], temp.shape[0]))

    for i in range(temp.shape[0]):
        absXsecA[:, i], DjA = DopplerBroadening(temp[i], nuTransA, Sj_TA[:, i], wavnumA)
        absXsecB[:, i], DjB = DopplerBroadening(temp[i], nuTransB, Sj_TB[:, i], wavnumB)

    t1 = time.time()
    deltat = t1 - t0
    logging.info('    >> Time elapsed:  ' + str(deltat) + ' sec')

    # -------------------------------------------------------------------------------------------
    # Incident and attenuated solar flux
    # -------------------------------------------------------------------------------------------
    logging.info('    Loading SORCE data')
    t0 = time.time()

    # Get SORCE incident flux data properly parsed -- SORCE data only starts at mjd 52775 (May 15, 2003)
    if mjd < 52775:
        mjd1 = 52775
    else:
        mjd1 = mjd
    date = mjdToYYYYMMDD(mjd1) + 0.5
    [flxDataLya, flxDataSRC, flxDataHar, flxDataA, flxDataB] = getSolFluxData(date)

    t1 = time.time()
    deltat = t1 - t0
    logging.info('    >> Time elapsed:  ' + str(deltat) + ' sec')

    logging.info('    Calculating attenuated solar flux')
    t0 = time.time()

    # Calculate attenuated solar flux spectra
    attSolFlx_Lya, absXsecLya = calcAttSolFlux(sza, altArray_km, densO2, absXsecLya, flxDataLya[0, 1])
    attSolFlx_SRC, absXsecSRC = calcAttSolFlux(sza, altArray_km, densO2, absXsecSRC, flxDataSRC[:, 1])
    attSolFlx_Har, absXsecHar = calcAttSolFlux(sza, altArray_km, densO3, absXsecHar, flxDataHar[:, 1])

    flxDataA = flxDataA[:, 1] * flxDataA[:, 0] ** 2 / 1e7  # convert to wavenumber bins (flx*lambda^2)
    attSolFlx_A, absXsecA = calcAttSolFlux(sza, altArray_km, densO2, absXsecA, flxDataA)

    flxDataB = flxDataB[:, 1] * flxDataB[:, 0] ** 2 / 1e7  # convert to wavenumber bins (flx*lambda^2)
    attSolFlx_B, absXsecB = calcAttSolFlux(sza, altArray_km, densO2, absXsecB, flxDataB)

    t1 = time.time()
    deltat = t1 - t0
    logging.info('    >> Time elapsed:  ' + str(deltat) + ' sec')

    # -------------------------------------------------------------------------------------------
    # Photochemical reaction rate profiles
    # -------------------------------------------------------------------------------------------
    logging.info('    Calculating g and J factors')

    wvlnSRC = flxDataSRC[:, 0]
    wvlnHar = flxDataHar[:, 0]

    gA, gB, JL, JS, JH = calc_gAndJfactors(wvlnSRC, wvlnHar, wavnumA, wavnumB, altArray_km,
                                           attSolFlx_Lya, attSolFlx_SRC, attSolFlx_Har, attSolFlx_A, attSolFlx_B,
                                           absXsecLya, absXsecSRC, absXsecHar, absXsecA, absXsecB)

    J2 = JL + JS
    J3 = JH

    # -------------------------------------------------------------------------------------------
    # Spectral Emission Weighting Function
    # -------------------------------------------------------------------------------------------
    logging.info('    Calculating spectral emission weighting function')
    t0 = time.time()

    # Define some constants
    h = 6.6256e-27  # Planck's constant [cm2*kg]
    c = 2.9979e10  # speed of light [cm/s]
    k = 1.3806e-16  # Boltzmann's constant [cm2*kg/s2*K]

    # Calculate the theoretical spectral VER by Doppler Broadening the transition VER values
    etaj_T = np.zeros(len(nuTransA))
    etaPrSpec = np.zeros((wavnumA.shape[0], temp.shape[0]))
    for i in range(len(temp)):
        for j in range(len(nuTransA)):
            etaj_T[j] = nuTransA[j] ** 2 * Sj_TA[j, i] * np.exp((-h * c * nuTransA[j]) / (k * temp[i]))
        etaPrSpec[:, i], _ = DopplerBroadening(temp[i], nuTransA, etaj_T, wavnumA)

    # Normalize to calculate the weighting function
    dnu = np.mean(np.diff(wavnumA))  # step size of the wavenumber array
    Wemi = np.zeros((wavnumA.shape[0], temp.shape[0]))
    for i in range(temp.shape[0]):
        Wemi[:, i] = etaPrSpec[:, i] / (np.sum(etaPrSpec[:, i]) * dnu)

    t1 = time.time()
    deltat = t1 - t0
    logging.info('    >> Time elapsed:  ' + str(deltat) + ' sec')

    # -------------------------------------------------------------------------------------------
    # Spectral Volume Emission Rate
    # -------------------------------------------------------------------------------------------
    logging.info('    Running photochemical model')
    # Find the density of excited O2A using the photochemical model
    densO2Apros = runPhotochemModel(temp, altArray_km, gA, gB, J2, J3, densO2, densO3, densN2, densO)
    densO2A = densO2Apros[:, 0] + densO2Apros[:, 1] + densO2Apros[:, 2] + densO2Apros[:, 3] + densO2Apros[:, 4]

    logging.info('    Calculating spectral VER')
    # Calculate the integrated VER profile
    Fc = 0.93  # Franck-Condon factor
    A1Sig = 0.085  # Einstein emission coefficient, [s-1]

    etaPro = Fc * A1Sig * densO2A

    # Calculate the actual spectral VER from the weighting function and the photochemical model results
    etaSpec = np.zeros((wavnumA.shape[0], temp.shape[0]))
    for i in range(temp.shape[0]):
        etaSpec[:, i] = Wemi[:, i] * etaPro[i]

    # -------------------------------------------------------------------------------------------
    # Radiative Transfer Source Function
    # -------------------------------------------------------------------------------------------
    # Just needs to be divided by a complete solid angle of 4pi to be used in SASKTRAN's model

    Jemi = etaSpec / (4 * np.pi)

    return Jemi, wavnumA
