from sasktran.legendre.wigner import WignerD
import numpy as np
from scipy.integrate import simpson


def compute_greek_coefficients_legendre(lm_p11: np.array, lm_p12: np.array, lm_p22: np.array, lm_p33: np.array, lm_p34: np.array, lm_p44: np.array,
                                        theta_grid=None):
    """
    Calculates the greek coefficients a1, a2, a3, a4, b1, b2 given the Legendre expansion of the phase function
    elements.

    Parameters
    ----------
    lm_p11 : np.array
        Legendre moments for P11
    lm_p12 : np.array
        Legendre moments for P12
    lm_p22 : np.array
        Legendre moments for P22
    lm_p33 : np.array
        Legendre moments for P33
    lm_p34 : np.array
        Legendre moments for P34
    lm_p44 : np.array
        Legendre moments for P44
    theta_grid : np.array, Optional
        Grid to perform the final integral over in degrees.  Should fully span 0 to 180

    Returns
    -------
    lm_a1 : np.array
        Greek coefficients for a1
    lm_a2 : np.array
        Greek coefficients for a2
    lm_a3 : np.array
        Greek coefficients for a3
    lm_a4 : np.array
        Greek coefficients for a4
    lm_b1 : np.array
        Greek coefficients for b1
    lm_b2 : np.array
        Greek coefficients for b2
    """
    if theta_grid is None:
        theta_grid = np.arange(0, 180.0001, 0.01)

    num_coeff = len(lm_p11)

    wigner0 = WignerD(0, 0)

    p11 = np.zeros(len(theta_grid))
    p12 = np.zeros_like(p11)
    p22 = np.zeros_like(p11)
    p33 = np.zeros_like(p11)
    p34 = np.zeros_like(p11)
    p44 = np.zeros_like(p11)

    for i in range(num_coeff):
        lpoly = wigner0.d(np.deg2rad(theta_grid), i)

        p11 += lpoly * lm_p11[i]
        p12 += lpoly * lm_p12[i]
        p22 += lpoly * lm_p22[i]
        p33 += lpoly * lm_p33[i]
        p34 += lpoly * lm_p34[i]
        p44 += lpoly * lm_p44[i]

    return compute_greek_coefficients(p11, p12, p22, p33, p34, p44, theta_grid, num_coeff)


def compute_greek_coefficients(p11: np.array, p12: np.array, p22: np.array, p33: np.array, p34: np.array, p44: np.array, angle_grid: np.array, num_coeff: int):
    """
    Calculates the greek coefficients a1, a2, a3, a4, b1, b2 given the Legendre expansion of the phase function
    elements.

    Parameters
    ----------
    p11 : np.array
        Phase function P11
    p12 : np.array
        Phase function P12
    p22 : np.array
        Phase function P22
    p33 : np.array
        Phase function P33
    p34 : np.array
        Phase function P34
    p44 : np.array
        Phase function P44
    angle_grid : np.array
        Angular grid the phase functions are specified on .  Should fully span 0 to 180
    num_coeff: int
        Maximum number of coefficients to return in the expansion.

    Returns
    -------
    lm_a1 : np.array
        Greek coefficients for a1
    lm_a2 : np.array
        Greek coefficients for a2
    lm_a3 : np.array
        Greek coefficients for a3
    lm_a4 : np.array
        Greek coefficients for a4
    lm_b1 : np.array
        Greek coefficients for b1
    lm_b2 : np.array
        Greek coefficients for b2
    """
    cos_theta = np.cos(np.deg2rad(angle_grid))
    angle_rad = np.deg2rad(angle_grid)

    wigner00 = WignerD(0, 0)
    wigner22 = WignerD(2, 2)
    wigner2m2 = WignerD(2, -2)
    wigner02 = WignerD(0, 2)

    lm_a1 = np.zeros(num_coeff)
    lm_a2 = np.zeros(num_coeff)
    lm_a3 = np.zeros(num_coeff)
    lm_a4 = np.zeros(num_coeff)
    lm_b1 = np.zeros(num_coeff)
    lm_b2 = np.zeros(num_coeff)

    for i in range(0, num_coeff):
        lpoly_00 = wigner00.d(angle_rad, i)
        lpoly_22 = wigner22.d(angle_rad, i)
        lpoly_2m2 = wigner2m2.d(angle_rad, i)
        lpoly_02 = wigner02.d(angle_rad, i)

        lm_a1[i] = simpson(p11 * lpoly_00, cos_theta) / simpson(lpoly_00 * lpoly_00, cos_theta)
        lm_a4[i] = simpson(p44 * lpoly_00, cos_theta) / simpson(lpoly_00 * lpoly_00, cos_theta)

        if i >= 2:
            lm_b1[i] = simpson(p12 * lpoly_02, cos_theta) / simpson(lpoly_02 * lpoly_02, cos_theta)
            lm_b2[i] = -simpson(p34 * lpoly_02, cos_theta) / simpson(lpoly_02 * lpoly_02, cos_theta)

            temp1 = simpson((p22 + p33) * lpoly_22, cos_theta) / simpson(lpoly_22 * lpoly_22, cos_theta)
            temp2 = simpson((p22 - p33) * lpoly_2m2, cos_theta) / simpson(lpoly_2m2 * lpoly_2m2, cos_theta)

            lm_a2[i] = (temp1 + temp2) / 2
            lm_a3[i] = (temp1 - temp2) / 2

    return lm_a1, lm_a2, lm_a3, lm_a4, lm_b1, lm_b2
