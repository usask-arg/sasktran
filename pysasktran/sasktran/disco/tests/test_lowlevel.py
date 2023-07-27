import sasktran.disco.lowlevel as lowlevel
import numpy as np


def test_coulsen_albedo_08():
    """
    Performs a single layer test from the Coulsen tables.  These tests are done in multiple places, this is
    mostly just a sanity check on the lowlevel interface inputs.
    """
    nstr = 40
    nlyr = 1
    nwavel = 1
    nderiv = 0
    nstokes = 3
    nlos = 6

    atmosphere = lowlevel.Atmosphere(nstr, nlyr, nwavel)
    config = lowlevel.Config(nstr, nwavel, nlyr, nstokes, 0, use_pseudo_spherical=False)
    weightingfunctions = lowlevel.WeightingFunctions(nstr, nlos, nstokes, nwavel, nderiv)
    viewing_geometry = lowlevel.ViewingGeometry(nlos)

    atmosphere.od[0, :] = 0.5

    atmosphere.ssa[0, :] = 1.0

    atmosphere.a1[0, 0, :] = 1
    atmosphere.a1[2, 0, :] = 0.5

    atmosphere.a2[2, 0, :] = 3

    atmosphere.a4[1, 0, :] = 3/2

    atmosphere.b1[2, 0, :] = np.sqrt(6.0) * 0.5

    # Set this to something really small to essentially disable the pesudo-spherical correction
    atmosphere.layer_boundary_altitudes[0] = 1000000

    atmosphere.albedo[0] = 0.8

    viewing_geometry.cos_sza = 0.2
    viewing_geometry.cos_vza[0] = 0.02
    viewing_geometry.cos_vza[1] = 0.4
    viewing_geometry.cos_vza[2] = 1.00
    viewing_geometry.cos_vza[3] = 0.02
    viewing_geometry.cos_vza[4] = 0.4
    viewing_geometry.cos_vza[5] = 1.00

    viewing_geometry.saa[3] = -np.pi/3
    viewing_geometry.saa[4] = -np.pi/3
    viewing_geometry.saa[5] = -np.pi/3

    output = lowlevel.calculate(atmosphere, config, weightingfunctions, viewing_geometry).xarray().isel(wavelength=0)['radiance'].values * np.pi

    true = np.array([0.47382125, -0.01553672, 0,
                     0.23059806, 0.01144320, 0,
                     0.13280858, 0.03755859, 0,
                     0.33343531, -0.15766132, 0.07365528,
                     0.18923236, -0.06041229, 0.05293867,
                     0.13280858, -0.01877930, 0.03252669]).reshape((6, 3)).T

    np.testing.assert_allclose(output, true, rtol=1e-3, atol=1e-4)


def test_multiple_layer():
    """
    Repeats the coulsen albedo test using multiple layers, ensuring that the answer is basically the same when
    multiple layers are used.
    """
    nstr = 40
    nlyr = 1
    nwavel = 1
    nderiv = 0
    nstokes = 3
    nlos = 6

    atmosphere = lowlevel.Atmosphere(nstr, nlyr, nwavel)
    config = lowlevel.Config(nstr, nwavel, nlyr, nstokes, 0, use_pseudo_spherical=False)
    weightingfunctions = lowlevel.WeightingFunctions(nstr, nlos, nstokes, nwavel, nderiv)
    viewing_geometry = lowlevel.ViewingGeometry(nlos)

    atmosphere.od[0, :] = 0.5

    atmosphere.ssa[0, :] = 1.0

    atmosphere.a1[0, 0, :] = 1
    atmosphere.a1[2, 0, :] = 0.5

    atmosphere.a2[2, 0, :] = 3

    atmosphere.a4[1, 0, :] = 3/2

    atmosphere.b1[2, 0, :] = np.sqrt(6.0) * 0.5

    # Not really used?
    atmosphere.layer_boundary_altitudes[0] = 0.00001

    atmosphere.albedo[0] = 0

    viewing_geometry.cos_sza = 0.2
    viewing_geometry.cos_vza[0] = 0.02
    viewing_geometry.cos_vza[1] = 0.4
    viewing_geometry.cos_vza[2] = 1.00
    viewing_geometry.cos_vza[3] = 0.02
    viewing_geometry.cos_vza[4] = 0.4
    viewing_geometry.cos_vza[5] = 1.00

    viewing_geometry.saa[3] = np.pi/3
    viewing_geometry.saa[4] = np.pi/3
    viewing_geometry.saa[5] = np.pi/3

    output = lowlevel.calculate(atmosphere, config, weightingfunctions, viewing_geometry).xarray()['radiance']

    nstr = 40
    nlyr = 10
    nwavel = 1
    nderiv = 0
    nstokes = 3
    nlos = 6

    atmosphere = lowlevel.Atmosphere(nstr, nlyr, nwavel)
    config = lowlevel.Config(nstr, nwavel, nlyr, nstokes, 0, use_pseudo_spherical=False)
    weightingfunctions = lowlevel.WeightingFunctions(nstr, nlos, nstokes, nwavel, nderiv)
    viewing_geometry = lowlevel.ViewingGeometry(nlos)

    atmosphere.od[:, :] = 0.5 / nlyr

    atmosphere.ssa[:, :] = 1.0

    atmosphere.a1[0, :, :] = 1
    atmosphere.a1[2, :, :] = 0.5

    atmosphere.a2[2, :, :] = 3

    atmosphere.a4[1, :, :] = 3/2

    atmosphere.b1[2, :, :] = np.sqrt(6.0) * 0.5

    atmosphere.layer_boundary_altitudes[:] = np.linspace(0.0001, 0.00001, nlyr)

    atmosphere.albedo[0] = 0

    viewing_geometry.cos_sza = 0.2
    viewing_geometry.cos_vza[0] = 0.02
    viewing_geometry.cos_vza[1] = 0.4
    viewing_geometry.cos_vza[2] = 1.00
    viewing_geometry.cos_vza[3] = 0.02
    viewing_geometry.cos_vza[4] = 0.4
    viewing_geometry.cos_vza[5] = 1.00

    viewing_geometry.saa[3] = np.pi/3
    viewing_geometry.saa[4] = np.pi/3
    viewing_geometry.saa[5] = np.pi/3

    output2 = lowlevel.calculate(atmosphere, config, weightingfunctions, viewing_geometry).xarray()['radiance']

    np.testing.assert_allclose(output.values, output2.values, atol=1e-7, rtol=1e-6)


def test_single_layer_wf_scalar():
    nstr = 16
    nlyr = 1
    nwavel = 1
    nderiv = 4 + nstr
    nstokes = 1
    nlos = 6

    atmosphere = lowlevel.Atmosphere(nstr, nlyr, nwavel)
    config = lowlevel.Config(nstr, nwavel, nlyr, nstokes, 0)
    weightingfunctions = lowlevel.WeightingFunctions(nstr, nlos, nstokes, nwavel, nderiv)
    viewing_geometry = lowlevel.ViewingGeometry(nlos)

    atmosphere.od[0, :] = 0.5

    atmosphere.ssa[0, :] = 0.8

    atmosphere.a1[0, 0, :] = 1
    atmosphere.a1[2, 0, :] = 0.5

    atmosphere.f[0, :] = 0.01

    atmosphere.albedo[0] = 0.2

    atmosphere.layer_boundary_altitudes[0] = 100000

    viewing_geometry.cos_sza = 0.2
    viewing_geometry.cos_vza[0] = 0.02
    viewing_geometry.cos_vza[1] = 0.4
    viewing_geometry.cos_vza[2] = 1.00
    viewing_geometry.cos_vza[3] = 0.02
    viewing_geometry.cos_vza[4] = 0.4
    viewing_geometry.cos_vza[5] = 1.00

    viewing_geometry.saa[3] = np.pi/3
    viewing_geometry.saa[4] = np.pi/3
    viewing_geometry.saa[5] = np.pi/3

    weightingfunctions.d_layerindex[:] = np.zeros(nderiv)

    weightingfunctions.d_od[0, 0] = 1
    weightingfunctions.d_ssa[1, 0] = 1
    weightingfunctions.d_f[2, 0] = 1

    weightingfunctions.d_albedo[3, 0] = 1

    for idx in range(4, nderiv):
        weightingfunctions.d_a1[idx, idx-4, 0] = 1

    output = lowlevel.calculate(atmosphere, config, weightingfunctions, viewing_geometry)

    wf = output._d_radiance.reshape((nlos*nstokes, nderiv))

    deriv_delta_fraction = 0.0001
    for idx in range(nderiv):
        deriv_delta = 0

        atmosphere.ssa[0, 0] += weightingfunctions.d_ssa[idx, 0] * deriv_delta_fraction
        deriv_delta += weightingfunctions.d_ssa[idx, 0] * deriv_delta_fraction

        atmosphere.od[0, 0] += weightingfunctions.d_od[idx, 0] * deriv_delta_fraction
        deriv_delta += weightingfunctions.d_od[idx, 0] * deriv_delta_fraction

        atmosphere.f[0, 0] += weightingfunctions.d_f[idx, 0] * deriv_delta_fraction
        deriv_delta += weightingfunctions.d_f[idx, 0] * deriv_delta_fraction

        atmosphere.albedo[0] += weightingfunctions.d_albedo[idx, 0] * deriv_delta_fraction
        deriv_delta += weightingfunctions.d_albedo[idx, 0] * deriv_delta_fraction

        for idy in range(nstr):
            atmosphere.a1[idy, 0, 0] += weightingfunctions.d_a1[idx, idy, 0] * deriv_delta_fraction
            deriv_delta += weightingfunctions.d_a1[idx, idy, 0] * deriv_delta_fraction

        output_p = lowlevel.calculate(atmosphere, config, weightingfunctions, viewing_geometry)

        numerical_wf = (output_p._radiance - output._radiance) / deriv_delta

        np.testing.assert_array_almost_equal(wf[:, idx], numerical_wf, 4)

        atmosphere.ssa[0, 0] -= weightingfunctions.d_ssa[idx, 0] * deriv_delta_fraction
        atmosphere.od[0, 0] -= weightingfunctions.d_od[idx, 0] * deriv_delta_fraction
        atmosphere.f[0, 0] -= weightingfunctions.d_f[idx, 0] * deriv_delta_fraction
        atmosphere.albedo[0] -= weightingfunctions.d_albedo[idx, 0] * deriv_delta_fraction

        for idy in range(nstr):
            atmosphere.a1[idy, 0, 0] -= weightingfunctions.d_a1[idx, idy, 0] * deriv_delta_fraction


def test_ss_phase_wf():
    nstr = 16
    nlyr = 1
    nwavel = 1
    nderiv = 1
    nstokes = 1
    nlos = 6

    atmosphere = lowlevel.Atmosphere(nstr, nlyr, nwavel, nlos, nstokes)
    config = lowlevel.Config(nstr, nwavel, nlyr, nstokes, 0, exact_ss=True)
    weightingfunctions = lowlevel.WeightingFunctions(nstr, nlos, nstokes, nwavel, nderiv)
    viewing_geometry = lowlevel.ViewingGeometry(nlos)

    atmosphere.od[0, :] = 0.5

    atmosphere.ssa[0, :] = 0.8

    atmosphere.a1[0, 0, :] = 1
    atmosphere.a1[2, 0, :] = 0.5

    atmosphere.f[0, :] = 0.01

    atmosphere.albedo[0] = 0.2

    atmosphere.ss_phase[0, :, 0, 0] = 0.1

    atmosphere.layer_boundary_altitudes[0] = 100000

    viewing_geometry.cos_sza = 0.2
    viewing_geometry.cos_vza[0] = 0.02
    viewing_geometry.cos_vza[1] = 0.4
    viewing_geometry.cos_vza[2] = 1.00
    viewing_geometry.cos_vza[3] = 0.02
    viewing_geometry.cos_vza[4] = 0.4
    viewing_geometry.cos_vza[5] = 1.00

    viewing_geometry.saa[3] = np.pi/3
    viewing_geometry.saa[4] = np.pi/3
    viewing_geometry.saa[5] = np.pi/3

    weightingfunctions.d_layerindex[:] = np.zeros(nderiv)

    weightingfunctions.d_ss_phase[0, 0, :, 0] = 1

    output = lowlevel.calculate(atmosphere, config, weightingfunctions, viewing_geometry)

    atmosphere.ss_phase[:] += 0.01

    output_p = lowlevel.calculate(atmosphere, config, weightingfunctions, viewing_geometry)

    wf = (output_p._radiance - output._radiance) / 0.01

    np.testing.assert_allclose(output._d_radiance, wf)
