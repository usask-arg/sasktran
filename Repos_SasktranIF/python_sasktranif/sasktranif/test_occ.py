"""
Tests that sktran_do modules is operating as expected.

Date: March 28, 2017
Geometry: TEMPO at (20, -100, 35786000), looking at latitudes [-38.51, 26.05], longitudes [-100, -100]
"""

import unittest, math
import numpy as np
import sasktranif.sasktranif as skif

import sys
def printf(format, *args):
    sys.stdout.write(format % args)



#------------------------------------------------------------------------------
#           MakeOpticalStateForOccultationHITRAN
#------------------------------------------------------------------------------

def MakeOpticalStateForOccultationHITRAN(  occengine : skif.ISKEngine, min_wavenum: float, max_wavenum : float):

    rayleigh         = skif.ISKOpticalProperty('Rayleigh' )
    msis90           = skif.ISKClimatology    ('MSIS90')
    co2numberdensity = skif.ISKClimatology    ('USERDEFINED_PROFILE')
    co2numberdensity.SetProperty('DoLogInterpolation', 1)

    co2profile =  np.array( [    0.000, 9.5620350469e+15,  1000.000, 8.5604676285e+15,  2000.000, 7.7062120091e+15,  3000.000, 6.9531991470e+15,  4000.000, 6.2702731320e+15,  5000.000, 5.6375862919e+15,
                    6000.000, 5.0651291274e+15,  7000.000, 4.4975838604e+15,  8000.000, 3.9468136861e+15,  9000.000, 3.4348048814e+15, 10000.000, 2.9871067830e+15, 11000.000, 2.5656416175e+15,
                   12000.000, 2.1874053365e+15, 13000.000, 1.8533816021e+15, 14000.000, 1.6023327829e+15, 15000.000, 1.3568375796e+15, 16000.000, 1.1279532788e+15, 17000.000, 9.7672446573e+14,
                   18000.000, 8.4173283897e+14, 19000.000, 7.1576699275e+14, 20000.000, 6.2070908062e+14, 21000.000, 5.2364297410e+14, 22000.000, 4.3181248841e+14, 23000.000, 3.7860567983e+14,
                   24000.000, 3.2428122678e+14, 25000.000, 2.7110791383e+14, 26000.000, 2.3526785090e+14, 27000.000, 2.0344493146e+14, 28000.000, 1.7304110039e+14, 29000.000, 1.4714113133e+14,
                   30000.000, 1.2544180466e+14, 31000.000, 1.0738346125e+14, 32000.000, 9.2442937053e+13, 33000.000, 8.0342242281e+13, 34000.000, 6.9455591820e+13, 35000.000, 5.9095214441e+13,
                   36000.000, 5.0374561563e+13, 37000.000, 4.3515754800e+13, 38000.000, 3.7794009046e+13, 39000.000, 3.2874895083e+13, 40000.000, 2.8685628465e+13, 41000.000, 2.4978923024e+13,
                   42000.000, 2.1682117851e+13, 43000.000, 1.8864809592e+13, 44000.000, 1.6431826141e+13, 45000.000, 1.4348899126e+13, 46000.000, 1.2595260698e+13, 47000.000, 1.1093125765e+13,
                   48000.000, 9.8376261311e+12, 49000.000, 8.8026864921e+12, 50000.000, 7.8993464447e+12, 51000.000, 7.0038829664e+12, 52000.000, 6.0771348455e+12, 53000.000, 5.2887296427e+12,
                   54000.000, 4.6787494256e+12, 55000.000, 4.1667051367e+12, 56000.000, 3.6751620506e+12, 57000.000, 3.1811011797e+12, 58000.000, 2.7604364326e+12, 59000.000, 2.4249492298e+12,
                   60000.000, 2.1420175118e+12, 61000.000, 1.8772791073e+12, 62000.000, 1.6195294613e+12, 63000.000, 1.3994285676e+12, 64000.000, 1.2229247260e+12, 65000.000, 1.0734951007e+12,
                   66000.000, 9.3270881894e+11, 67000.000, 7.9345730980e+11, 68000.000, 6.7795327304e+11, 69000.000, 5.9174431127e+11, 70000.000, 5.2173619614e+11, 71000.000, 4.5523334147e+11,
                   72000.000, 3.8840635314e+11, 73000.000, 3.3304529951e+11, 74000.000, 2.9045416707e+11, 75000.000, 2.5517516779e+11, 76000.000, 2.2127024526e+11, 77000.000, 1.8582366434e+11,
                   78000.000, 1.5596546276e+11, 79000.000, 1.3362547386e+11, 80000.000, 1.1541990113e+11, 81000.000, 9.8756976417e+10, 82000.000, 8.2629944315e+10, 83000.000, 6.8563739750e+10,
                   84000.000, 5.6814363571e+10, 85000.000, 4.6797966799e+10, 86000.000, 3.8795906044e+10, 87000.000, 3.2908654369e+10, 88000.000, 2.7811184596e+10, 89000.000, 2.2974282383e+10,
                   90000.000, 1.8716304570e+10, 91000.000, 1.5254396937e+10, 92000.000, 1.2548308770e+10, 93000.000, 1.0295593615e+10, 94000.000, 8.3338827301e+09, 95000.000, 6.6488536883e+09,
                   96000.000, 5.2936443303e+09, 97000.000, 4.2242029799e+09, 98000.000, 3.3594428424e+09, 99000.000, 2.6511281727e+09]).reshape( [100,2])

    co2numberdensity.SetProperty('Heights', co2profile[:,0])
    co2numberdensity.SetPropertyUserDefined('SKCLIMATOLOGY_CO2_CM3', co2profile[:,1])

    co2_opticalprops = skif.ISKOpticalProperty('HITRANCHEMICAL_CO2')
    co2_opticalprops.SetProperty( 'SetWavenumberRange', (min_wavenum, max_wavenum) )

    occengine.AddSpecies( 'SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3', msis90, rayleigh)
    occengine.AddSpecies( 'SKCLIMATOLOGY_CO2_CM3', co2numberdensity, co2_opticalprops  )


#------------------------------------------------------------------------------
#           Class QuickOCCTest
#------------------------------------------------------------------------------

class QuickOCCTest(unittest.TestCase):
    """
    Test bench to run sanity checks on SASKTRAN-DISCO. Note that the primary function of these tests if verifying the
    SKIF interface as well as the internal interface between the SASKTRAN framework and the DISCO kernel (eg. making
    sure species are being parameterized properly for the DISCO kernel). This test bench also is used to track
    continuity between releases (eg. track that radiances between releases are the same).

    The primary test bench for DISCO is the ``sasktran_disco_engine_tests`` which is a program that runs lots of edge-
    case tests, but only verifies that the DISCO kernel is operating properly.
    """

    def test_weighting_functions(self):

        min_wavenumber = 934.0                  # minimum wavenumber in cm-1
        max_wavenumber = 936.0                  # maximum wavenumber in cm-1
        res_wavenumber = 0.0005                 # wavenumber resolution in cm-1

        mjd = 54242.26386852                    # MJD("2007-05-22 06:19:58.24") Use ACE-FTS scan ace.sr20314 as an example.
        latitude = 68.91
        longitude = -79.65
        observeralt = 600000.0
        tanalts = np.array( [ 95.542934418655, 93.030998230911, 90.518486023880, 87.999366761185, 85.485855103470, 82.971916199661, 80.457603455521, 77.942962647415, 75.421955109573, 72.906806946732,
                     70.391479493118, 67.869934082962, 65.354396820999, 62.838825226761, 60.317182541824, 57.801700592972, 55.286336899734, 52.765050888992, 50.250070572830, 47.735359192825,
                     45.220966340042, 42.682148825007, 40.254586233282, 37.901745439957, 35.689252976524, 33.619203107470, 31.633878541417, 29.706157206720, 27.941217916525, 26.315136637345,
                     24.759740931714, 23.273057386890, 21.701357220703, 20.435203333687, 19.296175927224, 18.238125008002, 17.137857798933, 15.665431416870, 14.623809766528, 13.581115284387,
                     12.793781944543, 11.901170623281, 10.978181776555, 10.1851695349872, 9.4383271471788, 8.7424541473265, 8.0540969039894, 7.5483134223615, 7.0824804787830, 6.7903857771487,
                     6.3015475934096] )*1000.0

        shellheights = np.arange( 0.0, 100000.0, 1000.0 )
        wavenumber = np.arange( min_wavenumber, max_wavenumber, res_wavenumber)
        occengine  = skif.ISKEngine('OCC')
        occengine.SetWavelengths( 1.0E7/wavenumber)
        occengine.SetProperty( 'SetReferencePoint', ( latitude, longitude, mjd ) )
        occengine.SetProperty( 'SetRayTracingShells', shellheights  )

        for h in tanalts: occengine.SetProperty( 'AddLineOfSightFromTangentAlt', [h, observeralt] )             # Add a line of sight based upon the tangent altitude



        MakeOpticalStateForOccultationHITRAN( occengine, min_wavenumber, max_wavenumber )
        occengine.InitializeModel()
        [ok,extinction] = occengine.CalculateRadiance()



if __name__ == "__main__":
    tests = QuickOCCTest()
    unittest.main()

