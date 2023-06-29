
#include <skopticalproperties21.h>
/*
O2-O2 CIA Primary Reference:
----------------------------
Karman et al. 2019: Update of the HITRAN collission-induced absorption section. 
Icarus, 328 (2019), 160-175,
https://doi.org/10.1016/j.icarus.2019.02.034

We have only implemented the "Main" CIA folder and have not implemented the "Alternate" folders discussed in Karman et al. 2019

CIA data files downloaded from: https://hitran.org/cia/
and then edited to be compatible with C-data structures (inserting commas, curly braces and variable names only)
*/

#include "O2-O2_1150.000_1950.000_4002_193.4_8.767E-45_0.500_P-up-to-4atm_2.h"				// Region 1. Temperature dependent
#include "O2-O2_1150.000_1950.000_4002_206.8_8.282E-45_0.500_P-up-to-4atm_2.h"
#include "O2-O2_1150.000_1950.000_4002_218.6_7.854E-45_0.500_P-up-to-4atm_2.h"
#include "O2-O2_1150.000_1950.000_4002_220.8_7.882E-45_0.500_P-up-to-4atm_2.h"
#include "O2-O2_1150.000_1950.000_4002_229.4_7.580E-45_0.500_P-up-to-4atm_2.h"
#include "O2-O2_1150.000_1950.000_4002_229.6_7.704E-45_0.500_P-up-to-4atm_2.h"
#include "O2-O2_1150.000_1950.000_4002_240.0_7.332E-45_0.500_P-up-to-4atm_2.h"
#include "O2-O2_1150.000_1950.000_4002_249.4_7.124E-45_0.500_P-up-to-4atm_2.h"
#include "O2-O2_1150.000_1950.000_4002_270.7_6.849E-45_0.500_P-up-to-4atm_2.h"
#include "O2-O2_1150.000_1950.000_4002_297.5_6.580E-45_0.500_P-up-to-4atm_2.h"
#include "O2-O2_1150.000_1950.000_4002_297.8_6.468E-45_0.500_P-up-to-4atm_2.h"
#include "O2-O2_1150.000_1950.000_4002_320.5_6.370E-45_0.500_P-up-to-4atm_2.h"
#include "O2-O2_1150.000_1950.000_4002_330.8_6.461E-45_0.500_P-up-to-4atm_2.h"
#include "O2-O2_1150.000_1950.000_4002_344.9_6.422E-45_0.500_P-up-to-4atm_2.h"
#include "O2-O2_1150.000_1950.000_4002_353.4_6.470E-45_0.500_P-up-to-4atm_2.h"

#include "O2-O2_7450.132_8491.628_4254_296.0_2.131e-45_0.500_Updated-from-Mate-1999_10.h"	// Region 2. Temperature independent
#include "O2-O2_9091.0000_9596.0000_506_293.00_1.364e-45_-.999_a1-X0_27.h"					// Region 3. Temperature independent
#include "O2-O2_10512.000_11228.000_717_293.00_7.092e-47_-.999_a2-X0_29.h"					// Region 4. Temperature independent
#include "O2-O2_12600.120_13839.642_2572_296.0_2.636E-46_0.500_A-Band-P-20-200atm_4.h"		// Region 5. Temperature Independent
#include "O2-O2_14206.000_14898.000_693_293.00_2.494e-47_-.999_b1-X0-B-band_30.h"		// Region 6. Temperature Independent

#include "O2-O2_15290.000_16664.000_1375_203.00_8.645e-46_-.999_Double-transitions_32.h"	// Region 7. Temperature dependent
#include "O2-O2_15290.000_16664.000_1375_233.00_7.692e-46_-.999_Double-transitions_32.h"
#include "O2-O2_15290.000_16664.000_1375_287.00_7.111e-46_-.999_Double-transitions_32.h"
#include "O2-O2_15292.000_16664.000_1373_253.00_7.682e-46_-.999_Double-transitions_32.h"

#include "O2-O2_16645.000_29784.000_13140_293.00_7.821e-44_-.999_Double-transitions_32.h"	// Region 8. Temperature dependent
#include "O2-O2_16658.000_29748.000_13091_203.00_5.196e-42_-.999_Double-transitions_32.h"
#include "O2-O2_16668.000_29802.000_13135_273.00_6.522e-44_-.999_Double-transitions_32.h"
#include "O2-O2_16780.000_29757.000_12978_233.00_2.356e-43_-.999_Double-transitions_32.h"
#include "O2-O2_16791.000_29837.000_13047_253.00_2.946e-42_-.999_Double-transitions_32.h"

/*---------------------------------------------------------------------------
 * skOpticalProperties_O4_HitranEntry_NotDependent::skOpticalProperties_O4_HitranEntry_NotDependent 2020-01-02 */
/** **/
/*---------------------------------------------------------------------------*/
skOpticalProperties_O4_HitranEntry_NotDependent::skOpticalProperties_O4_HitranEntry_NotDependent( int regionid)
{
	bool	ok = false;

	switch (regionid)
	{
		case 2:	ok = ConfigureAsRegion2(); break;
		case 3: ok = ConfigureAsRegion3(); break;
		case 4: ok = ConfigureAsRegion4(); break;
		case 5: ok = ConfigureAsRegion5(); break;
		case 6: ok = ConfigureAsRegion6(); break;
	};
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_O4_HitranEntry_TempDependent::Constructor,Error creating cross-sections for spectral region %d", (int)regionid);
	}
}



/*---------------------------------------------------------------------------
 * skOpticalProperties_O4_HitranEntry_NotDependent::ConfigureAsRegion2 2020-01-02 */
/**
 *	Region 2 is from Mate et al. 1999
 *	Data points are 0.2449 cm-1 apart
 *  Resolution is about 0.05 cm-1
 *	Wavenumber Range = 7450 - 8491 cm-1
 *  Wavelength Range = 1178 - 1342 nm
  **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_O4_HitranEntry_NotDependent::ConfigureAsRegion2()
{
	bool ok;
	ok = AddAscendingWavenumberEntry( 296.0, &o2o2_7450_8491_296p0[0], 2*sizeof(o2o2_7450_8491_296p0[0]), &o2o2_7450_8491_296p0[1], 2*sizeof(o2o2_7450_8491_296p0[0]),  N_ELEMENTS(o2o2_7450_8491_296p0)/2 );
	SetInstrumentPSF_FWHMWavenumber		( 0.05 );
	SetInstrumentPointSpacingWavenumber	( 0.2449 );
	return ok;
}



/*---------------------------------------------------------------------------
 * skOpticalProperties_O4_HitranEntry_NotDependent::ConfigureAsRegion3 2020-01-02 */
/** 
 *	Region 3 is from Karman et al. 2018
 *	Data points are 1.0 cm-1 apart
 *  Resolution is about 0.01 cm-1 (Actual value could not be found in article)
 *	Wavenumber Range = 9091 - 9596 cm-1
 *  Wavelength Range = 1042 - 1100 nm
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_O4_HitranEntry_NotDependent::ConfigureAsRegion3()
{
	bool ok;
	ok = AddAscendingWavenumberEntry( 293.0, &o2o2_9091_956_293p0[0], 2*sizeof(o2o2_9091_956_293p0[0]), &o2o2_9091_956_293p0[1], 2*sizeof(o2o2_9091_956_293p0[0]),  N_ELEMENTS(o2o2_9091_956_293p0)/2 );
	SetInstrumentPSF_FWHMWavenumber		( 0.01);
	SetInstrumentPointSpacingWavenumber	( 1.0 );
	return ok;
}


/*---------------------------------------------------------------------------
 * skOpticalProperties_O4_HitranEntry_NotDependent::ConfigureAsRegion4 2020-01-02 */
/** 
 *	Region 4 is from Spiering and van der Zande, 2012 
 *	Data points are 1.0 cm-1 apart
 *  Resolution is better than 0.0143 cm-1
 *	Wavenumber Range = 10512 - 11228 cm-1
 *  Wavelength Range =   891 - 951   nm
**/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_O4_HitranEntry_NotDependent::ConfigureAsRegion4()
{
	bool ok;
	ok = AddAscendingWavenumberEntry( 293.0, &o2o2_10512_11228_293p0[0], 2*sizeof(o2o2_10512_11228_293p0[0]), &o2o2_10512_11228_293p0[1], 2*sizeof(o2o2_10512_11228_293p0[0]),  N_ELEMENTS(o2o2_10512_11228_293p0)/2 );
	SetInstrumentPSF_FWHMWavenumber		( 0.0143 );
	SetInstrumentPointSpacingWavenumber	( 1.0);

	return ok;
}

/*---------------------------------------------------------------------------
 * skOpticalProperties_O4_HitranEntry_NotDependent::ConfigureAsRegion5 2020-01-02 */
/** 
 *	Region 5 is from Tran et al, 2006 
 *	Data points are 1.0 cm-1 apart
 *  Resolution is better than 0.5 cm-1
 *	Wavenumber Range = 12600 - 13839 cm-1
 *  Wavelength Range =   722 - 793   nm
**/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_O4_HitranEntry_NotDependent::ConfigureAsRegion5()
{
	bool ok;
	ok = AddAscendingWavenumberEntry( 296.0, &o2o2_12600_13839_296p0[0], 2*sizeof(o2o2_12600_13839_296p0[0]), &o2o2_12600_13839_296p0[1], 2*sizeof(o2o2_12600_13839_296p0[0]),  N_ELEMENTS(o2o2_12600_13839_296p0)/2 );
	SetInstrumentPSF_FWHMWavenumber		( 0.5 );
	SetInstrumentPointSpacingWavenumber	( 1.0);
	return ok;
}


/*---------------------------------------------------------------------------
 * skOpticalProperties_O4_HitranEntry_NotDependent::ConfigureAsRegion6 2020-01-02 */
/**
 *	Region 5 is from Spiering et al, 2011 
 *	Data points are 1.0 cm-1 apart
 *  Resolution is around 0.07 cm-1
 *	Wavenumber Range = 14206 - 14898 cm-1
 *  Wavelength Range =   671 - 704   nm

**/
/*---------------------------------------------------------------------------*/
bool skOpticalProperties_O4_HitranEntry_NotDependent::ConfigureAsRegion6()
{
	bool ok;
	ok = AddAscendingWavenumberEntry( 293.0, &o2o2_14206_14898_293[0], 2*sizeof(o2o2_14206_14898_293[0]), &o2o2_14206_14898_293[1], 2*sizeof(o2o2_14206_14898_293[0]),  N_ELEMENTS(o2o2_14206_14898_293)/2 );
	SetInstrumentPSF_FWHMWavenumber		( 0.07 );
	SetInstrumentPointSpacingWavenumber	( 1.0);
	return ok;
}



/*---------------------------------------------------------------------------
 * skOpticalProperties_O4_HitranEntry_TempDependent::skOpticalProperties_O4_HitranEntry_TempDependent 2020-01-02 */
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_O4_HitranEntry_TempDependent::skOpticalProperties_O4_HitranEntry_TempDependent( int regionid)
{
	bool	ok = false;

	SetVacuumWavelengths(true);					// The wavelength tables are for air not Vacuum 
	SetQuietWavelengthTruncation(true);			// quietly set cross-sections outside the wavelength range to zero.

	switch (regionid)
	{
		case 1:	ok = ConfigureAsRegion1(); break;
		case 7: ok = ConfigureAsRegion7(); break;
		case 8: ok = ConfigureAsRegion8(); break;
	};
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_O4_HitranEntry_TempDependent::Constructor,Error creating cross-sections for spectral region %d", (int)regionid);
	}
}


/*---------------------------------------------------------------------------
 * skOpticalProperties_O4_HitranEntry_TempDependent::ConfigureAsRegion2 2020-01-02 */
/** 
 *	Region 1 is from Baranov 2004. 
 *	Data points are 1.0 cm-1 apart
 *  resolution is about 0.43 cm-1
 *	Wavenumber range = 1150 - 1950  cm-1, 
 *	Wavelength range =  5128 - 8695 nm .
**/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_O4_HitranEntry_TempDependent::ConfigureAsRegion1()
{
	bool ok = true;

	ClearEntries();

	ok = ok && AddAscendingWavenumberEntry( 193.4, &o2o2_1150_1950_193p4[0], 2*sizeof(o2o2_1150_1950_193p4[0]), &o2o2_1150_1950_193p4[1], 2*sizeof(o2o2_1150_1950_193p4[0]),  N_ELEMENTS(o2o2_1150_1950_193p4)/2 );
	ok = ok && AddAscendingWavenumberEntry( 206.8, &o2o2_1150_1950_206p8[0], 2*sizeof(o2o2_1150_1950_206p8[0]), &o2o2_1150_1950_206p8[1], 2*sizeof(o2o2_1150_1950_206p8[0]),  N_ELEMENTS(o2o2_1150_1950_206p8)/2 );
	ok = ok && AddAscendingWavenumberEntry( 218.6, &o2o2_1150_1950_218p6[0], 2*sizeof(o2o2_1150_1950_218p6[0]), &o2o2_1150_1950_218p6[1], 2*sizeof(o2o2_1150_1950_218p6[0]),  N_ELEMENTS(o2o2_1150_1950_218p6)/2 );
	ok = ok && AddAscendingWavenumberEntry( 220.8, &o2o2_1150_1950_220p8[0], 2*sizeof(o2o2_1150_1950_220p8[0]), &o2o2_1150_1950_220p8[1], 2*sizeof(o2o2_1150_1950_220p8[0]),  N_ELEMENTS(o2o2_1150_1950_220p8)/2 );
	ok = ok && AddAscendingWavenumberEntry( 229.4, &o2o2_1150_1950_229p4[0], 2*sizeof(o2o2_1150_1950_229p4[0]), &o2o2_1150_1950_229p4[1], 2*sizeof(o2o2_1150_1950_229p4[0]),  N_ELEMENTS(o2o2_1150_1950_229p4)/2 );
	ok = ok && AddAscendingWavenumberEntry( 229.6, &o2o2_1150_1950_229p6[0], 2*sizeof(o2o2_1150_1950_229p6[0]), &o2o2_1150_1950_229p6[1], 2*sizeof(o2o2_1150_1950_229p6[0]),  N_ELEMENTS(o2o2_1150_1950_229p6)/2 );
	ok = ok && AddAscendingWavenumberEntry( 240.0, &o2o2_1150_1950_240p0[0], 2*sizeof(o2o2_1150_1950_240p0[0]), &o2o2_1150_1950_240p0[1], 2*sizeof(o2o2_1150_1950_240p0[0]),  N_ELEMENTS(o2o2_1150_1950_240p0)/2 );
	ok = ok && AddAscendingWavenumberEntry( 249.4, &o2o2_1150_1950_249p4[0], 2*sizeof(o2o2_1150_1950_249p4[0]), &o2o2_1150_1950_249p4[1], 2*sizeof(o2o2_1150_1950_249p4[0]),  N_ELEMENTS(o2o2_1150_1950_249p4)/2 );
	ok = ok && AddAscendingWavenumberEntry( 240.7, &o2o2_1150_1950_270p7[0], 2*sizeof(o2o2_1150_1950_270p7[0]), &o2o2_1150_1950_270p7[1], 2*sizeof(o2o2_1150_1950_270p7[0]),  N_ELEMENTS(o2o2_1150_1950_270p7)/2 );
	ok = ok && AddAscendingWavenumberEntry( 297.5, &o2o2_1150_1950_297p5[0], 2*sizeof(o2o2_1150_1950_297p5[0]), &o2o2_1150_1950_297p5[1], 2*sizeof(o2o2_1150_1950_297p5[0]),  N_ELEMENTS(o2o2_1150_1950_297p5)/2 );
	ok = ok && AddAscendingWavenumberEntry( 297.8, &o2o2_1150_1950_297p8[0], 2*sizeof(o2o2_1150_1950_297p8[0]), &o2o2_1150_1950_297p8[1], 2*sizeof(o2o2_1150_1950_297p8[0]),  N_ELEMENTS(o2o2_1150_1950_297p8)/2 );
	ok = ok && AddAscendingWavenumberEntry( 320.5, &o2o2_1150_1950_320p5[0], 2*sizeof(o2o2_1150_1950_320p5[0]), &o2o2_1150_1950_320p5[1], 2*sizeof(o2o2_1150_1950_320p5[0]),  N_ELEMENTS(o2o2_1150_1950_320p5)/2 );
	ok = ok && AddAscendingWavenumberEntry( 330.8, &o2o2_1150_1950_330p8[0], 2*sizeof(o2o2_1150_1950_330p8[0]), &o2o2_1150_1950_330p8[1], 2*sizeof(o2o2_1150_1950_330p8[0]),  N_ELEMENTS(o2o2_1150_1950_330p8)/2 );
	ok = ok && AddAscendingWavenumberEntry( 344.9, &o2o2_1150_1950_344p9[0], 2*sizeof(o2o2_1150_1950_344p9[0]), &o2o2_1150_1950_344p9[1], 2*sizeof(o2o2_1150_1950_344p9[0]),  N_ELEMENTS(o2o2_1150_1950_344p9)/2 );
	ok = ok && AddAscendingWavenumberEntry( 353.4, &o2o2_1150_1950_353p4[0], 2*sizeof(o2o2_1150_1950_353p4[0]), &o2o2_1150_1950_353p4[1], 2*sizeof(o2o2_1150_1950_353p4[0]),  N_ELEMENTS(o2o2_1150_1950_353p4)/2 );
	Set_Temperature						(293.0);
	SetInstrumentPSF_FWHMWavenumber		( 0.43);
	SetInstrumentPointSpacingWavenumber	( 1.0 );
	return ok;
}

/*---------------------------------------------------------------------------
 * skOpticalProperties_O4_HitranEntry_TempDependent::ConfigureAsRegion5 2020-01-02 */
/**
 *	Region 7 is from Thalman and Volkamer 2013. 
 *	Data points are 0.2449 cm-1 apart
 *  resolution is about 0.68 nm at 630 nm => 17.13 cm-1 
 *	Wavenumber range = 15290 - 16664 cm-1, 
 *	Wavelength range =   600 -   654 nm .

**/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_O4_HitranEntry_TempDependent::ConfigureAsRegion7()
{
	bool ok = true;

	ClearEntries();
	ok = ok && AddAscendingWavenumberEntry( 203.0, &o2o2_15290_16664_203p0[0], 2*sizeof(o2o2_15290_16664_203p0[0]), &o2o2_15290_16664_203p0[1], 2*sizeof(o2o2_15290_16664_203p0[0]),  N_ELEMENTS(o2o2_15290_16664_203p0)/2 );
	ok = ok && AddAscendingWavenumberEntry( 233.0, &o2o2_15290_16664_233p0[0], 2*sizeof(o2o2_15290_16664_233p0[0]), &o2o2_15290_16664_233p0[1], 2*sizeof(o2o2_15290_16664_233p0[0]),  N_ELEMENTS(o2o2_15290_16664_233p0)/2 );
	ok = ok && AddAscendingWavenumberEntry( 253.0, &o2o2_15290_16664_253p0[0], 2*sizeof(o2o2_15290_16664_253p0[0]), &o2o2_15290_16664_253p0[1], 2*sizeof(o2o2_15290_16664_253p0[0]),  N_ELEMENTS(o2o2_15290_16664_253p0)/2 );
	ok = ok && AddAscendingWavenumberEntry( 287.0, &o2o2_15290_16664_287p0[0], 2*sizeof(o2o2_15290_16664_287p0[0]), &o2o2_15290_16664_287p0[1], 2*sizeof(o2o2_15290_16664_287p0[0]),  N_ELEMENTS(o2o2_15290_16664_287p0)/2 );

	Set_Temperature						(293.0);
	SetInstrumentPSF_FWHMWavenumber		(17.13);
	SetInstrumentPointSpacingWavenumber	(1.0 );
	return ok;
}


/*---------------------------------------------------------------------------
 * skOpticalProperties_O4_HitranEntry_TempDependent::ConfigureAsRegion8 2020-01-02 */
/**
 *	Region 8 is from Thalman and Volkamer 2013.
 *	Data points are 1.0 cm-1 apart
 *  Resolution is about 0.37 nm at 450 nm => 18.27 cm-1 
 *	Wavenumber Range = 16700 - 29800 cm-1
 *  Wavelength Range =   335 - 599 nm
**/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_O4_HitranEntry_TempDependent::ConfigureAsRegion8()
{
	bool ok = true;

	ClearEntries();
	ok = ok && AddAscendingWavenumberEntry( 203.0, &o2o2_16658_29748_203p0[0], 2*sizeof(o2o2_16658_29748_203p0[0]), &o2o2_16658_29748_203p0[1], 2*sizeof(o2o2_16658_29748_203p0[0]),  N_ELEMENTS(o2o2_16658_29748_203p0)/2 );
	ok = ok && AddAscendingWavenumberEntry( 233.0, &o2o2_16780_29757_233p0[0], 2*sizeof(o2o2_16780_29757_233p0[0]), &o2o2_16780_29757_233p0[1], 2*sizeof(o2o2_16780_29757_233p0[0]),  N_ELEMENTS(o2o2_16780_29757_233p0)/2 );
	ok = ok && AddAscendingWavenumberEntry( 253.0, &o2o2_16791_29837_253p0[0], 2*sizeof(o2o2_16791_29837_253p0[0]), &o2o2_16791_29837_253p0[1], 2*sizeof(o2o2_16791_29837_253p0[0]),  N_ELEMENTS(o2o2_16791_29837_253p0)/2 );
	ok = ok && AddAscendingWavenumberEntry( 273.0, &o2o2_16668_29802_273p0[0], 2*sizeof(o2o2_16668_29802_273p0[0]), &o2o2_16668_29802_273p0[1], 2*sizeof(o2o2_16668_29802_273p0[0]),  N_ELEMENTS(o2o2_16668_29802_273p0)/2 );
	ok = ok && AddAscendingWavenumberEntry( 293.0, &o2o2_16645_29784_293p0[0], 2*sizeof(o2o2_16645_29784_293p0[0]), &o2o2_16645_29784_293p0[1], 2*sizeof(o2o2_16645_29784_293p0[0]),  N_ELEMENTS(o2o2_16645_29784_293p0)/2 );
	
	Set_Temperature						(293.0);
	SetInstrumentPSF_FWHMWavenumber		(18.27);
	SetInstrumentPointSpacingWavenumber	( 1.0 );
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O4_Hermans2000::skOpticalProperties_O4_Hermans2000		2009-11-6*/
 /** **/
 /*---------------------------------------------------------------------------*/

skOpticalProperties_O4_Hitran2016::skOpticalProperties_O4_Hitran2016()
{
	AddEntry(  1150.0,  1950.0, new skOpticalProperties_O4_HitranEntry_TempDependent(1) );
	AddEntry(  7450.0,  8491.0, new skOpticalProperties_O4_HitranEntry_NotDependent (2) );
	AddEntry(  9091.0,  9596.0, new skOpticalProperties_O4_HitranEntry_NotDependent (3) );
	AddEntry( 10512.0, 11228.0, new skOpticalProperties_O4_HitranEntry_NotDependent (4) );
	AddEntry( 12600.0, 13839.0, new skOpticalProperties_O4_HitranEntry_NotDependent (5) );
	AddEntry( 14206.0, 14898.0, new skOpticalProperties_O4_HitranEntry_NotDependent (6) );
	AddEntry( 15290.0, 16664.0, new skOpticalProperties_O4_HitranEntry_TempDependent(7) );
	AddEntry( 16700.0, 29800.0, new skOpticalProperties_O4_HitranEntry_TempDependent(8) );
}


/*---------------------------------------------------------------------------
 * skOpticalProperties_O4_Hitran2016::~skOpticalProperties_O4_Hitran20162020-01-02 */
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_O4_Hitran2016::~skOpticalProperties_O4_Hitran2016()
{
}

