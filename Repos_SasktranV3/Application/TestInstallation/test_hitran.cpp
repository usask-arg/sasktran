#include <sasktran.h>



/*---------------------------------------------------------------------------
 *                    make_pre_set_wavenumbers                    2019-11-08 */
/** **/
/*---------------------------------------------------------------------------*/
static void  make_pre_set_wavenumbers( nx1dArray<double>* wavenum, double* min_wavenum, double* max_wavenum)
{
	size_t n =  size_t( (1680-1640)/0.01 +0.5);			// Set wavenumbers from 1640 nm to 1680 nm
	wavenum->SetSize( n );
	for (size_t i = 0; i < n; i++)
	{
		wavenum->at(i) = 1.0E7/(1640 + i*0.01);
	}
	*min_wavenum = wavenum->back();
	*max_wavenum = wavenum->front();
}

/*---------------------------------------------------------------------------
 *                          test_hitran                           2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/

bool test_hitran()
{
	double				min_wavenum;
    double				max_wavenum;
	double				absxs;
	double				extxs;
	double				scattxs;
	bool				xschanged;
	nx1dArray<double>	wavenum;
	bool				ok;
	FILE*				f = stdout; //fopen("xs.txt", "wt");

	skOpticalProperties_HitranChemical*		ch4_opticalprops;
	skClimatology_MSIS90*					msis90;														// MSIS90 atmospheric model to set rayleigh atmosphere
	GEODETIC_INSTANT						pt( 52.1315, -106.6335, 10.0, 	58290.75000000);			// Set a point close to the ground, at our office in Saskatoon

	fprintf(f,"Starting Test of HITRAN. This will briefly test the cross-section code near to the ground.\n");
	fprintf(f,"     (i)   Spectral Line Caches\n");
	fprintf(f,"     (ii)  Internal Partition Caches\n");
	fprintf(f,"     (iii) Pre-set wavelength optimization for the HR Engine\n");
	fprintf(f,"     NOTE THE CODE MAY BE SLOW THE FIRST TIME THROUGH IF IT HAS TO CREATE THE CACHES\n");

	make_pre_set_wavenumbers( &wavenum, &min_wavenum, &max_wavenum);

    ch4_opticalprops = new skOpticalProperties_HitranChemical("CO2");
	msis90           = new skClimatology_MSIS90;
	ch4_opticalprops->AddRef();
	ch4_opticalprops->SetWavenumberRange(min_wavenum, max_wavenum);
	msis90->AddRef();

	ok =       ch4_opticalprops->EnableCachedCrossSections( wavenum.UnsafeArrayBasePtr(), wavenum.size() );
	ok = ok && ch4_opticalprops->SetAtmosphericState(msis90);
	ok = ok && ch4_opticalprops->SetLocation(pt, &xschanged);

	for (size_t i = 0; i < wavenum.size(); i++)
	{
		ch4_opticalprops->CalculateCrossSections( wavenum[i], &absxs, &extxs, &scattxs );
	}
	ch4_opticalprops->Release();
	msis90->Release();
	fprintf(f,"Finished Hitran test. We just calculated %d cross-sections\n", (int)wavenum.size());
	return ok;
}