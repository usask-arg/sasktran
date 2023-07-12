#include <sasktran.h>



/*---------------------------------------------------------------------------
 *                    make_pre_set_wavenumbers                    2019-11-08 */
/** **/
/*---------------------------------------------------------------------------*/

static void  make_pre_set_wavenumbers( double startnm, double endnm, nx1dArray<double>* wavenum, double* min_wavenum, double* max_wavenum)
{
	size_t n =  size_t( (endnm-startnm) +0.5);			// Set wavenumbers from 1640 nm to 1680 nm
	wavenum->SetSize( n );
	for (size_t i = 0; i < n; i++)
	{
		wavenum->at(i) = 1.0E7/(startnm + i);
	}
	*min_wavenum = wavenum->back();
	*max_wavenum = wavenum->front();
}

/*---------------------------------------------------------------------------
 *                        test_fally_o2_o2                        2019-12-03 */
/** **/
/*---------------------------------------------------------------------------*/

static bool test_fally_o2_o2()
{
	double				min_wavenum;
    double				max_wavenum;
	double				absxs;
	double				extxs;
	double				scattxs;
	bool				xschanged;
	nx1dArray<double>	wavenum;
	bool				ok;
	bool				ok1;
	FILE*				f = stdout; //fopen("xs.txt", "wt");

	skOpticalProperties_O4_Fally2000*		o2_o2_opticalprops;
	skClimatology_MSIS90*					msis90;														// MSIS90 atmospheric model to set rayleigh atmosphere
	GEODETIC_INSTANT						pt( 0.0, 0.0, 25000.0, 	52393.3792987115);			// Set a point close to the ground, at our office in Saskatoon

	//fprintf(f,"Testing FReeman 1984 SO2.\n");

	make_pre_set_wavenumbers( 330, 670, &wavenum, &min_wavenum, &max_wavenum);

    o2_o2_opticalprops = new skOpticalProperties_O4_Fally2000(); 
	msis90           = new skClimatology_MSIS90;
	o2_o2_opticalprops->AddRef();
	msis90->AddRef();

	ok = true;
	ok = ok && o2_o2_opticalprops->SetAtmosphericState(msis90);
	ok = ok && o2_o2_opticalprops->SetLocation(pt, &xschanged);

	for (size_t i = 0; i < wavenum.size(); i++)
	{
		ok1 = o2_o2_opticalprops->CalculateCrossSections( wavenum[i], &absxs, &extxs, &scattxs );
		ok = ok && ok1;
	}
	o2_o2_opticalprops->Release();
	msis90->Release();
	fprintf(f,"O2/O2 Status [%s]: Fally 2000 test. Calculated %d cross-sections.\n", (const char*)(ok ? "OK" : "FAIL"), (int)wavenum.size()  );
	return ok;
}


/*---------------------------------------------------------------------------
 *                     test_hitran2016_o2_o2                      2020-01-02 */
/** **/
/*---------------------------------------------------------------------------*/

static bool test_hitran2016_o2_o2()
{
	double				min_wavenum;
    double				max_wavenum;
	double				absxs;
	double				extxs;
	double				scattxs;
	bool				xschanged;
	nx1dArray<double>	wavenum;
	bool				ok;
	bool				ok1;
	FILE*				f = stdout; //fopen("xs.txt", "wt");

	skOpticalProperties_O4_Hitran2016*		o2_o2_opticalprops;
	skClimatology_MSIS90*					msis90;														// MSIS90 atmospheric model to set rayleigh atmosphere
	GEODETIC_INSTANT						pt( 0.0, 0.0, 25000.0, 	52393.3792987115);			// Set a point close to the ground, at our office in Saskatoon

	//fprintf(f,"Testing FReeman 1984 SO2.\n");

	make_pre_set_wavenumbers( 336, 8000, &wavenum, &min_wavenum, &max_wavenum);

    o2_o2_opticalprops = new skOpticalProperties_O4_Hitran2016(); 
	msis90           = new skClimatology_MSIS90;
	o2_o2_opticalprops->AddRef();
	msis90->AddRef();

	ok = true;
	ok = ok && o2_o2_opticalprops->SetAtmosphericState(msis90);
	ok = ok && o2_o2_opticalprops->SetLocation(pt, &xschanged);

	for (size_t i = 0; i < wavenum.size(); i++)
	{
		double wv = wavenum[i];
		ok1 = o2_o2_opticalprops->CalculateCrossSections( wv, &absxs, &extxs, &scattxs );
		ok = ok && ok1;
		//fprintf(f, "%10.3f  %20.7e\n", (double)wv, (double)absxs);
	}
	o2_o2_opticalprops->Release();
	msis90->Release();
	fprintf(f,"O2/O2 Status [%s]: Hitran 2016 test. Calculated %d cross-sections.\n", (const char*)(ok ? "OK" : "FAIL"), (int)wavenum.size()  );
	return ok;
}


/*---------------------------------------------------------------------------
 *                        test_fally_o2_o2                        2019-12-03 */
/** **/
/*---------------------------------------------------------------------------*/
static bool test_thalman_o2_o2()
{
	double				min_wavenum;
    double				max_wavenum;
	double				absxs;
	double				extxs;
	double				scattxs;
	bool				xschanged;
	nx1dArray<double>	wavenum;
	bool				ok;
	bool				ok1;
	FILE*				f = stdout; //fopen("xs.txt", "wt");

	skOpticalProperties_O4_Thalman2013*		o2_o2_opticalprops;
	skClimatology_MSIS90*					msis90;														// MSIS90 atmospheric model to set rayleigh atmosphere
	GEODETIC_INSTANT						pt( 0.0, 0.0, 25000.0, 	52393.3792987115);			// Set a point close to the ground, at our office in Saskatoon

	//fprintf(f,"Testing FReeman 1984 SO2.\n");

	make_pre_set_wavenumbers( 330, 670, &wavenum, &min_wavenum, &max_wavenum);

    o2_o2_opticalprops = new skOpticalProperties_O4_Thalman2013(); 
	msis90           = new skClimatology_MSIS90;
	o2_o2_opticalprops->AddRef();
	msis90->AddRef();

	ok = true;
	ok = ok && o2_o2_opticalprops->SetAtmosphericState(msis90);
	ok = ok && o2_o2_opticalprops->SetLocation(pt, &xschanged);

	for (size_t i = 0; i < wavenum.size(); i++)
	{
		ok1 = o2_o2_opticalprops->CalculateCrossSections( wavenum[i], &absxs, &extxs, &scattxs );
		ok = ok && ok1;
	}
	o2_o2_opticalprops->Release();
	msis90->Release();
	fprintf(f,"O2/O2 Status [%s]: Thalman 2013 test. Calculated %d cross-sections.\n", (const char*)(ok ? "OK" : "FAIL"), (int)wavenum.size()  );
	return ok;
}


/*---------------------------------------------------------------------------
 *                          test_hitran                           2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/

bool test_o2_o2()
{
	bool ok1;
	bool ok2;

	ok1 = test_hitran2016_o2_o2();
	ok1 = test_fally_o2_o2();
	ok2 = test_thalman_o2_o2();
	return ok1;
}
