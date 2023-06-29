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
 *                        test_freeman_so2                        2019-12-02 */
/** **/
/*---------------------------------------------------------------------------*/
static bool test_freeman_so2()
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

	skOpticalProperties_SO2_Freeman1984*	so2_opticalprops;
	skClimatology_MSIS90*					msis90;														// MSIS90 atmospheric model to set rayleigh atmosphere
	GEODETIC_INSTANT						pt( 0.0, 0.0, 25000.0, 	52393.3792987115);			// Set a point close to the ground, at our office in Saskatoon

	//fprintf(f,"Testing FReeman 1984 SO2.\n");

	make_pre_set_wavenumbers( 170, 240, &wavenum, &min_wavenum, &max_wavenum);

    so2_opticalprops = new skOpticalProperties_SO2_Freeman1984(); 
	msis90           = new skClimatology_MSIS90;
	so2_opticalprops->AddRef();
	msis90->AddRef();

	ok = true;
	ok = ok && so2_opticalprops->SetAtmosphericState(msis90);
	ok = ok && so2_opticalprops->SetLocation(pt, &xschanged);

	for (size_t i = 0; i < wavenum.size(); i++)
	{
		ok1 = so2_opticalprops->CalculateCrossSections( wavenum[i], &absxs, &extxs, &scattxs );
		ok = ok && ok1;
	}
	so2_opticalprops->Release();
	msis90->Release();
	fprintf(f,"SO2 Status [%s]: Freeman 1984 test. Calculated %d cross-sections.\n", (const char*)(ok ? "OK" : "FAIL"), (int)wavenum.size()  );
	return ok;
}


/*---------------------------------------------------------------------------
 *                          test_hitran                           2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/

static bool test_bogumil_so2()
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

	skOpticalProperties_SO2_Bogumil2003*	so2_opticalprops;
	skClimatology_MSIS90*					msis90;														// MSIS90 atmospheric model to set rayleigh atmosphere
	GEODETIC_INSTANT						pt( 0.0, 0.0, 25000.0, 	52393.3792987115);			// Set a point close to the ground, at our office in Saskatoon

	//fprintf(f,"Testing Bogumil 2003 SO2.\n");

	make_pre_set_wavenumbers( 235, 400, &wavenum, &min_wavenum, &max_wavenum);

    so2_opticalprops = new skOpticalProperties_SO2_Bogumil2003(); 
	msis90           = new skClimatology_MSIS90;
	so2_opticalprops->AddRef();
	msis90->AddRef();

	ok = true;
	ok = ok && so2_opticalprops->SetAtmosphericState(msis90);
	ok = ok && so2_opticalprops->SetLocation(pt, &xschanged);

	for (size_t i = 0; i < wavenum.size(); i++)
	{
		ok1 = so2_opticalprops->CalculateCrossSections( wavenum[i], &absxs, &extxs, &scattxs );
		ok = ok && ok1;
	}
	so2_opticalprops->Release();
	msis90->Release();
	fprintf(f,"SO2 Status [%s]: Bogumil 2003 test. Calculated %d cross-sections.\n", (const char*)(ok ? "OK" : "FAIL"), (int)wavenum.size() );
	return ok;
}




/*---------------------------------------------------------------------------
 *                         test_rufus_so2                         2019-12-02 */
/** **/
/*---------------------------------------------------------------------------*/
static bool test_rufus_so2()
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

	skOpticalProperties_SO2_Rufus2003*	so2_opticalprops;
	skClimatology_MSIS90*					msis90;														// MSIS90 atmospheric model to set rayleigh atmosphere
	GEODETIC_INSTANT						pt( 0.0, 0.0, 25000.0, 	52393.3792987115);			// Set a point close to the ground, at our office in Saskatoon

	//fprintf(f,"Testing Rufus 2003 SO2.\n");

	make_pre_set_wavenumbers( 210, 350, &wavenum, &min_wavenum, &max_wavenum);

    so2_opticalprops = new skOpticalProperties_SO2_Rufus2003(); 
	msis90           = new skClimatology_MSIS90;
	so2_opticalprops->AddRef();
	msis90->AddRef();

	ok = true;
	ok = ok && so2_opticalprops->SetAtmosphericState(msis90);
	ok = ok && so2_opticalprops->SetLocation(pt, &xschanged);

	for (size_t i = 0; i < wavenum.size(); i++)
	{
		ok1 = so2_opticalprops->CalculateCrossSections( wavenum[i], &absxs, &extxs, &scattxs );
		ok = ok && ok1;
	}
	so2_opticalprops->Release();
	msis90->Release();
	fprintf(f,"SO2 Status [%s]: Rufus 2003 test. Calculated %d cross-sections.\n", (const char*)(ok ? "OK" : "FAIL"), (int)wavenum.size() );
	return ok;
}

/*---------------------------------------------------------------------------
 *                         test_rufus_so2                         2019-12-02 */
/** **/
/*---------------------------------------------------------------------------*/
static bool test_vandaele2009_so2()
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

	skOpticalProperties_SO2_Vandaele2009*	so2_opticalprops;
	skClimatology_MSIS90*					msis90;														// MSIS90 atmospheric model to set rayleigh atmosphere
	GEODETIC_INSTANT						pt( 0.0, 0.0, 25000.0, 	52393.3792987115);			// Set a point close to the ground, at our office in Saskatoon

	//fprintf(f,"Testing Vandaele 2009 SO2.\n");

	make_pre_set_wavenumbers( 225, 420, &wavenum, &min_wavenum, &max_wavenum);

    so2_opticalprops = new skOpticalProperties_SO2_Vandaele2009(); 
	msis90           = new skClimatology_MSIS90;
	so2_opticalprops->AddRef();
	msis90->AddRef();

	ok = true;
	ok = ok && so2_opticalprops->SetAtmosphericState(msis90);
	ok = ok && so2_opticalprops->SetLocation(pt, &xschanged);

	for (size_t i = 0; i < wavenum.size(); i++)
	{
		ok1 = so2_opticalprops->CalculateCrossSections( wavenum[i], &absxs, &extxs, &scattxs );
		ok = ok && ok1;
	}
	so2_opticalprops->Release();
	msis90->Release();
	fprintf(f,"SO2 Status [%s]: Vandaele 2009 test. Calculated %d cross-sections. \n", (const char*)(ok ? "OK" : "FAIL") , (int)wavenum.size());
	return ok;
}

/*---------------------------------------------------------------------------
 *                          test_hitran                           2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/

bool test_so2()
{
	bool ok1;
	bool ok2;
	bool ok3;
	bool ok4;

	ok1 = test_freeman_so2();
	ok2 = test_bogumil_so2();
	ok3 = test_rufus_so2();
	ok4 = test_vandaele2009_so2();
	return ok1 && ok2;
}
