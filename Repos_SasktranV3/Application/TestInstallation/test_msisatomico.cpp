#include <sasktran.h>

 bool test_msisatomicO()
 {
	skClimatology_MSIS90*			msis     = new skClimatology_MSIS90();
	double						     mjd      = 52393.3792987115;
	double							 f[2]     =  {0.0, 1.0};
	double                           h[2]     =  {10000.0, 15000.0};
	nx1dArray<double>				 farray(2, f);
	nx1dArray<double>				 harray(2, h);
	bool							 ok;



	nx1dArray<double>		heightm;
	nx1dArray<double>		v;
	GEODETIC_INSTANT		pt;
	size_t					numbad;

	msis->AddRef();
	msis->CachedMSIS().AddSpecies(SKCLIMATOLOGY_O_CM3);
	msis->CachedMSIS().SetMaxHeightKMS(250.0);

	pt.latitude = 52.0;
	pt.longitude = 106.0;
	pt.heightm   = 0.0;
	pt.mjd       = mjd;

	heightm.Indgen(1000);
	heightm *= 250.0;
	v.SetSize( heightm.size() );

	ok = msis->GetHeightProfile( SKCLIMATOLOGY_O_CM3, pt, heightm.ArrayBasePtr(), (int)heightm.size(), v.UnsafeArrayBasePtr(), true, &numbad );
	
	msis->Release();
	return ok;
 }
