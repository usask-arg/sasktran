#include <sasktran.h>

 bool test_linearcombo()
 {
	skClimatology_LabowOzoneVMR*    labow    = new skClimatology_LabowOzoneVMR();
	skClimatology_Constant*		    constant = new skClimatology_Constant();
	skClimatologyLinearCombination* climate  = new skClimatologyLinearCombination();
	double						     mjd      = 52393.3792987115;
	double							 f[2]     =  {0.0, 1.0};
	double                           h[2]     =  {10000.0, 15000.0};
	nx1dArray<double>				 farray(2, f);
	nx1dArray<double>				 harray(2, h);
	bool							 ok;

	constant->SetConstantValue( 8.0E15);
	climate->SetFirstClimatology(labow);
	climate->SetSecondClimatology(constant);
	climate->SetHeightProfileCoeffsOfFirstProperty(harray, farray);

	nx1dArray<double>		heightm;
	nx1dArray<double>		v;
	GEODETIC_INSTANT		pt;
	size_t					numbad;


	pt.latitude = 52.0;
	pt.longitude = 106.0;
	pt.heightm   = 0.0;
	pt.mjd       = mjd;

	heightm.Indgen(160);
	heightm *= 250.0;
	v.SetSize( heightm.size() );

	ok = climate->GetHeightProfile( SKCLIMATOLOGY_O3_CM3, pt, heightm.ArrayBasePtr(), (int)heightm.size(), v.UnsafeArrayBasePtr(), true, &numbad );
	return ok;
 }
