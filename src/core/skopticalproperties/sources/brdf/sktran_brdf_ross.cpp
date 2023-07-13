#include <skopticalproperties21.h>
#include <sources/sasktranif_opticalimpl/skbrdf_stubs.h>



/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Ross::ISKBrdf_Stub_Ross		 2017-03-08      */
/** **/
/*---------------------------------------------------------------------------*/

ISKBrdf_Stub_Ross_Kernel::ISKBrdf_Stub_Ross_Kernel(SKTRAN_BRDF_Ross_Kernel* rossbrdf)
	: ISKBrdf_Stub_Base(rossbrdf)
{
	m_rossbrdf = rossbrdf;
}


/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Ross::~ISKBrdf_Stub_Ross		 2017-03-08      */
/** **/
/*---------------------------------------------------------------------------*/

ISKBrdf_Stub_Ross_Kernel::~ISKBrdf_Stub_Ross_Kernel(){}

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Ross::SKTRAN_BRDF_Ross()		 2017-03-08*/
/** **/
/*---------------------------------------------------------------------------*/
SKTRAN_BRDF_Ross_Kernel::SKTRAN_BRDF_Ross_Kernel(){}


/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Ross::cos_primed_scattering_angle()		 2017-03-08*/
/** **/
/*---------------------------------------------------------------------------*/
double SKTRAN_BRDF_Ross_Kernel::cos_primed_scattering_angle(double cs, double cv, double cp) const
{
	double ss = sqrt(std::max(0.0, 1.0 - cs*cs));
	double sv = sqrt(std::max(0.0, 1.0 - cv*cv));
	double ceta = cs*cv + ss*sv*cp;
	return ceta;
}



/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Ross_Thick::BRDF()		 2017-03-08*/
/** **/
/*---------------------------------------------------------------------------*/
bool SKTRAN_BRDF_RossThick_Kernel::BRDF(double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* brdf) const
{
	CheckCosines(MU_in, MU_out, COSDPHI, "SKTRAN_BRDF_RossThick_Kernel::BRDF");

	double cs = MU_in;
	double cv = MU_out;
	double cp = -COSDPHI; // take the negative to match Sasktran's azimuth convention

	double ceta = std::max(-1.0, std::min(1.0, SKTRAN_BRDF_Ross_Kernel::cos_primed_scattering_angle(cs, cv, cp)));
	double seta = sqrt(std::max(0.0, 1.0 - ceta*ceta));
	double eta = nxmath::DegreesToRadians(nxmath::acosd(ceta));
	
	*brdf = (((nxmath::PiOver2 - eta) * ceta + seta) / (cs + cv)) / nxmath::Pi - 0.25; // scale down by pi to match sasktran
	return NXFINITE(*brdf);
}

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Ross_Thin::BRDF()		 2017-03-08*/
/** **/
/*---------------------------------------------------------------------------*/
bool SKTRAN_BRDF_RossThin_Kernel::BRDF(double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* brdf) const
{
	CheckCosines(MU_in, MU_out, COSDPHI, "SKTRAN_BRDF_RossThin_Kernel::BRDF");

	double cs = MU_in;
	double cv = MU_out;
	double cp = -COSDPHI; // take the negative to match Sasktran's azimuth convention

	double ceta = std::max(-1.0, std::min(1.0, SKTRAN_BRDF_Ross_Kernel::cos_primed_scattering_angle(cs, cv, cp)));
	double seta = sqrt(std::max(0.0, 1.0 - ceta*ceta));
	double eta = nxmath::DegreesToRadians(nxmath::acosd(ceta));

	*brdf = ((nxmath::PiOver2 - eta) * ceta + seta) / (cs * cv * nxmath::Pi) - 0.5; // scale down by pi to match sasktran
	return NXFINITE(*brdf);
}