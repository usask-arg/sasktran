#include <skopticalproperties21.h>
#include <sources/sasktranif_opticalimpl/skbrdf_stubs.h>



/*-----------------------------------------------------------------------------
 *					ISKBrdf_Stub_Snow_Kokhanovsky2012::ISKBrdf_Stub_Snow_Kokhanovsky2012		 2016- 12- 12*/
/** **/
/*---------------------------------------------------------------------------*/

ISKBrdf_Stub_Snow_Kokhanovsky2012::ISKBrdf_Stub_Snow_Kokhanovsky2012( SKTRAN_BRDF_Snow_Kokhanovsky2012* snowbrdf)
	                          : ISKBrdf_Stub_Base( snowbrdf)
{
	m_snowbrdf = snowbrdf;
}


/*-----------------------------------------------------------------------------
 *					ISKBrdf_Stub_Snow_Kokhanovsky2012::~ISKBrdf_Stub_Snow_Kokhanovsky2012		 2016- 12- 12*/
/** **/
/*---------------------------------------------------------------------------*/

ISKBrdf_Stub_Snow_Kokhanovsky2012::~ISKBrdf_Stub_Snow_Kokhanovsky2012()
{}

/*-----------------------------------------------------------------------------
 *					ISKBrdf_Stub_Snow_Kokhanovsky2012::SetPropertyArray		 2016- 12- 12*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKBrdf_Stub_Snow_Kokhanovsky2012::SetPropertyArray( const char* propertyname, const double* value, int numpoints)
{
	nxString					name( propertyname);
	bool						ok = m_snowbrdf != nullptr;

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"ISKBrdf_Stub_Snow_Kokhanovsky2012::SetProperty, the internal C++ object is undefined. Thats not good");
	}
	else
	{
		if (name == "BRDFParameters")
		{
			ok = (numpoints == 2);
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKBrdf_Stub_Snow_Kokhanovsky2012::SetProperty(BRDFParameters), accepts only 2 parameters. You passed in %d parameters", (int)numpoints);
			}
			else
			{
				double L  = value[0];
				double M  = value[1];
		
				ok = m_snowbrdf->SetBRDFParameters( L, M);
				if (!ok)
				{
					nxLog::Record( NXLOG_WARNING,"ISKBrdf_Stub_Snow_Kokhanovsky2012::SetProperty(BRDFParameters), failed to set parameters to L=%e, M=%e", (double)L, (double)M);
				}
			}
		}
		else
		{
			ok = ISKBrdf_Stub_Base::SetPropertyArray( propertyname,value, numpoints);
		}
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_BRDF_Snow_Kokhanovsky2012::SKTRAN_BRDF_Snow_Kokhanovsky2012		 2016- 12- 12*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_BRDF_Snow_Kokhanovsky2012::SKTRAN_BRDF_Snow_Kokhanovsky2012()
{
	m_L = 3.6;				// Default values taken from figure 1 of Kokhanovsky et al. 2012, derived from a site in central Greenland in May 2006
	m_M = 5.5E-08;			// Kokhanovsky recommnds that these are 2 free parameters which should be fitted from the data observations
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_BRDF_Snow_Kokhanovsky2012::K0		 2016- 12- 12*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_BRDF_Snow_Kokhanovsky2012::K0( double mu) const
{
	return (3.0/7.0)*( 1.0 + 2.0*mu);
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_BRDF_Snow_Kokhanovsky2012::SetBRDFParameters		 2016- 12- 12*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_BRDF_Snow_Kokhanovsky2012::SetBRDFParameters( double L, double M )
{
	m_L = L;
	m_M = M;
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_BRDF_Snow_Kokhanovsky2012::p		 2016- 12- 12*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_BRDF_Snow_Kokhanovsky2012::p( double thetadegrees) const
{
	return (11.1*exp( -0.087*thetadegrees) + 1.1*exp( -0.014*thetadegrees));
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_BRDF_Snow_Kokhanovsky2012::R0		 2016- 12- 12*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_BRDF_Snow_Kokhanovsky2012::R0( double mus, double muv, double theta) const
{
	const double    a = 1.247;
	const double	b = 1.186;
	const double	c = 5.157;

	return ( a + b*(mus+muv) + c*mus*muv + p(theta))/(4.0*(mus+muv));
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_BRDF_Snow_Kokhanovsky2012::Chi		 2016- 12- 12*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_BRDF_Snow_Kokhanovsky2012::Chi( double wavelen_nm ) const
{
	return m_iceri.RefractiveIndex(1.0E07/wavelen_nm).imag();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_BRDF_Snow_Kokhanovsky2012::BRDF		 2016- 12- 12*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_BRDF_Snow_Kokhanovsky2012::BRDF( double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* brdf) const
{
	CheckCosines(MU_in, MU_out, COSDPHI, "SKTRAN_BRDF_Snow_Kokhanovsky2012::BRDF");

	double			mus = MU_in;
	double			muv = MU_out;
	double			ss  = sqrt( 1 - mus*mus);
	double			sv  = sqrt( 1 - muv*muv);
	double			cost  = std::max( -1.0, std::min( 1.0, -mus*muv + ss*sv*COSDPHI) );
	double			theta = nxmath::acosd(cost);
	double			alpha;
	double			gamma;
	double			r0;

	gamma  = 4 * nxmath::Pi * (Chi(wavelennm) + m_M) / wavelennm;
	alpha  = sqrt(gamma * m_L);
	r0     = R0(mus, muv, theta);
	*brdf = r0 * (exp(-alpha * K0(mus) * K0(muv) / r0)) / nxmath::Pi; // scaled down by pi to match sasktran
	return true;
}

