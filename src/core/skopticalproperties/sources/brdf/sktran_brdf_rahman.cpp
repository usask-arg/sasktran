#include <skopticalproperties21.h>
#include <sources/sasktranif_opticalimpl/skbrdf_stubs.h>



/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Rahman::ISKBrdf_Stub_Rahman		 2017-07-28*/
/** **/
/*---------------------------------------------------------------------------*/

ISKBrdf_Stub_Rahman::ISKBrdf_Stub_Rahman(SKTRAN_BRDF_Rahman* rahmanbrdf)
	: ISKBrdf_Stub_Base(rahmanbrdf)
{
	m_rahmanbrdf = rahmanbrdf;
}


/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Rahman::~ISKBrdf_Stub_Rahman		 2017-07-28*/
/** **/
/*---------------------------------------------------------------------------*/

ISKBrdf_Stub_Rahman::~ISKBrdf_Stub_Rahman()
{}

/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Rahman::SetPropertyArray		 2017-07-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKBrdf_Stub_Rahman::SetPropertyArray(const char* propertyname, const double* value, int numpoints)
{
	nxString					name(propertyname);
	bool						ok = m_rahmanbrdf != nullptr;

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "ISKBrdf_Stub_Rahman::SetProperty, the internal C++ object is undefined. Thats not good");
	}
	else
	{
		if (name == "BRDFParameters")
		{
			ok = (numpoints == 3);
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "ISKBrdf_Stub_Rahman::SetProperty(BRDFParameters), accepts only 3 parameters. You passed in %d parameters", (int)numpoints);
			}
			else
			{
				double rho0 = value[0]; // characterizes overall intensity of reflection (but isn't single scatter albedo) - rho0 >= 0
				double theta = value[1]; // characerizes the amount of forward (0 <= theta <= 1) or backward (-1 <= theta <= 0) scattering
				double k = value[2]; // characterizes anisotropy in the surface
				

				ok = m_rahmanbrdf->SetBRDFParameters(rho0, theta, k);
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING, "ISKBrdf_Stub_Rahman::SetProperty(BRDFParameters), failed to set parameters to rho0=%e, theta=%e, k=%e", (double)rho0, (double)theta, (double)k);
				}
			}
		}
		else
		{
			ok = ISKBrdf_Stub_Base::SetPropertyArray(propertyname, value, numpoints);
		}
	}
	return ok;
}



/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Rahman::SKTRAN_BRDF_Rahman		 2017-07-28*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_BRDF_Rahman::SKTRAN_BRDF_Rahman()
{
	m_rho0 = std::numeric_limits<double>::quiet_NaN();
	m_theta = std::numeric_limits<double>::quiet_NaN();
	m_k = std::numeric_limits<double>::quiet_NaN();
}


/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Rahman::SetBRDFParameters		 2017-07-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_BRDF_Rahman::SetBRDFParameters(double rho0, double theta, double k)
{
	m_rho0 = rho0;
	m_theta = theta;
	m_k = k;
	return true;
}


/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Rahman::BRDF		 2017-07-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_BRDF_Rahman::BRDF(double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* brdf) const
{
	bool ok = NXFINITE(m_rho0) && NXFINITE(m_theta) && NXFINITE(m_k);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_BRDF_Rahman::BRDF, one or more of the 3 BRDF parameters is NaN. Cannot calculate BRDF until you define value values for 3 paraemeters");
		*brdf = std::numeric_limits<double>::quiet_NaN();
	}
	else
	{
		CheckCosines(MU_in, MU_out, COSDPHI, "SKTRAN_BRDF_Rahman::BRDF");

		// basic sines and cosines
		double cs = MU_in;
		double cv = MU_out;
		double ss = sqrt(1 - cs * cs);
		double sv = sqrt(1 - cv * cv);
		double ts = ss / cs;
		double tv = sv / cv;
		double cp = -COSDPHI; // take the negative to match Sasktran's azimuth convention

		double ceta = cs * cv + ss * sv * cp; // Rahman eqn (4)
		double G = sqrt(ts * ts + tv * tv - 2 * ts * tv * cp); // Rahman eqn (6)
		double R = (1.0 - m_rho0) / (1.0 + G); // Rahman eqn (5)
		double F = (1.0 - m_theta * m_theta) / pow(1.0 + m_theta * m_theta + 2.0 * m_theta * ceta, 1.5); // Rahman eqn (3)
		*brdf = (m_rho0 * F * (1.0 + R) * pow(cs * cv * (cs + cv), m_k - 1.0)) / nxmath::Pi; // Rahman eqn (2) - scale down by pi to match sasktran
	}
	return ok;
}

