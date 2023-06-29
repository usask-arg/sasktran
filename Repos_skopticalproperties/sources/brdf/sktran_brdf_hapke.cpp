#include <skopticalproperties21.h>
#include <sources/sasktranif_opticalimpl/skbrdf_stubs.h>



/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Hapke::ISKBrdf_Stub_Hapke		 2017-07-28*/
/** **/
/*---------------------------------------------------------------------------*/

ISKBrdf_Stub_Hapke::ISKBrdf_Stub_Hapke(SKTRAN_BRDF_Hapke* hapkebrdf)
	: ISKBrdf_Stub_Base(hapkebrdf)
{
	m_hapkebrdf = hapkebrdf;
}


/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Hapke::~ISKBrdf_Stub_Hapke		 2017-07-28*/
/** **/
/*---------------------------------------------------------------------------*/

ISKBrdf_Stub_Hapke::~ISKBrdf_Stub_Hapke()
{}

/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Hapke::SetPropertyArray		 2017-07-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKBrdf_Stub_Hapke::SetPropertyArray(const char* propertyname, const double* value, int numpoints)
{
	nxString					name(propertyname);
	bool						ok = m_hapkebrdf != nullptr;

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "ISKBrdf_Stub_Hapke::SetProperty, the internal C++ object is undefined. Thats not good");
	}
	else
	{
		if (name == "BRDFParameters")
		{
			ok = (numpoints == 3);
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "ISKBrdf_Stub_Hapke::SetProperty(BRDFParameters), accepts only 3 parameters. You passed in %d parameters", (int)numpoints);
			}
			else
			{
				double omega = value[0];
				double delta = value[1]; 
				double b0 = value[2]; 


				ok = m_hapkebrdf->SetBRDFParameters(omega, delta, b0);
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING, "ISKBrdf_Stub_Hapke::SetProperty(BRDFParameters), failed to set parameters to omega=%e, delta=%e, b0=%e", (double)omega, (double)delta, (double)b0);
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
*					SKTRAN_BRDF_Hapke::SKTRAN_BRDF_Hapke		 2017-07-28*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_BRDF_Hapke::SKTRAN_BRDF_Hapke()
{
	m_omega = std::numeric_limits<double>::quiet_NaN();
	m_delta = std::numeric_limits<double>::quiet_NaN();
	m_b0 = std::numeric_limits<double>::quiet_NaN();
}


/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Hapke::SetBRDFParameters		 2017-07-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_BRDF_Hapke::SetBRDFParameters(double omega, double delta, double b0)
{
	m_omega = omega;
	m_delta = delta;
	m_b0 = b0;
	return true;
}


/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Hapke::BRDF		 2017-07-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_BRDF_Hapke::BRDF(double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* brdf) const
{
	bool ok = NXFINITE(m_omega) && NXFINITE(m_delta) && NXFINITE(m_b0);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_BRDF_Hapke::BRDF, one or more of the 3 BRDF parameters is NaN. Cannot calculate BRDF until you define value values for 3 paraemeters");
		*brdf = std::numeric_limits<double>::quiet_NaN();
	}
	else
	{
		CheckCosines(MU_in, MU_out, COSDPHI, "SKTRAN_BRDF_Hapke::BRDF");

		double cs = MU_in;
		double cv = MU_out;
		double ss = sqrt(1 - cs * cs);
		double sv = sqrt(1 - cv * cv);
		double ts = ss / cs;
		double tv = sv / cv;
		double cp = -COSDPHI; // take the negative to match Sasktran's azimuth convention

		double ceta = std::max(-1.0, std::min(1.0, cs * cv + ss * sv * cp));
		double eta = acos(ceta);
		double gamma = sqrt(1 - m_omega);
		double Tv = (1 + 2 * cv) / (1 + 2 * cv * gamma);
		double Ts = (1 + 2 * cs) / (1 + 2 * cs * gamma);
		double B = (m_b0 * m_delta) / (m_delta + tan(0.5 * eta));
		double P = 1 + 0.5 * ceta;
		double Rir = m_omega / (4 * (cs + cv));
		*brdf = Rir * ((1 + B) * P + Ts*Tv - 1) / nxmath::Pi; // Spurr (A.10) - scale down by pi to match sasktran
	}
	return ok;
}

