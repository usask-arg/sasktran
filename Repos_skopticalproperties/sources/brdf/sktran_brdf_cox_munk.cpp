#include <skopticalproperties21.h>
#include <sources/sasktranif_opticalimpl/skbrdf_stubs.h>



/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Cox_Munk::ISKBrdf_Stub_Cox_Munk		 2017-07-27*/
/** **/
/*---------------------------------------------------------------------------*/

ISKBrdf_Stub_Cox_Munk::ISKBrdf_Stub_Cox_Munk(SKTRAN_BRDF_CoxMunk* waterbrdf)
	: ISKBrdf_Stub_Base(waterbrdf)
{
	m_waterbrdf = waterbrdf;
}


/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Cox_Munk::~ISKBrdf_Stub_Cox_Munk		 2017-07-27*/
/** **/
/*---------------------------------------------------------------------------*/

ISKBrdf_Stub_Cox_Munk::~ISKBrdf_Stub_Cox_Munk()
{}

/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Cox_Munk::SetPropertyArray		 2017-07-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKBrdf_Stub_Cox_Munk::SetPropertyArray(const char* propertyname, const double* value, int numpoints)
{
	nxString					name(propertyname);
	bool						ok = m_waterbrdf != nullptr;

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "ISKBrdf_Stub_Cox_Munk::SetProperty, the internal C++ object is undefined. Thats not good");
	}
	else
	{
		if (name == "BRDFParameters")
		{
			ok = (numpoints == 2);
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "ISKBrdf_Stub_Cox_Munk::SetProperty(BRDFParameters), accepts only 2 parameters. You passed in %d parameters", (int)numpoints);
			}
			else
			{
				double w = value[0]; // wind speed, m/s
				double n = value[1]; // water refractive index

				ok = m_waterbrdf->SetBRDFParameters(w, n);
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING, "ISKBrdf_Stub_Cox_Munk::SetProperty(BRDFParameters), failed to set parameters to wind_speed=%e, index_of_refraction=%e", (double)w, (double)n);
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
*					SKTRAN_BRDF_Cox_Munk::SKTRAN_BRDF_Cox_Munk		 2017-07-27*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_BRDF_CoxMunk::SKTRAN_BRDF_CoxMunk()
{
	m_wind_speed = std::numeric_limits<double>::quiet_NaN();
	m_index_of_refraction = std::numeric_limits<double>::quiet_NaN();
}


/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Cox_Munk::SetBRDFParameters		 2017-07-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_BRDF_CoxMunk::SetBRDFParameters(double wind_speed, double index_of_refraction)
{
	m_wind_speed = wind_speed;
	m_index_of_refraction = index_of_refraction;
	return true;
}


/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Cox_Munk::BRDF		 2017-07-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_BRDF_CoxMunk::BRDF(double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* brdf) const
{
	bool ok = NXFINITE(m_index_of_refraction) && NXFINITE(m_wind_speed);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_BRDF_Kernel_CoxMunk::BRDF, one or more of the 2 BRDF parameters is NaN. Cannot calculate BRDF until you define value values for 2 paraemeters");
		*brdf = std::numeric_limits<double>::quiet_NaN();
	}
	else
	{
		// force appropriate cosine values
		CheckCosines(MU_in, MU_out, COSDPHI, "SKTRAN_BRDF_CoxMunk::BRDF");

		// basic sines and cosines
		double cs = MU_in;
		double cv = MU_out;
		double ss = sqrt(1 - cs * cs);
		double sv = sqrt(1 - cv * cv);
		double cp = -COSDPHI; // take the negative to match Sasktran's azimuth convention

		// cosine of the angle of incidence/reflection
		double ceta = cs * cv + ss * sv * cp;
		double eta = nxmath::acosd(ceta);
		double lambda = nxmath::cosd(eta / 2.0);

		// Fresnel reflection terms
		double m2 = m_index_of_refraction*m_index_of_refraction;
		double c = sqrt(m2 + lambda * lambda - 1);
		double r_plus = (m2 * lambda - c) / (m2 * lambda + c);
		double r_minus = (lambda - c) / (lambda + c);
		double R = (r_plus * r_plus + r_minus * r_minus) / 2.0; // Spurr (A.19)

		// water facet orientation probability term
		double sigma2 = 0.003 + 0.00512 * m_wind_speed; // wind speed in m/s
		double gamma = (cs + cv) / (2 * lambda);
		double tau = tan(nxmath::PiOver2 - asin(gamma));
		double P = exp(-tau * tau / (2 * sigma2)) / (4 * cs * cv * nxmath::Pi * sigma2 * pow(gamma, 4)); // Spurr (A.20)

		*brdf = P * R / nxmath::Pi; // Spurr (A.18) - scale down by pi to match Sasktran
	}
	return ok;
}

