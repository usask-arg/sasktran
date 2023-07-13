#include <skopticalproperties21.h>
#include <sources/sasktranif_opticalimpl/skbrdf_stubs.h>

/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Li::ISKBrdf_Stub_Li		 2017-03-06      */
/** **/
/*---------------------------------------------------------------------------*/

ISKBrdf_Stub_Li_Kernel::ISKBrdf_Stub_Li_Kernel(SKTRAN_BRDF_Li_Kernel* librdf)
	: ISKBrdf_Stub_Base(librdf)
{
	m_librdf = librdf;
}


/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Li::~ISKBrdf_Stub_Li		 2017-03-06      */
/** **/
/*---------------------------------------------------------------------------*/

ISKBrdf_Stub_Li_Kernel::~ISKBrdf_Stub_Li_Kernel()
{}

/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Li::SetPropertyArray		 2017-03-06      */
/** **/
/*---------------------------------------------------------------------------*/

bool ISKBrdf_Stub_Li_Kernel::SetPropertyArray(const char* propertyname, const double* value, int numpoints)
{
	nxString					name(propertyname);
	bool						ok = m_librdf != nullptr;

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "ISKBrdf_Stub_Li::SetProperty, the internal C++ object is undefined. Thats not good");
	}
	else
	{
		if (name == "BRDFParameters")
		{
			ok = (numpoints == 2);
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "ISKBrdf_Stub_Li::SetProperty(BRDFParameters), accepts only 2 parameters. You passed in %d parameters", (int)numpoints);
			}
			else
			{
				double br = value[0];
				double hb = value[1];

				ok = m_librdf->SetBRDFParameters(br, hb);
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING, "ISKBrdf_Stub_Snow_Kokhanovsky2012::SetProperty(BRDFParameters), failed to set parameters to b/r=%e, h/b=%e", (double)br, (double)hb);
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
*					SKTRAN_BRDF_Li::SKTRAN_BRDF_Li()		 2017-03-06*/
/** **/
/*---------------------------------------------------------------------------*/
SKTRAN_BRDF_Li_Kernel::SKTRAN_BRDF_Li_Kernel()
{
	m_br = std::numeric_limits<double>::quiet_NaN();
	m_hb = std::numeric_limits<double>::quiet_NaN();
}

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Li::SetBRDFParameters()		 2017-03-06*/
/** **/
/*---------------------------------------------------------------------------*/
bool SKTRAN_BRDF_Li_Kernel::SetBRDFParameters(double br, double hb)
{
	m_br = br;
	m_hb = hb;
	return true;
}

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Li::overlap()		 2017-03-06*/
/** **/
/*---------------------------------------------------------------------------*/
double SKTRAN_BRDF_Li_Kernel::overlap(double s, double v, double cp) const
{
	double ss = nxmath::sind(s);
	double sv = nxmath::sind(v);
	double cs = nxmath::cosd(s);
	double cv = nxmath::cosd(v);
	double ts = ss / cs;
	double tv = sv / cv;
	double sp2 = 1.0 - cp*cp;

	// equation 35 and 50
	double D2 = ts*ts + tv*tv - 2 * ts*tv*cp;
	// equation 34 and 49
	double ct = std::max(-1.0, std::min(1.0, m_hb*sqrt(D2 + ts*ts*tv*tv*sp2)*(cs*cv) / (cs + cv)));
	double t = nxmath::DegreesToRadians(nxmath::acosd(ct));
	double st = sqrt(std::max(0.0, 1 - ct*ct));

	double overlap = (1 / nxmath::Pi)*(t - st*ct)*(cs + cv) / (cs*cv);
	return overlap;
}

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Li::cos_primed_scattering_angle()		 2017-03-06*/
/** **/
/*---------------------------------------------------------------------------*/
double SKTRAN_BRDF_Li_Kernel::cos_primed_scattering_angle(double s, double v, double cp) const
{
	double ss = nxmath::sind(s);
	double sv = nxmath::sind(v);
	double cs = nxmath::cosd(s);
	double cv = nxmath::cosd(v);

	double ceta = cs*cv + ss*sv*cp;
	return ceta;
}


/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Li::primed_angle()		 2017-03-06*/
/** **/
/*---------------------------------------------------------------------------*/
double SKTRAN_BRDF_Li_Kernel::primed_angle(double mu) const
{
	double theta = nxmath::acosd(std::max(-1.0, std::min(1.0, mu)));
	double thetap = nxmath::atan2d(m_br*nxmath::sind(theta), nxmath::cosd(theta));
	return thetap;
}

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_LiSparse::BRDF()		 2017-03-06*/
/** **/
/*---------------------------------------------------------------------------*/
bool SKTRAN_BRDF_LiSparse_Kernel::BRDF(double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* brdf) const
{
	bool ok = NXFINITE(m_br) && NXFINITE(m_hb);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_BRDF_LiSparse::BRDF, one or more of the 2 BRDF parameters is NaN. Cannot calculate BRDF until you define value values for 2 paraemeters");
		*brdf = std::numeric_limits<double>::quiet_NaN();
	}
	else
	{
		CheckCosines(MU_in, MU_out, COSDPHI, "SKTRAN_BRDF_LiSparse_Kernel::BRDF");

		double secs = 1.0 / MU_in;
		double secv = 1.0 / MU_out;
		double cp = -COSDPHI; // take the negative to match Sasktran's azimuth convention

		double s = SKTRAN_BRDF_Li_Kernel::primed_angle(MU_in);
		double v = SKTRAN_BRDF_Li_Kernel::primed_angle(MU_out);

		double O = SKTRAN_BRDF_Li_Kernel::overlap(s, v, cp);
		double ceta = SKTRAN_BRDF_Li_Kernel::cos_primed_scattering_angle(s, v, cp);

		*brdf = (O - secs - secv + 0.5*secv*(1 + ceta)) / nxmath::Pi; // scale down by pi to match sasktran
	}
	return ok;
}

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_LiDense::BRDF()		 2017-03-06*/
/** **/
/*---------------------------------------------------------------------------*/
bool SKTRAN_BRDF_LiDense_Kernel::BRDF(double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* brdf) const
{
	bool ok = NXFINITE(m_br) && NXFINITE(m_hb);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_BRDF_LiDense::BRDF, one or more of the 2 BRDF parameters is NaN. Cannot calculate BRDF until you define value values for 2 paraemeters");
		*brdf = std::numeric_limits<double>::quiet_NaN();
	}
	else
	{
		CheckCosines(MU_in, MU_out, COSDPHI, "SKTRAN_BRDF_LiDense_Kernel::BRDF");

		double secs = 1.0 / MU_in;
		double secv = 1.0 / MU_out;
		double cp = -COSDPHI; // take the negative to match Sasktran's azimuth convention

		double s = SKTRAN_BRDF_Li_Kernel::primed_angle(MU_in);
		double v = SKTRAN_BRDF_Li_Kernel::primed_angle(MU_out);

		double O = SKTRAN_BRDF_Li_Kernel::overlap(s, v, cp);
		double ceta = SKTRAN_BRDF_Li_Kernel::cos_primed_scattering_angle(s, v, cp);

		*brdf = (((1 + ceta) * secv) / (secv + secs - O) - 2) / nxmath::Pi; // scale down by pi to match sasktran
	}
	return ok;
}

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_LiSparseReciprocal::BRDF()		 2017-07-31*/
/** **/
/*---------------------------------------------------------------------------*/
bool SKTRAN_BRDF_LiSparseReciprocal_Kernel::BRDF(double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* brdf) const
{
	bool ok = NXFINITE(m_br) && NXFINITE(m_hb);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_BRDF_LiDense::BRDF, one or more of the 2 BRDF parameters is NaN. Cannot calculate BRDF until you define value values for 2 paraemeters");
		*brdf = std::numeric_limits<double>::quiet_NaN();
	}
	else
	{
		CheckCosines(MU_in, MU_out, COSDPHI, "SKTRAN_BRDF_LiSparseReciprocal_Kernel::BRDF");

		double secs = 1.0 / MU_in;
		double secv = 1.0 / MU_out;
		double cp = -COSDPHI; // take the negative to match Sasktran's azimuth convention

		double s = SKTRAN_BRDF_Li_Kernel::primed_angle(MU_in);
		double v = SKTRAN_BRDF_Li_Kernel::primed_angle(MU_out);

		double O = SKTRAN_BRDF_Li_Kernel::overlap(s, v, cp);
		double ceta = SKTRAN_BRDF_Li_Kernel::cos_primed_scattering_angle(s, v, cp);

		*brdf = (O - secs - secv + 0.5*secs*secv*(1 + ceta)) / nxmath::Pi; // multiplying the last term by secs changes this from LiSparse to LiSparseReciprocal
																		   // scale down by pi to match sasktran
	}
	return ok;
}