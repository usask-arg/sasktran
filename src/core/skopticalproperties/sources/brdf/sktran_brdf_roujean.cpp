#include <skopticalproperties21.h>
#include <sources/sasktranif_opticalimpl/skbrdf_stubs.h>
#include <map>
#include <mutex>


/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Roujean::ISKBrdf_Stub_Roujean		 2016- 12- 19*/
/** **/
/*---------------------------------------------------------------------------*/

ISKBrdf_Stub_Roujean::ISKBrdf_Stub_Roujean(SKTRAN_BRDF_Roujean* roujean)
	: ISKBrdf_Stub_Base(roujean)
{
	m_roujeanbrdf = roujean;
}


/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Roujean::~ISKBrdf_Stub_Roujean		 2016- 12- 19*/
/** **/
/*---------------------------------------------------------------------------*/

ISKBrdf_Stub_Roujean::~ISKBrdf_Stub_Roujean()
{
}



/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Roujean::SetPropertyArray		 2016- 12- 19*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKBrdf_Stub_Roujean::SetPropertyArray(const char* propertyname, const double* values, int npts)
{
	nxString					name(propertyname);
	bool						ok = m_roujeanbrdf != nullptr;

	if (name == "BRDFParameters")
	{
		ok = ok && (npts == 3);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "ISKBrdf_Stub_Roujean::SetProperty(BRDFParameters), must have exactly 3 elements in the array, You passed in %d elements.", (int)npts);
		}
		else
		{
			double k0 = values[0];
			double k1 = values[1];
			double k2 = values[2];
			ok = m_roujeanbrdf->SetBRDFParameters(k0, k1, k2);
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "ISKBrdf_Stub_Roujean::SetProperty(BRDFParameters), there were problems setting BRDF with parameters, k0=%e, k1 = %e, k2 = %e", (double)k0, (double)k1, (double)k2);

			}
		}
	}
	else
	{
		ok = ISKBrdf_Stub_Base::SetPropertyArray(propertyname, values, npts);
	}
	return ok;

}


/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Roujean::SetPropertyScalar		 2016- 12- 20*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKBrdf_Stub_Roujean::SetPropertyScalar(const char* propertyname, double value)
{
	nxString					name(propertyname);
	bool						ok = m_roujeanbrdf != nullptr;

	if (ok)
	{
		if (name == "SetPredefinedParameters")
		{
			int			v = (int)value;
			nxString	pname;

			switch (v)
			{
			case  1: pname = "PLOWED FIELD VIS";     break;
			case  2: pname = "ANNUAL GRASS VIS";     break;
			case  3: pname = "HARD WHEAT VIS";       break;
			case  4: pname = "STEPPE VIS";           break;
			case  5: pname = "CORN VIS";             break;
			case  6: pname = "ORCHARD GRASS VIS";    break;
			case  7: pname = "IRRIGATED WHEAT VIS";  break;
			case  8: pname = "PINEFOREST VIS";       break;
			case  9: pname = "DECIDUOUS FOREST VIS"; break;
			case 10: pname = "SOYBEAN VIS";		  break;
			case 11: pname = "GRASS LAWN VIS";		  break;
			case 21: pname = "PLOWED FIELD NIR";	  break;
			case 22: pname = "ANNUAL GRASS NIR";	  break;
			case 23: pname = "HARD WHEAT NIR";       break;
			case 24: pname = "STEPPE NIR";           break;
			case 25: pname = "CORN NIR";             break;
			case 26: pname = "ORCHARD GRASS NIR";	  break;
			case 27: pname = "IRRIGATED WHEAT NIR";  break;
			case 28: pname = "PINEFOREST NIR";       break;
			case 29: pname = "DECIDUOUS FOREST NIR"; break;
			case 30: pname = "SOYBEAN NIR";          break;
			case 31: pname = "GRASS LAWN NIR";	      break;
			default: ok = false;
				nxLog::Record(NXLOG_WARNING, "ISKBrdf(Roujean) Property SetPredefinedParameters does not recognise value (%d) as a valid parameter set", (int)v);
				break;
			};
			if (ok)
			{
				ok = m_roujeanbrdf->LoadPredefinedParameters(pname);
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING, "ISKBrdf(Roujean) Property SetPredefinedParameters. There were errors selecting parameter (%d)= (%s)", (int)v, (const char*)pname);
				}
			}
		}
		else
		{
			ok = ISKBrdf_Stub_Base::SetPropertyScalar(propertyname, value);
		}

	}
	return ok;
}



/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Roujean::SKTRAN_BRDF_Roujean		 2016- 12- 8*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_BRDF_Roujean::SKTRAN_BRDF_Roujean()
{
	m_k0 = std::numeric_limits<double>::quiet_NaN();
	m_k1 = std::numeric_limits<double>::quiet_NaN();
	m_k2 = std::numeric_limits<double>::quiet_NaN();

}

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Roujean::LoadPredefinedParameters		 2016- 12- 8*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_BRDF_Roujean::LoadPredefinedParameters(const char* paramname)
{
	bool	ok;

	ok = CheckPredefinedSurfaces();
	if (ok)
	{
		nxString	str(paramname);
		str.MakeUpper();
		auto iter = m_predefinedsurfaces.find(std::string((const char*)str));
		ok = !(iter == m_predefinedsurfaces.end());
		if (ok)
		{
			m_k0 = (*iter).second.k0;
			m_k1 = (*iter).second.k1;
			m_k2 = (*iter).second.k2;
		}
		else
		{
			nxLog::Record(NXLOG_WARNING, "KTRAN_BRDF_Roujean::LoadPredefinedParameters, Could not find a pre-defined entry for <%s>", (const char*)str);
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Roujean::SetBRDFParameters		 2016- 12- 8*/
/** Sets the three parameters of the Roujean BRDF model to the user defined
*	values;
**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_BRDF_Roujean::SetBRDFParameters(double k0, double k1, double k2)
{
	bool	ok;

	m_k0 = k0;
	m_k1 = k1;
	m_k2 = k2;
	ok = NXFINITE(m_k0) && NXFINITE(m_k1) && NXFINITE(m_k2);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_BRDF_Roujean::SetBRDFParameters, you cannot use NaN as a value for any of the three Roujean BRDF paremters. The internals values are now corrupted.");
	}
	return ok;

}

std::map< std::string, SKTRAN_BRDF_Roujean_entry>	SKTRAN_BRDF_Roujean::m_predefinedsurfaces;

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Roujean::CheckPredefinedSurfaces		 2016- 12- 8*/
/** Implements the parameters published in Table 1 of Roujean et al. **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_BRDF_Roujean::CheckPredefinedSurfaces()
{
	static std::mutex					syncmutex;
	std::lock_guard< std::mutex>		lock(syncmutex);

	if (m_predefinedsurfaces.size() == 0)
	{
		m_predefinedsurfaces.insert(value_type("PLOWED FIELD VIS", SKTRAN_BRDF_Roujean_entry(24.3, 7.3, 64.2, 630.0)));
		m_predefinedsurfaces.insert(value_type("ANNUAL GRASS VIS", SKTRAN_BRDF_Roujean_entry(34.9, 4.4, 37.7, 630.0)));
		m_predefinedsurfaces.insert(value_type("HARD WHEAT VIS", SKTRAN_BRDF_Roujean_entry(27.3, 5.2, 26.9, 630.0)));
		m_predefinedsurfaces.insert(value_type("STEPPE VIS", SKTRAN_BRDF_Roujean_entry(26.6, 5.0, 5.9, 630.0)));
		m_predefinedsurfaces.insert(value_type("CORN VIS", SKTRAN_BRDF_Roujean_entry(8.4, 0.6, 0.1, 630.0)));
		m_predefinedsurfaces.insert(value_type("ORCHARD GRASS VIS", SKTRAN_BRDF_Roujean_entry(7.9, 1.2, 9.0, 630.0)));
		m_predefinedsurfaces.insert(value_type("IRRIGATED WHEAT VIS", SKTRAN_BRDF_Roujean_entry(5.2, 0.5, 27.3, 630.0)));
		m_predefinedsurfaces.insert(value_type("PINEFOREST VIS", SKTRAN_BRDF_Roujean_entry(3.7, 0.0, 13.3, 630.0)));
		m_predefinedsurfaces.insert(value_type("DECIDUOUS FOREST VIS", SKTRAN_BRDF_Roujean_entry(3.0, 0.0, 8.7, 630.0)));
		m_predefinedsurfaces.insert(value_type("SOYBEAN VIS", SKTRAN_BRDF_Roujean_entry(3.2, 0.0, 8.4, 630.0)));
		m_predefinedsurfaces.insert(value_type("GRASS LAWN VIS", SKTRAN_BRDF_Roujean_entry(4.8, 0.0, 10.2, 630.0)));
		m_predefinedsurfaces.insert(value_type("PLOWED FIELD NIR", SKTRAN_BRDF_Roujean_entry(28.8, 8.5, 74.8, 915.0)));
		m_predefinedsurfaces.insert(value_type("ANNUAL GRASS NIR", SKTRAN_BRDF_Roujean_entry(45.2, 5.3, 50.3, 915.0)));
		m_predefinedsurfaces.insert(value_type("HARD WHEAT NIR", SKTRAN_BRDF_Roujean_entry(37.3, 3.3, 80.2, 915.0)));
		m_predefinedsurfaces.insert(value_type("STEPPE NIR", SKTRAN_BRDF_Roujean_entry(35.6, 5.6, 21.7, 915.0)));
		m_predefinedsurfaces.insert(value_type("CORN NIR", SKTRAN_BRDF_Roujean_entry(27.2, 0.0, 28.5, 915.0)));
		m_predefinedsurfaces.insert(value_type("ORCHARD GRASS NIR", SKTRAN_BRDF_Roujean_entry(26.5, 1.5, 43.0, 915.0)));
		m_predefinedsurfaces.insert(value_type("IRRIGATED WHEAT NIR", SKTRAN_BRDF_Roujean_entry(42.1, 0.0, 121.0, 915.0)));
		m_predefinedsurfaces.insert(value_type("PINEFOREST NIR", SKTRAN_BRDF_Roujean_entry(28.2, 1.7, 24.3, 915.0)));
		m_predefinedsurfaces.insert(value_type("DECIDUOUS FOREST NIR", SKTRAN_BRDF_Roujean_entry(40.0, 4.0, 29.5, 915.0)));
		m_predefinedsurfaces.insert(value_type("SOYBEAN NIR", SKTRAN_BRDF_Roujean_entry(52.8, 1.0, 46.0, 915.0)));
		m_predefinedsurfaces.insert(value_type("GRASS LAWN NIR", SKTRAN_BRDF_Roujean_entry(36.3, 0.0, 56.4, 915.0)));
	}
	return true;
}


/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Roujean::BRDF					2016- 12- 8*/
/** Calculates the scalar BRDF of the Roujean parametrization using equation 2,
*	equation 8 and equation 10 of the Roujean et al paper.
*
*	\param MU_in
*	Cosine of zenith angle \f$\theta_s\f$ of incoming ray. This is  \f$\cos(\theta_s)\f$
*	in the Roujean paper.
*
*	\param MU_out
*	Cosine of zenith angle \f$\theta_v\f$ of outbound ray. This is  \f$\cos(\theta_v)\f$
*	in the Roujean paper.
*
*	\param COSDPHI
*	Cosine of delta azimuth angle. This is \f$cos(\Delta\phi)\f$ in the Roujean paper.

*	Reference:
*	Roujean, J.-L., M. Leroy, and P.-Y. Deschamps (1992), A bidirectional
*	reflectance model of the Earth's surface for the correction of remote
*	sensing data, J. Geophys. Res., 97(D18), 20455–20468, doi:10.1029/92JD01411.
**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_BRDF_Roujean::BRDF(double /*wavelennm*/, const GEODETIC_INSTANT& /*pt*/, double MU_in, double MU_out, double COSDPHI, double* brdf) const
{
	bool			ok;
	double			xmus;
	double			xmuv;
	double			fi;
	double			cosfi;
	double			sinfi;
	double			sinxmus;
	double			sinxmuv;
	double			coseta;
	double			sineta;
	double			ttv;
	double			tts;
	double			eta;
	double			f2;
	double			f1;

	ok = NXFINITE(m_k0) && NXFINITE(m_k1) && NXFINITE(m_k2);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_BRDF_Roujean::BRDF, one or more of the 3 BRDF parameters is NaN. Cannot calculate BRDF until you define value values for 3 paraemeters");
		*brdf = std::numeric_limits<double>::quiet_NaN();
	}
	else
	{
		// Make sure we are in the range 0.017 to +1 (ie 89 degrees to 0 degrees) for all of the incoming and outgoing cosines
		// Make sure cosine of delta azimuth is between -1 and + 1
		CheckCosines(MU_in, MU_out, COSDPHI, "SKTRAN_BRDF_Roujean::BRDF");
		
		xmus = MU_in;
		xmuv = MU_out;
		cosfi = -COSDPHI;	// take the negative to match Sasktran's azimuth convention
		fi = acos(cosfi);	// Get phi in the range 0 to pi
		sinxmus = sqrt(1.0 - xmus*xmus);
		sinxmuv = sqrt(1.0 - xmuv*xmuv);
		sinfi = sqrt(1.0 - cosfi*cosfi);
		tts = sinxmus / xmus;
		ttv = sinxmuv / xmuv;

		coseta = xmus*xmuv + sinxmus*sinxmuv*cosfi;																// Get the cosine of the scattering angle (presented immediately after equation 6 in the paper)
		coseta = std::max(-1.0, std::min(1.0, coseta));															// Make sure we dont have floating point roundoff issues
		sineta = sqrt(1.0 - coseta*coseta);
		eta = acos(coseta);																						// Get eta in the range 0 to pi
		f2 = 4.0 / (3.0*nxmath::Pi*(xmus + xmuv))*((nxmath::PiOver2 - eta)*coseta + sineta) - (1.0 / 3.0);		// Equation 8 of Roujean et al
		f1 = 0.5 / nxmath::Pi*((nxmath::Pi - fi)*cosfi + sinfi)*tts*ttv - (tts + ttv + sqrt(tts*tts + ttv*ttv - 2 * tts*ttv*cosfi)) / nxmath::Pi;	// Equation 2 of Roujean et al.
		*brdf = 0.01*(m_k0 + m_k1*f1 + m_k2*f2) / nxmath::Pi;													// scale from percentages (0-100) to decimal (0-1)
																												// scale down by pi to match sasktran
	}
	return ok;
}

/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Roujean_Kernel::ISKBrdf_Stub_Roujean		 2017-08-02*/
/** **/
/*---------------------------------------------------------------------------*/

ISKBrdf_Stub_Roujean_Kernel::ISKBrdf_Stub_Roujean_Kernel(SKTRAN_BRDF_Roujean_Kernel* roujean)
	: ISKBrdf_Stub_Base(roujean)
{
	m_roujeankernel = roujean;
}

ISKBrdf_Stub_Roujean_Kernel::~ISKBrdf_Stub_Roujean_Kernel() {}

SKTRAN_BRDF_Roujean_Kernel::SKTRAN_BRDF_Roujean_Kernel() {}

bool SKTRAN_BRDF_Roujean_Kernel::BRDF(double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* brdf) const
{
	double			xmus;
	double			xmuv;
	double			fi;
	double			sinxmus;
	double			sinxmuv;
	double			cosfi;
	double			sinfi;
	double			ttv;
	double			tts;
	double			f1;

	// Make sure we are in the range 0.017 to +1 (ie 89 degrees to 0 degrees) for all of the incoming and outgoing cosines
	// Make sure cosine of delta azimuth is between -1 and + 1
	CheckCosines(MU_in, MU_out, COSDPHI, "SKTRAN_BRDF_Roujean_Kernel::BRDF");

	xmus = MU_in;
	xmuv = MU_out;
	cosfi = -COSDPHI;	// take the negative to match Sasktran's azimuth convention

	fi = acos(cosfi);
	sinxmus = sqrt(1.0 - xmus*xmus);
	sinxmuv = sqrt(1.0 - xmuv*xmuv);
	sinfi = sqrt(1.0 - cosfi*cosfi);
	tts = sinxmus / xmus;
	ttv = sinxmuv / xmuv;

	f1 = 0.5 / nxmath::Pi*((nxmath::Pi - fi)*cosfi + sinfi)*tts*ttv - (tts + ttv + sqrt(tts*tts + ttv*ttv - 2 * tts*ttv*cosfi)) / nxmath::Pi;	// Equation 2 of Roujean et al.
	*brdf = f1 / nxmath::Pi;																													// scale down by pi to match sasktran
	
	return NXFINITE(*brdf);
}



#if 0
/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF::BRDF		 2016- 11- 24*/
/** **/
/*---------------------------------------------------------------------------*/

bool  SKTRAN_BRDF::BRDF(const SKTRAN_SourceTermQueryObject_Base& incoming,
	const SKTRAN_SourceTermQueryObject_Base& outgoing, SKTRAN_Stokes_NC* source) const
{
	double costhetain;
	double costhetaout;
	const HELIODETIC_POINT&      point = incoming.GetPoint();
	const HELIODETIC_UNITVECTOR& upH = point.LocalZenith();
	HELIODETIC_UNITVECTOR		 localdirs[3];								// Gets local direction unit vectors: north, west, up
	const HELIODETIC_UNITVECTOR& indirH = incoming.GetLookAway();
	const HELIODETIC_UNITVECTOR& outdirH = outgoing.GetLookAway();

	nxVector					 up(upH.X(), upH.Y(), upH.Z());
	nxVector					 indir(indirH.X(), indirH.Y(), indirH.Z());
	nxVector					 outdir(outdirH.X(), outdirH.Y(), outdirH.Z());
	nxVector					 north(localdirs[0].X(), localdirs[0].Y(), localdirs[0].Z());	// This not true geographic north but is the north which is towards the sun direction 
	nxVector					 west(localdirs[1].X(), localdirs[1].Y(), localdirs[1].Z());	// NOt true geographic west bu 90 degrees from the sun (ie either dawn or dusk not sure which one)
	;

	costhetain = (up & indir);
	costhetaout = (up & outdir);
	NXASSERT((costhetain >= 0.0));
	NXASSERT((costhetaout >= 0.0));

	point.LocalUnitVectors(localdirs, N_ELEMENTS(localdirs));
	double inx = indir & north;
	double iny = indir & west;
	double inazi = nxmath::atan2d(iny, inx);
	double outx = outdir & north;
	double outy = outdir & west;
	double outazi = nxmath::atan2d(outy, outx);

	return true;
}


/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Lambertian::SKTRAN_BRDF_Lambertian		 2016- 11- 24*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_BRDF_Lambertian::SKTRAN_BRDF_Lambertian()
{
	m_albedo = 0.3;
}

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Lambertian::BRDF		 2016- 11- 24*/
/** **/
/*---------------------------------------------------------------------------*/

bool  SKTRAN_BRDF_Lambertian::BRDF(const SKTRAN_SourceTermQueryObject_Base& incoming,
	const SKTRAN_SourceTermQueryObject_Base& outgoing,
	SKTRAN_Stokes_NC*						 source) const
{
	const HELIODETIC_POINT&      point = incoming.GetPoint();
	const HELIODETIC_UNITVECTOR& up = point.LocalZenith();
	const HELIODETIC_UNITVECTOR& indir = incoming.GetLookAway();

	double costhetain = (up & indir);
	source->SetTo(m_albedo*costhetain / nxmath::Pi);
	return true;
}



/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_ThetaInOutAndDeltaAzi::SKTRAN_BRDF_ThetaInOutAndDeltaAzi		 2016- 11- 24*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_BRDF_ThetaInOutAndDeltaAzi::SKTRAN_BRDF_ThetaInOutAndDeltaAzi()
{
}

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_ThetaInOutAndDeltaAzi::~SKTRAN_BRDF_ThetaInOutAndDeltaAzi		 2016- 11- 24*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_BRDF_ThetaInOutAndDeltaAzi::~SKTRAN_BRDF_ThetaInOutAndDeltaAzi()
{
}

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_ThetaInOutAndDeltaAzi::BRDF		 2016- 11- 24*/
/** **/
/*---------------------------------------------------------------------------*/

bool  SKTRAN_BRDF_ThetaInOutAndDeltaAzi::BRDF(const SKTRAN_SourceTermQueryObject_Base& incoming,
	const SKTRAN_SourceTermQueryObject_Base& outgoing,
	SKTRAN_Stokes_NC* source) const
{
	bool   ok;
	double costhetain;
	double costhetaout;
	const HELIODETIC_POINT&      point = incoming.GetPoint();
	const HELIODETIC_UNITVECTOR& upH = point.LocalZenith();
	HELIODETIC_UNITVECTOR		 localdirs[3];								// Gets local direction unit vectors: north, west, up
	const HELIODETIC_UNITVECTOR& indirH = incoming.GetLookAway();
	const HELIODETIC_UNITVECTOR& outdirH = outgoing.GetLookAway();

	nxVector					 up(upH.X(), upH.Y(), upH.Z());
	nxVector					 indir(indirH.X(), indirH.Y(), indirH.Z());
	nxVector					 outdir(outdirH.X(), outdirH.Y(), outdirH.Z());
	nxVector					 north(localdirs[0].X(), localdirs[0].Y(), localdirs[0].Z());	// This not true geographic north but is the north which is towards the sun direction 
	nxVector					 west(localdirs[1].X(), localdirs[1].Y(), localdirs[1].Z());	// NOt true geographic west bu 90 degrees from the sun (ie either dawn or dusk not sure which one)


	costhetain = (up & indir);
	costhetaout = (up & outdir);
	NXASSERT((costhetain >= 0.0));
	NXASSERT((costhetaout >= 0.0));

	point.LocalUnitVectors(localdirs, N_ELEMENTS(localdirs));
	double inx = indir & north;
	double iny = indir & west;
	double inazi = nxmath::atan2d(iny, inx);
	double outx = outdir & north;
	double outy = outdir & west;
	double outazi = nxmath::atan2d(outy, outx);
	double	deltaazi = nxmath::inrange((outazi - inazi), 360.0);

	ok = Brdf_impl(costhetain, costhetaout, deltaazi, source);
	return ok;
}

#endif

