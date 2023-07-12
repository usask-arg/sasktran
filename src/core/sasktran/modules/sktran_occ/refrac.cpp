#include "skoccultation.h"
//#include "include/sktran_occ_internals.h"
//#include "include/sktran_curvedrays_ace.h"

/*-----------------------------------------------------------------------------
 *					skRTRefractiveIndex_ACEFTSProfile::skRTRefractiveIndex_ACEFTSProfile		2014-1-28*/
/** **/
/*---------------------------------------------------------------------------*/

skRTRefractiveIndex_ACEFTSProfile::skRTRefractiveIndex_ACEFTSProfile()
{
	m_atmosphericstate = NULL;
	m_raytracingshells = NULL;
	m_Q       = std::numeric_limits<double>::quiet_NaN();
}


/*-----------------------------------------------------------------------------
 *					skRTRefractiveIndex_ACEFTSProfile::~skRTRefractiveIndex_ACEFTSProfile		2014-1-28*/
/** **/
/*---------------------------------------------------------------------------*/

skRTRefractiveIndex_ACEFTSProfile::~skRTRefractiveIndex_ACEFTSProfile()
{
	ReleaseResources();
}

/*-----------------------------------------------------------------------------
 *					skRTRefractiveIndex_ACEFTSProfile::ReleaseResources		2014-1-28*/
/** **/
/*---------------------------------------------------------------------------*/

void skRTRefractiveIndex_ACEFTSProfile::ReleaseResources()
{
	if (m_atmosphericstate != NULL) m_atmosphericstate->Release();
//	if (m_raytracingshells != NULL) m_raytracingshells->Release();
	m_atmosphericstate = NULL;
//	m_raytracingshells = NULL;
	m_Q                = std::numeric_limits<double>::quiet_NaN();
}

/*-----------------------------------------------------------------------------
 *					skRTRefractiveIndex_ACEFTSProfile::Initialize		2014-1-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTRefractiveIndex_ACEFTSProfile::Initialize( skClimatology* atmospheric_raytracing_state, std::shared_ptr<const SKTRAN_GridDefRayTracingShells_V21>	raytracingshells  )
{
	if ( atmospheric_raytracing_state != NULL) atmospheric_raytracing_state->AddRef();
	if ( m_atmosphericstate           != NULL) m_atmosphericstate->Release();
	m_atmosphericstate = atmospheric_raytracing_state;
	m_raytracingshells = raytracingshells;
	return true;
}
	
/*-----------------------------------------------------------------------------
 *					skRTRefractiveIndex_ACEFTSProfile::INTERP		2014-1-27*/
/** INTERPOLATE LOG PRESSURE ON A HEIGHT SCALE and return Refractive index (R)
 *	and pressure (P), height Z.  
 **/
/*---------------------------------------------------------------------------*/

bool skRTRefractiveIndex_ACEFTSProfile::INTERP( double h_meters, double* userR ) const
{
	bool	ok;
	double	P;
	double	R;
	double  T;
	double	D;

	ok =  NXFINITE(m_Q);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skRTRefractiveIndex_ACEFTSProfile::INTERP, The refractive index code is not properly initialized");
		R = std::numeric_limits<double>::quiet_NaN() ;
	}
	else
	{
		P = m_P.Interpolate( h_meters, std::numeric_limits<double>::quiet_NaN() );
		T = m_T.Interpolate( h_meters, std::numeric_limits<double>::quiet_NaN() );
		ok = NXFINITE(P) && !NXFINITE(T);
		if (!ok)
		{
			R = std::numeric_limits<double>::quiet_NaN() ;
		}
		else
		{
			P = exp(P);
			D = 353.0*P/T;														// PERFECT GAS LAW to get the number density (or is it molar density
			R = NXFINITE(D) ? sqrt( (1.0 +2.0*m_Q*D)/( 1.0 - m_Q*D) ) : 1.0;	// LORENZIAN-LORENZIAN LAW (n**2-1)/(n**2+2) = 4.pi/3 N .alpha,   m_Q = alpha, D = 4.pi/3 N
			NXTRACE_ONCEONLY(firsttime, ("skRTRefractiveIndex_ACEFTSProfile::INTERP, calculating refractivity rather than refractive index would allow a more accurate interpolation in height on a log scale\n"));
		}
	}
	*userR = R;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skRTRefractiveIndex_ACEFTSProfile::SetProfileLocation		2014-1-31*/
/** Sets the location of the Refractive index profile. This forces the
 *	refractive index object to recalculate its internal height profile
 **/
/*---------------------------------------------------------------------------*/

bool skRTRefractiveIndex_ACEFTSProfile::SetProfileLocationAndWavenumber( const GEODETIC_INSTANT& apoint, double NUS)
{
	std::vector<double>	PA;
	std::vector<double>	TA;
	std::vector<double>	HA;
	double				P;
	double				T;
	GEODETIC_INSTANT	point = apoint;
	size_t				i;
	bool				ok1;
	bool				ok;


	ok = ( m_raytracingshells != NULL ) && (m_atmosphericstate != NULL);
	ok = ok && SetWavenumber(NUS);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skRTRefractiveIndex_ACEFTSProfile::SetProfileLocation, the class is not properly initialized. We cannot do the spline interpolation./ Thats not good");
	}
	else
	{
		PA.reserve( m_raytracingshells->NumShells() );
		TA.reserve( m_raytracingshells->NumShells() );
		HA.reserve ( m_raytracingshells->NumShells() );

#if defined(REPLICATE_ACE)
		for (i = 0; i < m_raytracingshells->NumCells(); i++)
		{
			point.heightm = 0.5*(m_raytracingshells->At(i) + m_raytracingshells->At(i+1));						// To replicate ACE use the middle of the shells, keeps consistency with ACE code to within a few microns
#else
		for (i = 0; i < m_raytracingshells->NumShells(); i++)													// But if not replicating ACE its better
		{																										// To use the shell boundary. But consistency with ACE total path length is about 0.5 meters 
			point.heightm = m_raytracingshells->At(i);
#endif
			ok1 =         m_atmosphericstate->GetParameter( SKCLIMATOLOGY_PRESSURE_PA,   point,  &P,  false);	// Get the pressure from the user supplied climatology
			ok1 =  ok1 && m_atmosphericstate->GetParameter( SKCLIMATOLOGY_TEMPERATURE_K, point,  &T,  false);	// Get the Temperature from the user supplied climatology
			if (ok1 && NXFINITE(P) && NXFINITE(T))
			{
				P /= 101325.0;																					// Convert Pressure in Pascals to atmospheres
				PA.push_back( log(P) );
				TA.push_back( T );
				HA.push_back(point.heightm);
			}
		}
		ok = ok && m_P.Configure( HA, PA);
		ok = ok && m_T.Configure( HA, TA);
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skRTRefractiveIndex_ACEFTSProfile::SetProfileLocation, There were errors creating the refractive index interpolation. Thats not good");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayRefracted_ThomPepSim::ComputeRefractionOfLayers		2014-1-27*/
/** Copute the refractive index of each layer using pressure and temperature.
 *  Filippenko (1982, PASP, 94, 715). 
**/
/*---------------------------------------------------------------------------*/

bool skRTRefractiveIndex_ACEFTSProfile::SetWavenumber(	double NUS )
{
	double	T = NUS*1.0E-4;							// 1/LAMBDA (MICRONS)
	double	A, B, C, AA;
    
#if defined(REPLICATE_ACE)
	nxLog::Record(NXLOG_WARNING,"skRTRefractiveIndex_ACEFTSProfile::ComputeRefractionOfLayers,The refractive index calculation is only for nu = 2400, not the requested wavenumber %12.5f", (double)NUS);
	T    = 2400.0*1.0E-4 ;										// 1/LAMBDA (MICRONS)  
#endif
	T    = T*T;
	A    = 64.328;												// CONSTANTS USED IN
	B    = 29498.1/(146.0 - T);									// CALCULATION OF
	C    = 255.4/(41.0 - T);									// REFRACTIVE INDEX. I'm guessing this is some form of Edlens dry air refractive index formula.
	AA   = 1.0E-6*(A+B+C);
	m_Q  = AA*( AA + 2.0 )/( 1.225014 * (AA*(AA + 2.0) + 3.0));
	return true;
}


/*-----------------------------------------------------------------------------
 *					skRTRefractiveIndex_ACEFTSProfile::InterpolateRefractiveIndex		2014-1-31*/
/** **/
/*---------------------------------------------------------------------------*/

double skRTRefractiveIndex_ACEFTSProfile::RefractiveIndex( double h_meters ) const
{
	double r;

	if ( h_meters < m_raytracingshells->front() ) h_meters =  m_raytracingshells->front();
	if ( h_meters > m_raytracingshells->back()  ) h_meters =  m_raytracingshells->back();
	INTERP( h_meters, &r );
	if (!NXFINITE(r))
	{
		r = 1.0;
	}
	return r;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_RayRefracted_TrajectoryData::ReserveSpace		2014-2-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayRefracted_TrajectoryData::ReserveSpace( size_t n )
{
	m_trajectorypoints.resize(0);
	m_trajectorypoints.reserve(n);
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayRefracted_TrajectoryData::PushBack		2014-2-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayRefracted_TrajectoryData::PushBack( double r, double angle, double pathlength)
{
	pointentry	entry;
	entry.r = r;
	entry.angle = angle;
	entry.length = pathlength;
	m_trajectorypoints.push_back( entry );
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayRefracted_TrajectoryData::MirrorPoints		2014-2-4*/
/** Mirror the points on to the set of data points. This is useful when
 *	tracing a ray through the atmosphere. Exclude the very last point**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayRefracted_TrajectoryData::MirrorPoints()
{
	iterator	iter;
	iterator	nextiter;
	iterator	firstpoint;
	pointentry	entry;

	iter       = m_trajectorypoints.end();
	firstpoint = m_trajectorypoints.begin();
	++firstpoint;
	nextiter   = iter;
	--nextiter;												// Skip one so we dont duplicate the tangent point

	while (!(nextiter == firstpoint))
	{
		--iter;
		--nextiter;
		entry.r     = (*nextiter).r;
		entry.angle = (*iter).angle;
		entry.length= (*iter).length;
		m_trajectorypoints.push_back(entry);
	}
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayRefracted_TrajectoryData::MirrorPoints		2014-2-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayRefracted_TrajectoryData::SumAnglesAndLengths()
{
	iterator	iter;
	pointentry	entry;

	iter = m_trajectorypoints.begin();
	entry = *iter;
	while (!(iter == m_trajectorypoints.end()))
	{
		entry.r       = (*iter).r;
		entry.angle  += (*iter).angle;
		entry.length += (*iter).length;
		(*iter) = entry;
		++iter;
	}
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayRefracted_TrajectoryData::ConvertTo3DLocation		2014-2-4*/
/** **/
/*---------------------------------------------------------------------------*/

HELIODETIC_VECTOR SKTRAN_RayRefracted_TrajectoryData::ConvertTo3DLocation( const HELIODETIC_VECTOR& observer, const HELIODETIC_POINT& observerpt, const HELIODETIC_UNITVECTOR& look, const_iterator& iter )
{
	double				r1;
	double				r2;
	double				phi;
	double				l;
	double				sintheta;
	double				costheta;
	nxVector			lookv( look.X(),     look.Y(),     look.Z() );
	nxVector			obs  ( observer.X(), observer.Y(), observer.Z() );
	nxVector			lookc;
	nxVector			horiz;
	nxVector			vert;
	HELIODETIC_VECTOR	trajectorypoint;

	r1  = observerpt.Radius();						// Use the cosine rule to get the straight line distance from 
	r2  = (*iter).r;								// the observer to the trajectory point
	phi = (*iter).angle;
	l = r1*r1 + r2*r2 - 2.0*r1*r2*cos(phi);
	if (l < 0.0)
	{
		if ( l < -0.00000001) nxLog::Record(NXLOG_WARNING, "SKTRAN_RayRefracted_TrajectoryData::ConvertTo3DLocation, distance from observer  to point does not obey coisne rule. Wow thats weird.");
		l = 0.0;
	}
	l = sqrt(l);															// Get distance to the trajectory point.

	sintheta = r2*sin(phi)/l;												// Use the sine rule to get the straight line nadir angle at the observer
	costheta = sqrt( 1.0 - sintheta*sintheta);								// Pythagorus gives us the cosine
	horiz    = lookv.ComponentPerpendicularTo( obs ).UnitVector();			// Get the horizontal looking component in the plane of the look vector
	vert     = obs.UnitVector();											// Get the vertical looking unitvector
	lookc    = (horiz*sintheta) - (vert*costheta);							// Get the look vector from the observer to the curved ray point
	obs     +=  (lookc*l);														// Go from the observer to the trajectory point
#if defined(NXDEBUG)
	double deflection;
	double rnew;
	deflection = lookc.AngleTo( lookv);
	rnew = obs.Magnitude();
	NXASSERT( deflection < 10 );												// we should never see deflections of more than 0.5 degrees or so
	NXASSERT( fabs( rnew - r2 ) < 0.001 );										// our final radius should be the same as the original 
#endif
	trajectorypoint.SetCoords( obs.X(), obs.Y(), obs.Z() );
	return trajectorypoint;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayRefracted_ThomPepSim::SKTRAN_RayRefracted_ThomPepSim		2014-1-28*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_RayRefracted_ThomPepSim::SKTRAN_RayRefracted_ThomPepSim()
{
//	m_raytracingshells = NULL;
	m_zerokmradius     = std::numeric_limits<double>::quiet_NaN();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayRefracted_ThomPepSim::~SKTRAN_RayRefracted_ThomPepSim		2014-1-28*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_RayRefracted_ThomPepSim::~SKTRAN_RayRefracted_ThomPepSim()
{
//	if (m_raytracingshells != NULL) m_raytracingshells->Release();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayRefracted_ThomPepSim::Configure		2014-1-31*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayRefracted_ThomPepSim::Configure( std::shared_ptr< const SKTRAN_GridDefRayTracingShells_V21> raytracingshells, std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords)
{

//	if ( raytracingshells   != NULL) raytracingshells->AddRef();
//	if ( m_raytracingshells != NULL) m_raytracingshells->Release();
	m_raytracingshells = raytracingshells;
	m_coordinates      = coords;

	if (m_coordinates.get() != nullptr)  m_zerokmradius = m_coordinates->AltitudeToRadius( 0.0 );
	else                                 m_zerokmradius = std::numeric_limits<double>::quiet_NaN();

	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayRefracted_ThomPepSim::SetRILocationAndWavenumber		2014-2-1*/
/** Sets up the location and wavenumber used for the refractive index 
 *	profile.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayRefracted_ThomPepSim::SetRILocationAndWavenumber( const GEODETIC_INSTANT& point, double NUS )
{
	return m_RI.SetProfileLocationAndWavenumber( point, NUS);
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayRefracted_ThomPepSim::EarthRadius		2014-1-29*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_RayRefracted_ThomPepSim::EarthRadius( double TANLAT )
{
	double			T1;
	const double	eqrad = 6378137.0;
	double			earthrd;

	T1 = nxmath::DegreesToRadians( TANLAT);
	//    earthrd = eqrad / SQRT ( COS(t1)**2 + 1.0067395D0 * SIN(t1)**2 )
	earthrd = eqrad/ sqrt(1.0 + 0.006739497*nxmath::sqr(sin(T1)));
	return earthrd;
}

/*-----------------------------------------------------------------------------
 *					ComputeApproximateTangentHeight		2014-1-23*/
/** Use Bougers Law at the observer to approximate the radius of the
 *	tangent point by assuming the tangent point occurs when NADANG = 90 (thats always true) and the refractive index at the tangent point is 1.0 (not so true)
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayRefracted_ThomPepSim::ComputeApproximateTangentRadius( double HOBS, double RSPEC, double NADANG, double * RT) const
{
	double	RDXOB  =1.0;												// (RT).  ASSUME FOR THE FIRST ITERATION THAT REFTAN=1.000

	RDXOB = m_RI.RefractiveIndex( HOBS);					// IF THE OBSERVER IS INSIDE THE ATMOS. FIND REFRAC. INDEX AT THAT HEIGHT, Returns RDXOB and P
	*RT   =  RSPEC*sin(NADANG)*RDXOB;									// Use Bougers law, RSPEC = Observer radius. RT is a conserved quatity and close to the radius io fthe tangent point. SET REFTAN*RT=R(SPACECRAFT)*SIN(NADANG)*REFRACTIVE INDEX AT OBSERVER
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayRefracted_ThomPepSim::SetAtmosphericState		2014-1-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayRefracted_ThomPepSim::SetRIAtmosphericState( skClimatology* atmosphericstate )
{
	bool	ok;

	ok = (m_raytracingshells != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_RayRefracted_ThomPepSim::SetAtmosphericState, m_raytracingshells not defined. Please call Configure before calling SetAtmopshericState");
	}
	else
	{
		ok    = m_RI.Initialize( atmosphericstate, m_raytracingshells );
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayRefracted_ThomPepSim::FindGlobalTangentPoint		2014-1-28*/
/** Finds the global tangent point and the minimum radius. The tangent point
 *	can be below the groundand the code relies upon the refractive index code to
 *	truncate refractive index profiles at the top and bottom of the atmosphere
 *	so this code works properly. This code usually converges very quickly for
 *	a sensible atmosphere.
 *
 *	\param RT
 *		On input contains the approximate radius of the tangent point using straight line geometry (refractive index = 1.0)
 *		from a call to #ComputeApproximateTangentRadius. Upon exit contains the new estimate of the tangent point radius
 *		in meters. This height can be well below the ground. In principle it can also be above the atmosphere
 *		but this option is usually filtered out by earlier code.
 *
 *	\param minradius
 *		Returns the minimum radius of a ray in the atmosphere. This will either be identical to
 *		the tangent radius of it will equal the ground radius.
 *
 *	\param groundishit
 *		Returns true if the groundishit.
 *	
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayRefracted_ThomPepSim::FindGlobalTangentPoint(  double *RT, double* minradius, bool* groundishit) const
{

	double				Z;
	double				REFTAN;
	double				RT1;
	double				RT2;
	double				DIFF;
	double				ERR;
	bool				converged = false;
	size_t				icount;
	bool				ok = true;
	double				groundheight;


	ERR          = 1.0E-05;												// Converge to within 10 microns				
	RT2          = *RT;
	Z            = RadiusToAltitude( *RT );					// KNOW APPROX. RADIUS OF TANGENT HEIGHT (RT) - GET PRESSURE AND TEMPERATURE AND GET NEW REFTAN (RT) FOR NEXT ITERATION
	groundheight = m_raytracingshells->front();							// Get the radius of the ground
	icount = 0;
	while ( !converged && (icount < 50))
	{
		REFTAN = m_RI.RefractiveIndex( Z );
		RT1    = (*RT)/REFTAN;										// CALCULATE NEW TANGENT HEIGHT using Bougers Law applied to tangent point (ie sin(angle) = 1.0)
		DIFF   = fabs(RT2-RT1);
		if (icount >= 50)
		{
				RT1 = 0.5*(RT1+RT2);
		}
		else
		{
			RT2  = RT1;
			Z    = RadiusToAltitude(RT2);
			icount++;
		}
		converged = (DIFF <= ERR);
	}
	*RT = RT1;

	if (Z <= groundheight)
	{
		*minradius   = AltitudeToRadius(groundheight);
		*groundishit = true;
	}
	else
	{
		*minradius   = RT1;
		*groundishit = false;
	}
	return converged;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayRefracted_ThomPepSim::IntegrateCurvedPathLength		2014-1-28*/
/**	Integrate along along a curved ray from one radius to another.
 *	For a description of the integration method used see:
 *	"Ray tracing in a refracting spherically symmetric atmosphere": THOMPSON, PEPIN, & SIMON,
 *	J.OPT.SOC.AM., 72, 1498-1501, 1982.
 **/
/*---------------------------------------------------------------------------*/

double SKTRAN_RayRefracted_ThomPepSim::IntegrateCurvedPathLengthOld(  double RT, double REFTAN,  double RJ, double RJ1, int RINC) const
{
	double		START;
	double		FIN;
	double		DELX;
	double		PL;
	int			I1;
	double		X;
	double		Z1;
	double		RIDEX;
	double		A1;
	double		B1;
	double		B2;
	double		B3;
	double		C1;
	double		F;
	double		G;
	double		DELTA;



	START = sqrt(RJ - RT);
	FIN   = sqrt(RJ1- RT);
	DELX  = fabs(START-FIN)/((double)(RINC));
	PL    = 0.0;

	for (I1=0; I1 < RINC; I1++)									// Integrate path length across this shell from one radius to the next the summation from START to FIN 
	{															// so integrate																		
		X  = START + DELX*(I1+0.5);								// from the lower radius to the upper radius
		Z1 = X*X + RadiusToAltitude(RT);			// CALCULATE the height of this radius
		RIDEX = m_RI.RefractiveIndex( Z1 );							// calculate the refractive index
      
		A1 = sqrt( X*X + 2.0*RT);
		B2 = (RIDEX-REFTAN)/RIDEX*RT/(X*X);
		B3 = 1.0 + B2;
		B1 = sqrt( fabs(B3) );							// IF B3 < 0 THEN HAVE AN IMAGINARY SITUATION
		if (B3 < -1.0E-06)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_RayRefracted_ThomPepSim::IntegrateCurvedPathLength, the B3 coefficient is badly negative. Thats a problem. Dont trust answers");
		}
		C1 = sqrt(X*X + (RIDEX + REFTAN)/RIDEX*RT);
		F  = A1+B1*C1;																		// Equation 22 of Thompson, Pepin and Simon
		G  = A1*B1*C1;																		// Equation 23 of Thompson, Pepin and Simon
		PL = PL + ((RIDEX+REFTAN)/RIDEX)*((REFTAN -RIDEX)/RIDEX)*(X*X + RT)/(X*X*F*G);		// Equation 21 of Thompson, Pepin and Simon
	}
	DELTA = sqrt(RJ1*RJ1 - RT*RT) - sqrt(RJ*RJ - RT*RT);
	PL   *= 2.0*RT*RT*DELX;																	// Multiplier term inside integral of Equation 21 of Thompson, Pepin and Simon
	PL   += DELTA;																			// Constant term at end of Equation 21 of Thompson, Pepin and Simon
	return PL;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayRefracted_ThomPepSim::FindGlobalTangentPoint		2014-1-28*/
/** Returns the radial height of the curved ray, whose curved ray constant
 *	is "k", at the location where sin(theta) is sintheta
 **/
/*---------------------------------------------------------------------------*/

double SKTRAN_RayRefracted_ThomPepSim::RadialHeightOfCosTheta( double k, double sintheta ) const
{

	double				z;
	double				r;
	double				n;
	double				RT1;
	double				RT2;
	double				DIFF;
	double				ERR = 1.0E-05;											// Converge to within 10 microns;
	bool				converged = false;
	size_t				icount;

	r    = k/sintheta;																// Approximate radius given Refractive index is 1. This will systematically overshoot by a little bit.
	RT2   = r;
	z     = RadiusToAltitude( r );												// KNOW APPROX. RADIUS OF TANGENT HEIGHT (RT) - GET PRESSURE AND TEMPERATURE AND GET NEW REFTAN (RT) FOR NEXT ITERATION
	icount = 0;
	while ( !converged && (icount < 50))
	{
		n    = m_RI.RefractiveIndex( z );
		RT1  = k/(n*sintheta);
		DIFF = fabs(RT2-RT1);
		if (icount >= 50)
		{
			RT1 = 0.5*(RT1+RT2);
		}
		else
		{
			RT2  = RT1;
			z    = RadiusToAltitude(RT2);
			icount++;
		}
		converged = (DIFF <= ERR);
	}
	r = RT1;
	return r;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayRefracted_ThomPepSim::IntegrateCurvedPathLength_Technique2		2014-2-3*/
/** Integrates the angular distance and teh curved path length from a lower
 *	radius to an upper radius. This integration follows  the algorithm
 *	outlines in word file "\docs\Ray tracing in a spherically symmetric atmosphere.docx".
 *	The technique seems to give answers within 1 or 2 mm of the "old" Thompson technique,
 *	IntegrateCurvedPathLengthOld, even on *	60 km path lengths.
 *
 *	This function has the advantage that it also calculates the angular distance
 *	of the curved segment as well as the path length. This is necessary to work out the
 *	3-D location of the ray.
 *
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayRefracted_ThomPepSim::IntegrateCurvedPathLength(  double rt, double nt,  double r1, double r2,  size_t numsteps, double*L, double *PSI) const
{
	double		k;
	double		n1;
	double		n2;
	double		theta1;
	double		theta2;
	double		deltheta;
	double		thetastart;
	double		thetaend;
	double		theta;
	double		n;
	double		dndr;
	double		costheta;
	double		denom;
	double		dl;
	double		dpsi;
	double		r;
	double		z;
	double		n2costheta;
	
	k         = rt*nt;																			// Get the curved ray constant "k".
	n1        = m_RI.RefractiveIndex ( RadiusToAltitude(r1));									// Get the refractive index of the lower radius
	n2        = m_RI.RefractiveIndex ( RadiusToAltitude(r2));									// Get the refractive index of the upper radius
	theta1    =  acos( k/(n1*r1) );																// get theta for the lower radius
	theta2    =  acos( k/(n2*r2) );																// get theta for the upper radius
	deltheta  =  (theta2 - theta1)/numsteps;													// get the change in theta
	*L        = 0.0;																			// reset the path length integral
	*PSI      = 0.0;																			// reset the angular distance integral

	thetastart = theta1;																		// integrate from theta at the lower radius
	for (size_t i=0; i < numsteps; i++)															// Integrate to the upper radius theta in numsteps
	{																							// so
		thetaend   = thetastart + deltheta;														// Get the end theta angle of this section
		theta      = 0.5*( thetaend + thetastart );												// get theta in the middle of the section
		costheta   = cos(theta);																// get the cosine
		r          = RadialHeightOfCosTheta( k, costheta);										// Find the radius of this ray and value of cos(theta)
		z          = RadiusToAltitude( r );														// Convert to an altitude
		n          = m_RI.RefractiveIndex      ( z );											// and get the refractive index in the middle of the section
		dndr       = (m_RI.RefractiveIndex( z+100) -  m_RI.RefractiveIndex( z-100))/200.0;		// Get the refractive index gradient, this could be improved by using the refractivity
		n2costheta = n*n*costheta;																// calculate a few
		denom      = (k*dndr + n2costheta);														// intermediate terms
		dl         = n*k*deltheta/( costheta*denom );											// calculate the change in path length
		dpsi       = n2costheta*deltheta/denom;													// calculate the change in angular distance
		*L        += dl;																		// sum up the path lengths
		*PSI      += dpsi;																		// sum up the angular distances.
		thetastart = thetaend;																	// and the end theta becomes the new start theta.
	}
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayRefracted_ThomPepSim::TraceRayInNadir		2014-2-7*/
/** Traces rays from the observer down to the ground**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayRefracted_ThomPepSim::TraceRayInNadir(  double									robs,					//!< radius of the observer
													   SKTRAN_RayRefracted_TrajectoryData*		trajectory) const		//!< Object to store the ray trajectory parameters
{
	size_t				nlevel;
	double				range;
	size_t				i;
	double				r;
	bool				ok;
	bool				ok1;
	double				hobs;

	hobs = RadiusToAltitude(robs);
	ok   =  m_raytracingshells->IndexOfPointBelow( hobs, &nlevel);							// Get the first height below the observer
	if (!ok)																				// if the observer is at or below the ground
	{																						// then there is nothing to trace.
		ok = trajectory->ReserveSpace(0);												
	}
	else
	{
		ok   = ok && trajectory->ReserveSpace( nlevel+3 );									// reserve space for the trajectory
		ok   = ok && trajectory->PushBack    ( robs,  0.0, 0.0 );							// The first point is the location of the observer
		for ( i=nlevel+1; i > 0; )
		{
			--i;
			r     =  AltitudeToRadius( m_raytracingshells->At(i));							// Get the top of the Jth layer
			range = (robs-r); 
			NXASSERT(range >= 0.0);
			ok1   = trajectory->PushBack(r, 0.0, range);
			ok = ok && ok1;
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayRefracted_ThomPepSim::TraceRayInNadir		2014-2-7*/
/** Traces rays from the observer up to the top of the atmosphere **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayRefracted_ThomPepSim::TraceRayInZenith(  double									robs,					//!< radius of the observer
													   SKTRAN_RayRefracted_TrajectoryData*		trajectory) const		//!< Object to store the ray trajectory parameters
{
	size_t				nlevel;
	size_t				maxlevel;
	double				range;
	size_t				i;
	double				r;
	bool				ok;
	bool				ok1;
	double				hobs;
	SKTRAN_GridIndex	startidx;

	hobs     = RadiusToAltitude(robs);									// Get the altitude of the observer
	maxlevel = m_raytracingshells->NumShells();										// get the maximum number of shells
	ok = m_raytracingshells->IndexOfPointEqualOrAbove( hobs, &startidx);			// Find the shell above of equal to to the observer
	if (!ok)																		// Observer is above the atmosphere
	{																				// so there 
		ok = trajectory->ReserveSpace(0);											// is nothing to trace
	}
	else
	{
		range = m_raytracingshells->At(startidx)-robs;								// see if we have the "Equals condition
		NXASSERT( range >= 0);														// Quick reality check, should always be true
		if (range < 0.0001) ++startidx;												// If we are equal then skip the first level
		nlevel = maxlevel - startidx;
		ok   = ok && trajectory->ReserveSpace( nlevel+3 );							// reserve space for the trajectory
		ok   = ok && trajectory->PushBack    ( robs,  0.0, 0.0 );					// The first point is the location of the observer
		for ( i=startidx; i < maxlevel; i++ )
		{
			r     = AltitudeToRadius( m_raytracingshells->At(i));			// Get the top of the Jth layer
			range = (r-robs); 
			NXASSERT(range >= 0.0);
			ok1   = trajectory->PushBack(r, 0.0, range);
			ok    = ok && ok1;
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayRefracted_ThomPepSim::TraceRayOutsideAtmosphere		2014-1-29*/
/** Traces a ray through a spherical atmosphere. The ray tracing is very dependent
 *	upon Bouger's law, n.r.sin(theta) = k and this law becomes full of singularities
 *	if the ray is exactly straight down (theta = 0.0) . The code works well for
 *	all situations where the rays are not exactly straight down. The straight down
 *	condition is simply a straight line from the observer to the ground.
**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayRefracted_ThomPepSim::TraceRayOutsideAtmosphere(  double										robs,					//!< radius of the observer
																 double										rtoa,					//!< radius of the top of the atmosphere
																 double										nadang,					//!< the angle of the ray (from the observer) from the nadir at the observer	
																 double										rt,						//!< The tangent height of the ray
																 double										minradius,				//!< Minimum radius of the ray (either ground or tangent point)
																 bool										groundishit,			//!< true if the ground is hit
																 SKTRAN_RayRefracted_TrajectoryData*		trajectory) const		//!< Object to store the ray trajectory parameters
{
	double		R20km;
	double		nt;
	double		arg;
	double		alpha1;
	double		psi;
	double		dpsi;
	double		range;
	size_t		nlevel;
	size_t		numsteps;
	bool		aboveminradius;
	size_t		i;
	size_t		j;
	double		RJ1;
	double		RJ;
//	double		PL;
	double		PL2;
	bool		ok;
	bool		ok1;

	if (fabs(nadang) < 1.0E-6)																// Are we looking straight down, within 1 millionth of a degree
	{																						// Yup
		TraceRayInNadir(  robs, trajectory);											// So trace the ray down to the ground
	}
	else
	{
		nlevel  = m_raytracingshells->NumCells();
		R20km   = AltitudeToRadius(20000.0);												// Get the radius of 20 km
		nt      = m_RI.RefractiveIndex( RadiusToAltitude(rt) );								// Get the refractive index atthetangent point
		arg     = robs*sin(nadang)/rtoa;													// from law of sines with refractive index equal to 1 get the sin( Nadir angle at the top of the atmosphere
		alpha1  = asin(arg);																// Get the ANGLE OF INCIDENCE
		psi     = alpha1 - nadang;															// Angular distance of observer to entry point into atmosphere
		arg     = robs*robs + rtoa*rtoa - 2.0*robs*rtoa*cos(psi);							// cosine rule, square of distance from observer to entry point in the atmosphere
		if (arg < 0.0) arg=0.0;
		range = sqrt(arg);																	//  Distance from observer to entry point into the atmosphere.

		ok =       trajectory->ReserveSpace( 3*nlevel );									// reserve space for the trajectory
		ok = ok && trajectory->PushBack    ( robs,  0.0, 0.0 );								// The first point is the location of the observer
		ok = ok && trajectory->PushBack    ( rtoa,  psi, range);							// The second point is the entrance to the top of the atmosphere.

		numsteps  = 10;																						// use 10 point integration at higher altitudes
		aboveminradius = true;																
		for ( i=0; (i < nlevel) && aboveminradius; i++)
		{
			j   = nlevel-i;																					// Start at the top of the atmosphere
			RJ1 = AltitudeToRadius( m_raytracingshells->At(j));												// Get the top of the Jth layer
			RJ  = AltitudeToRadius( m_raytracingshells->At(j-1));											// Get the bottom of the Jth layer.
			if (RJ <= R20km) numsteps=20;																	// Increase number of steps for the lower atmosphere.
			aboveminradius = ((RJ - minradius) > 0.00001);													// Is the bottom of our layer below the tangent altitude
			if (!aboveminradius)																			// If we are at or below the tangent point
			{																								// then
				RJ        = minradius;																		// the bottom radius of this layer is the tangent layer
				numsteps  = 50;																				// and use 50 steps 
			}
//			PL  =        IntegrateCurvedPathLengthOld(  rt, nt,  RJ, RJ1, (int)numsteps);
			ok1 =        IntegrateCurvedPathLength   (  rt, nt,  RJ, RJ1,  numsteps, &PL2, &dpsi);
			ok1 = ok1 && trajectory->PushBack(RJ,  dpsi, PL2);	
			ok  = ok && ok1;
		}
		if (!groundishit) ok = ok && trajectory->MirrorPoints();						// Mirror the points from the tangent point out of the atmosphere if we dont hit the ground
		ok = ok && trajectory->SumAnglesAndLengths();
	}
	return true;
}



/*-----------------------------------------------------------------------------
 *					REFRAC		2014-1-23*/
/** This routine originated from the ACE-FTS refrac.F ray tracing code. This code
 *	has been shown to be very close to the ACE-FTS implementation. The code calculates
 *	The curved ray trajectory for a spherically symmetric atmosphere.
 *
 *	\param HOBS
 *	Height of the observer/ray origin above the geoid in meters
 *
 *	\param ZANGLE
 *	The zenith angle of the ray away from teh observer in degrees
 *
 *	\param trajectory
 *	returns the trajectory data as an array of points in an intermediate spherical shell system. The system
 *	returns the angular distance of the trajectory point from the observer and the radius of the
 *	point. The user can then use trigonometry to calculate the points in a 3-D system.
**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayRefracted_ThomPepSim::REFRAC (	double									HOBS,			// The height of the observer in meters
											    double									ZANGLE,			// The zenith angle of the ray away from the observer in degrees
												SKTRAN_RayRefracted_TrajectoryData*		trajectory) const
{
	bool			ok;
	double			RT;
	double			ZENANG;
	double			NADANG;
	double			ZMAX;
	double			RSPEC;
	double			RTOTL;
	double			maxnadir;
	double			minradius;
	bool			groundishit;

	if (ZANGLE > 180.0) ZANGLE = 360.0-ZANGLE;													// Ensure zenith angle is in range -180 to + 180
	ZENANG = nxmath::DegreesToRadians(ZANGLE);													// ZENITH ANGLE in radians
	NADANG = nxmath::DegreesToRadians( 180.0-ZANGLE);											// NADIR ANGLE in radians
	ZMAX   = m_raytracingshells->HighestShell();												// Height of top of the atmosphere.
	RSPEC  = AltitudeToRadius( HOBS );															// OBSERVER DIST. FROM CENTER of Earth (spherical earth approx)
    RTOTL  = AltitudeToRadius( ZMAX );															// Distance. To top OF atmosphere
	
	if (HOBS > ZMAX)																			// Observer is outside atmosphere
	{																							// then
		maxnadir = asin(RTOTL/RSPEC);															// Get the maximum nadir angle to the top of the atmosphere
		if (NADANG >= maxnadir)																	// If our nadir angle is greater
		{																						// then this line of sight does not hit the atmosphere so
			ok = trajectory->ReserveSpace(0);													// clear the trajectory
		}																						// and that is that
		else																					// otherwise
		{																						// this ray goes into the amosphere
			ok =       ComputeApproximateTangentRadius( HOBS, RSPEC, NADANG, &RT);				// Get a nominal starting guess tangent location
			ok = ok && FindGlobalTangentPoint( &RT, &minradius, &groundishit );					// Iteratively find the global tanegnt point	
			NXTRACE(("Osculating sphere coords, Ray tangent altitude = %f12.3, Observer altitude = %12.f\n", (double)RT-m_zerokmradius, (double)HOBS ));
			ok = ok && TraceRayOutsideAtmosphere( RSPEC, RTOTL, NADANG, RT, minradius, groundishit, trajectory);		// and find tehtrajectory points
		}
	}
	else
	{
		trajectory->ReserveSpace(0);											// clear the trajectory
		ok = false;
		nxLog::Record(NXLOG_WARNING, "SKTRAN_RayRefracted_ThomPepSim::REFRAC, Does not yet support ray tracing with observer inside atmosphere");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayRefracted_ThomPepSim::FindObserverAndLookFromTangentAltitude		 2014- 5- 22*/
/** A piece of code to help set up observer and lines of sight for the ACE-FTS system.
 *	The ACE-FTS, and maybe other occulters, determine tangent altitude from pressure broadening
 *	of spectra. We then have to back trace to find the observer locations and look direction.
 *	The nadir angle at the obsever is relatively easy as it follows Bougers law. We then use 
 *	straight line geometry to find the approximate observer location and then get the look vector.
 *	The observers location is then translated back to the regular geoid from the osculating sphere.
 *
 *	The final observers location is not completely consistent with the requested tangent point location
 *	as we use straight line geometry to find the observers location. This will result in the curved ray
 *	trajectory drawn in the occulation engine having a slight angular displacement although it will
 *	have an identical tangent altitude. Given that our goal in the occultation engine is to only support
 *	spherically symmetric conditions this will not generate any problems.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayRefracted_ThomPepSim::FindObserverAndLookFromTangentAltitude( double tanalt, double observeralt, const nxVector& usersun, nxVector* observer, nxVector* look)
{
	nxVector	refpt;
	nxVector	refptup;
	nxVector	refpthz;
	nxVector	obs;
	nxVector	up;
	nxVector	hz;
	nxVector	sun(usersun);
	PlanetSun	sunpos;
	double		rt;
	double		robs;
	double		k;
	double		nadang;
	double		nt;
	double		l;
	double		costheta;
	double		sintheta;
	bool		ok;

	ok    = ( m_coordinates != nullptr ) && (m_raytracingshells != nullptr) && (observeralt > tanalt) && sun.IsValid();
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_RayRefracted_ThomPepSim::FindObserverAndLookFromTangentAltitude, Either the ray tracer is uninitialized, the sun is invalid or the observer's altitude (%g) is greater than the tangent altitude (%g)", (double)observeralt, (double)tanalt);
		look->SetInvalid();
		observer->SetInvalid();
	}
	else
	{

		nt     = m_RI.RefractiveIndex( tanalt );							// Get the refractive index at the tangent point
		rt     = AltitudeToRadius( tanalt      );							// Get the radius of the tangent point
		k      = nt*rt;														// Get the ray constant at the tangent point sin(nadang)  = 1.0
		robs   = AltitudeToRadius( observeralt );							// Get the radius of the observer
		nadang = 180.0 - nxmath::asind( k/robs );							// Get the nadir angle at the observer (it will be greater than 90.0)
		l      = sqrt( robs*robs - rt*rt);									// Get the (approximate)  straight line distance from tangent point to observer.

		refpt   = m_coordinates->ReferencePointUnitVector() * rt;			// Get the location of the reference point in the osculating sphere system
		refptup = m_coordinates->ReferencePointUnitVector();				// Get the up vector at reference point in osculating sphere system
		refpthz = sun.ComponentPerpendicularTo( refptup ).UnitVector();		// Get the component of the sun vector perpendicular to "up" at the tangent point
		obs     = refpt - l*refpthz;										// Get the nominal location of observer in osculating sphere system.

		costheta = rt/robs;													// Theta = (approx) angular distance from reference point to observer
		sintheta = l/robs;
		up = costheta*refptup - sintheta*refpthz;							// Get the "up" unit vector at the observers location
		hz = sintheta*refptup + costheta*refpthz;							// Get the horizontal component at the observer's location
		*look     = nxmath::cosd(nadang)*up + nxmath::sind(nadang)*hz;		// Get the look vector at the observers location
		*observer = m_coordinates->TranslateOsculatingSphereToGeoid( obs );		 
	}
	return ok;
}