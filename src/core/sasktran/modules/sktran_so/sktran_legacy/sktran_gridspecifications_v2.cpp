#include "../sasktranv21_internals.h"
#include "../sasktran_legacy21.h"



/*-----------------------------------------------------------------------------
 *					SKTRANSO_SpecificationsUser_Legacy::SKTRANSO_SpecificationsUser_Legacy		2008-1-11*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_SpecificationsUser_Legacy::SKTRANSO_SpecificationsUser_Legacy()
                                    : m_diffusespecs(&m_rayregionmanager)
{
	m_atmosphericemissions_enabled = false;
	ConfigureEvenSpacedShells( 0.0, 1000.0, 100000.0 );
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_SpecificationsUser_Legacy::~SKTRANSO_SpecificationsUser_Legacy		2008-1-11*/
/** Release all of the usual stuff**/
/*---------------------------------------------------------------------------*/

SKTRANSO_SpecificationsUser_Legacy::~SKTRANSO_SpecificationsUser_Legacy()
{
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_SpecificationsUser_Legacy::ConfigureEvenSpacedShells		2010-2-10*/
/** Configures even spaced spherical shells for the radiative transfer calculations.
 *	This actually configures the ray tracing shell boundaries, the diffuse point
 *	locations and the location of the optical properties
 *
 *	The diffuse points are located in the mid-point of the ray tracing shells
  *	and the optical properties are located at both the mid point and shell boundaries.
*/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_SpecificationsUser_Legacy::ConfigureEvenSpacedShells( double minshellheight_meters, double shellwidth, double maxshell)
{
	bool						ok;
	nx1dArray<double>			shellalts;
	nx1dArray<double>			diffusealts;
	nx1dArray<double>			optpropheights;
	size_t						npts;

	npts                          =  (size_t) ((maxshell - minshellheight_meters)/shellwidth + 1);
	shellalts.Indgen(npts)        *=    shellwidth;							// Ray Shell Altitudes at 1 km boundaries to 100 km, 0.0, 1.0, 2.0, 3.0 ... 100.0
	diffusealts.Indgen(npts)	  *=    shellwidth;							// Diffuse Altitude    at 0.0 and then 0.5, 1.5, 2.5, 3.5, 99.5
	diffusealts		             -= 0.5*shellwidth;							// Offset to the 1/2 km point
	diffusealts[0]		          = 0.0;							// Make sure the bottom point is at 0 km

	optpropheights.Indgen(2*npts) *= (0.5*shellwidth);							// Optical properties Shell altitudes at 500 m intervals too 100.5 kms

	shellalts		+= minshellheight_meters;
	diffusealts		+= minshellheight_meters;
	optpropheights  += minshellheight_meters;

	ok    =       m_rayregionmanager.SetGroundAltitude(minshellheight_meters);
	ok    = ok && m_raytracingspecs.ConfigureRayTracingShellAlts     ( shellalts.UnsafeArrayBasePtr(),        shellalts.size()   );
	ok    = ok && m_diffusespecs.ConfigureDiffuseAltitudeResolution  ( diffusealts.UnsafeArrayBasePtr(),      diffusealts.size() );
	ok    = ok && m_opticalpropspecs.ConfigureOpticalPropertyShells	 ( optpropheights.UnsafeArrayBasePtr(),   optpropheights.size() );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRANSO_SpecificationsUser_Legacy::ConfigureEvenSpacedShells, Error configuring the default specifications. This is a problem");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_SpecificationsUser_Legacy::ConfigureUserDefinedShells		2013-6-20*/
/** Configures the shell specifications using the shells defined by the user. The
 *	shellalts parameter must pass in the altitudes of the shells in ascending order
 *	in meters.
 *	The diffuse points are placed in the middle of the shells and the optical
 *	properties are placed in the middle and on the shell
 */
/*---------------------------------------------------------------------------*/

bool SKTRANSO_SpecificationsUser_Legacy::ConfigureUserDefinedShells( const std::vector<double>& shellalts )
{
	size_t					numshells;
	double					minshellheight;
	double					maxshellheight;
	std::vector<double>		diffusealts;
	std::vector<double>		optpropheights;		
	double					lastalt;
	size_t					idx;
	bool					ok;

	ok = true;
	lastalt = shellalts.front() -1.0;
	for (idx = 0; idx < shellalts.size(); idx++)
	{
		ok = ok && (shellalts.at(idx) > lastalt);
		lastalt = shellalts.at(idx);
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRANSO_SpecificationsUser_Legacy::ConfigureUserDefinedShells, The user defined ray tracing shells must be unique values in ascending order");
	}
	else
	{

		numshells      = shellalts.size();
		minshellheight = shellalts.front();
		maxshellheight = shellalts.back();

		diffusealts.assign    ( numshells,      0.0 );
		optpropheights.assign ( 2*numshells-1, 0.0 );

		diffusealts.at(0) = minshellheight;
		for (idx = 1; idx < numshells; idx++ ) 
		{
			diffusealts.at(idx) = 0.5*(shellalts.at(idx-1) + shellalts.at(idx));
		}

		for (idx = 0; idx < (numshells-1); idx++)
		{
			optpropheights.at(2*idx) = shellalts.at(idx);
			optpropheights.at(2*idx+1) = 0.5*(shellalts.at(idx+1) + shellalts.at(idx));
		}
		optpropheights.at(2*idx) = shellalts.back();


		ok    =       m_rayregionmanager.SetGroundAltitude(minshellheight);
		ok    = ok && m_rayregionmanager.SetLowerBoundAltitude( minshellheight);
		ok    = ok && m_rayregionmanager.SetUpperBoundAltitude( maxshellheight);
		ok    = ok && m_raytracingspecs.ConfigureRayTracingShellAlts     ( &shellalts.front(),        shellalts.size()   );
		ok    = ok && m_diffusespecs.ConfigureDiffuseAltitudeResolution  ( &diffusealts.front(),      diffusealts.size() );
		ok    = ok && m_opticalpropspecs.ConfigureOpticalPropertyShells	 ( &optpropheights.front(),   optpropheights.size() );

		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRANSO_SpecificationsUser_Legacy::ConfigureUserDefinedShells, Error configuring the user defined shell specifications. This is a problem");
		}
	}
	return ok;



}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_SpecificationsUser_Legacy::CheckSpecsAreProperlyDefined		2010-5-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_SpecificationsUser_Legacy::CheckSpecsAreProperlyDefined() const
{
	bool ok;

	ok =    m_diffusespecs.IsProperlyDefined()
		 && m_rayregionmanager.IsProperlyDefined();
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRANSO_SpecificationsUser_Legacy::UpdateUndefinedParametersFromLinesOfSight		2010-6-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_SpecificationsUser_Legacy::UpdateUndefinedParametersFromLinesOfSight( const SKTRAN_LineOfSightArray_V21& linesofsight )
{
	return m_rayregionmanager.UpdateUndefinedParametersFromLinesOfSight(linesofsight);
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_SpecificationsUser_Legacy::GetCoordinateSystem		2010-5-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_SpecificationsUser_Legacy::CreateCoordinateSystem( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>* usercoords, double groundht, double toaht ) const
{
	bool	ok;

	ok = m_rayregionmanager.MakeCoordinateSystem( usercoords, groundht, toaht);
	return ok;;
}




