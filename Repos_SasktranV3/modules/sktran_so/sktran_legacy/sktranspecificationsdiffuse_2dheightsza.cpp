#include "../sasktranv21_internals.h"
#include "../sasktran_legacy21.h"
#include "sktran_legacy_internals.h"

#if defined(min)
#undef min
#endif

#if defined(max)
#undef max
#endif


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_Diffuse_Legacy		2010-5-19*/
/** This is the class used to specify the diffuse profiles. This class is stored
 *	internally and is created by #SKTRAN_SpecsUser_Diffuse_Legacy. This is done
 *	a) so the user does not have to worry about lifetime and object management of
 *	the specifications and b) so we can tweak any user settings before making the
 *	the actual settings used by the model.
 *
 *	Ground Points update
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_SpecsInternal_Diffuse_Legacy : public SKTRAN_SpecsInternal_Diffuse_V21
{
	private:
		double													m_maxdiffusealtitudealongincomingrays;
		SKTRAN_GridDefScatterAngle_V21*							m_gridScatterAngle;					//!< The grid used for caching scattering angle properties.
		SKTRAN_GridDefCosSZA_V21*								m_solartransmissionSzaGrid;			//!< The CosSZA grid used to calculate solar transmission to all elements of all lines of sight
		SKTRAN_GridDefSLON_V21*									m_solartransmissionSlonGrid;
		SKTRAN_GridDefCosSZA_V21*								m_diffuseSzaGrid;					//!< The CosSZA grid used to calculate solar transmission to all elements of all lines of sight
		SKTRAN_GridDefSLON_V21*									m_diffuseSlonGrid;
		SKTRAN_UnitSphere_V2*									m_unitsphere;						//!< The unit sphere used for outgoing radiances
//		std::vector< SKTRAN_GridDefDiffuseIncomingZenith_V2* >	m_incomingzenithangle;				//!< An array [numdiffuseheights] of diffuse zenith angles as a function of heights from m_diffuseradii
		std::vector< SKTRAN_UnitSphereLatLonGrid* >				m_incomingunitspheres;				//!< An array [numdiffuseheights] of incoming unit spheres
//		SKTRAN_GridDefDiffuseIncomingAzimuth_V2*				m_incomingazimuthangle;				//!< The incoming azimuth angles usined for rays coming into a point
		SKTRAN_GridDefDiffuseHeights_V21*						m_diffuseradii;						//!< The diffuse radii used to calculate the location of diffuse profile shells
		SKTRAN_TableRayLOSFactory_Legacy*						m_lossinglescatterfactory;			//!< The factory object used to create the internal single scatter tables used by observer LOS rays

	private:
		void													ReleaseResources	();

	public:
																SKTRAN_SpecsInternal_Diffuse_Legacy();
		virtual												   ~SKTRAN_SpecsInternal_Diffuse_Legacy();

		virtual const SKTRAN_GridDefScatterAngle_V21*			ScattterAngleGrid				() const override { return m_gridScatterAngle;}

		virtual const SKTRAN_GridDefCosSZA_V21*					DiffuseProfileCosSZA			( ) const override;
		virtual const SKTRAN_GridDefSLON_V21*					DiffuseProfileSLON				( ) const override;
		virtual const SKTRAN_GridDefDiffuseHeights_V21*			DiffuseHeights					( size_t profileidx) const override;
		virtual double 											MaxDiffuseAltitude				( size_t profileindex) const override;

		virtual const SKTRAN_GridDefCosSZA_V21*					SolarTransmissionCosSZA			( ) const override;
		virtual const SKTRAN_GridDefSLON_V21*					SolarTransmissionSLON			( ) const override;
		virtual const SKTRAN_TableRayLOSFactory*				LOSSingleScatterTableFactory	( ) const override;
	    
//		virtual const SKTRAN_GridDefDiffuseIncomingZenith_V2*	IncomingZenith					( size_t profileindex, size_t heightindex ) const;
//	    virtual const SKTRAN_GridDefDiffuseIncomingAzimuth_V2*	IncomingAzimuth					( size_t profileindex, size_t heightindex ) const;
		virtual const SKTRAN_UnitSphereLatLonGrid*				IncomingUnitSphere				( size_t profileindex, size_t heightindex ) const override;

		virtual const SKTRAN_UnitSphere_V2*						OutboundUnitSphere				( size_t profileindex, size_t heightindex ) const override;
		virtual bool											FirstTimeInitializeDiffuseObjects (	) override;

		virtual bool											CreateEmptyDiffuseTables		( SKTRANSO_TableDiffusePoints**	   diffusepointstable,
																								  SKTRANSO_TableSolarTransmission**  solartranstable,
																								  SKTRANSO_TableGroundPointDiffuse** groundpointstable,
																								  SKTRANSO_InternalEmissionPropertiesTable**		 emissionstable) const override;


	friend class SKTRAN_SpecsUser_Diffuse_Legacy;
};



/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_Diffuse_Legacy::SKTRAN_SpecsInternal_Diffuse_Legacy		2010-5-19*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_SpecsInternal_Diffuse_Legacy::SKTRAN_SpecsInternal_Diffuse_Legacy()
{
	m_gridScatterAngle         = NULL;					// The grid used for caching scattering angle properties.
	m_solartransmissionSzaGrid = NULL;					// The CosSZA grid used to calculate solar transmission to all elements of all lines of sight
	m_solartransmissionSlonGrid= NULL;					// The CosSZA grid used to calculate solar transmission to all elements of all lines of sight
	m_diffuseSzaGrid           = NULL;
	m_diffuseSlonGrid          = NULL;
//	m_incomingazimuthangle     = NULL;					// The incoming azimuth angles usined for rays coming into a point
	m_unitsphere               = NULL;
	m_lossinglescatterfactory  = NULL;
	m_maxdiffusealtitudealongincomingrays = -99999.0;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_Diffuse_Legacy::~SKTRAN_SpecsInternal_Diffuse_Legacy		2010-5-21*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_SpecsInternal_Diffuse_Legacy::~SKTRAN_SpecsInternal_Diffuse_Legacy()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_Diffuse_Legacy::ReleaseResources		2010-5-19*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_SpecsInternal_Diffuse_Legacy::ReleaseResources()
{
	size_t			numheights;

	numheights   = m_incomingunitspheres.size();
	for (size_t hidx = 0; hidx < numheights; hidx++)
	{
//		if (m_incomingzenithangle[hidx] != NULL)  m_incomingzenithangle[hidx]->Release();
		if (m_incomingunitspheres[hidx] != NULL)  m_incomingunitspheres[hidx]->Release();
	}
//	m_incomingzenithangle.clear();
	m_incomingunitspheres.clear();

	if (m_gridScatterAngle         != NULL) m_gridScatterAngle->Release();
	if (m_solartransmissionSzaGrid != NULL) m_solartransmissionSzaGrid->Release();
	if (m_solartransmissionSlonGrid!= NULL) m_solartransmissionSlonGrid->Release();
//	if (m_incomingazimuthangle     != NULL) m_incomingazimuthangle->Release();
	if (m_diffuseradii             != NULL) m_diffuseradii->Release();
	if (m_diffuseSzaGrid           != NULL) m_diffuseSzaGrid->Release();
	if (m_diffuseSlonGrid          != NULL) m_diffuseSlonGrid->Release();
	if (m_unitsphere               != NULL) m_unitsphere->Release();
	if (m_lossinglescatterfactory  != NULL) m_lossinglescatterfactory->Release();

	m_gridScatterAngle         = NULL;
	m_solartransmissionSzaGrid = NULL;
	m_solartransmissionSlonGrid= NULL;
//	m_incomingazimuthangle     = NULL;
	m_diffuseradii             = NULL;
	m_diffuseSzaGrid           = NULL;
	m_diffuseSlonGrid          = NULL;
	m_unitsphere               = NULL;
	m_lossinglescatterfactory  = NULL;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_Diffuse_Legacy::DiffuseProfileCosSZA		2010-4-28*/
/** **/
/*---------------------------------------------------------------------------*/

const SKTRAN_GridDefCosSZA_V21* SKTRAN_SpecsInternal_Diffuse_Legacy::DiffuseProfileCosSZA () const
{
	bool	ok;

	ok = (m_diffuseSzaGrid != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_SpecsUser_Diffuse_Legacy::DiffuseProfileSLON, the configuration is dirty or is NULL");
	}
	return m_diffuseSzaGrid;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_Diffuse_Legacy::DiffuseProfileSLON		2010-4-28*/
/** **/
/*---------------------------------------------------------------------------*/

const SKTRAN_GridDefSLON_V21* SKTRAN_SpecsInternal_Diffuse_Legacy::DiffuseProfileSLON () const
{
	bool	ok;

	ok = (m_diffuseSlonGrid != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_SpecsUser_Diffuse_Legacy::DiffuseProfileSLON, the configuration is dirty or is NULL");
	}
	return m_diffuseSlonGrid;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_Diffuse_Legacy::GetSolarTransCosSzaGrid		2008-1-15*/
/** **/
/*---------------------------------------------------------------------------*/

const SKTRAN_GridDefCosSZA_V21*	SKTRAN_SpecsInternal_Diffuse_Legacy::SolarTransmissionCosSZA( ) const
{
	bool	ok;

	ok = (m_solartransmissionSzaGrid != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_SpecsInternal_Diffuse_Legacy::SolarTransmissionCosSZA, the configuration is dirty or is NULL");
	}
	return m_solartransmissionSzaGrid;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_Diffuse_Legacy::GetSolarTransCosSzaGrid		2008-1-15*/
/** **/
/*---------------------------------------------------------------------------*/

const SKTRAN_GridDefSLON_V21*	SKTRAN_SpecsInternal_Diffuse_Legacy::SolarTransmissionSLON( ) const
{
	bool	ok;

	ok = (m_solartransmissionSlonGrid != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_SpecsUser_Diffuse_Legacy::SolarTransmissionSLON, the configuration is dirty or is NULL");
	}
	return m_solartransmissionSlonGrid;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_Diffuse_Legacy::GetDiffuseHeightsGrid		2008-1-15*/
/** **/
/*---------------------------------------------------------------------------*/

const SKTRAN_GridDefDiffuseHeights_V21*	SKTRAN_SpecsInternal_Diffuse_Legacy::DiffuseHeights( size_t profileidx) const
{
	bool	ok;

	ok = (m_diffuseradii != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_SpecsInternal_Diffuse_Legacy::DiffuseHeights, the configuration is dirty or is NULL");
	}
	return m_diffuseradii;			// Ignore the profile index 
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_Diffuse_Legacy::MaxDiffuseAltitude		2010-4-28*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_SpecsInternal_Diffuse_Legacy::MaxDiffuseAltitude( size_t profileindex) const
{
	return m_maxdiffusealtitudealongincomingrays;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_Diffuse_Legacy::IncomingZenith		2010-4-28*/
/** **/
/*---------------------------------------------------------------------------*/

const SKTRAN_UnitSphereLatLonGrid* SKTRAN_SpecsInternal_Diffuse_Legacy::IncomingUnitSphere( size_t profileindex, size_t heightindex ) const
{
	bool	ok;

	ok = (heightindex < m_incomingunitspheres.size() );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_SpecsInternal_Diffuse_Legacy::IncomingZenith, the configuration is dirty, pointers are NULL or heightindex is out of range");
		if (heightindex >= m_diffuseradii->NumAltitudes()) heightindex = 0;
	}
	return (m_incomingunitspheres[heightindex]);
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_Diffuse_Legacy::IncomingZenith		2010-4-28*/
/** **/
/*---------------------------------------------------------------------------*/

/*
const SKTRAN_GridDefDiffuseIncomingZenith_V2* SKTRAN_SpecsInternal_Diffuse_Legacy::IncomingZenith( size_t profileindex, size_t heightindex ) const
{
	bool	ok;

	ok = (heightindex < m_incomingunitspheres.size() );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_SpecsInternal_Diffuse_Legacy::IncomingZenith, the configuration is dirty, pointers are NULL or heightindex is out of range");
		if (heightindex >= m_diffuseradii->NumAltitudes()) heightindex = 0;
	}
	return (m_incomingunitspheres[heightindex]);
}
*/

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_Diffuse_Legacy::IncomingAzimuth		2010-4-28*/
/** **/
/*---------------------------------------------------------------------------*/

/*
const SKTRAN_GridDefDiffuseIncomingAzimuth_V2*	 SKTRAN_SpecsInternal_Diffuse_Legacy::IncomingAzimuth( size_t profileindex, size_t heightindex ) const
{
	bool	ok;

	ok = (m_incomingazimuthangle != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_SpecsInternal_Diffuse_Legacy::IncomingAzimuth, the configuration is dirty or pointers are NULL");
	}
	return m_incomingazimuthangle;
}
*/


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_Diffuse_Legacy::OutboundUnitSphere		2010-4-28*/
/** **/
/*---------------------------------------------------------------------------*/

const SKTRAN_UnitSphere_V2* SKTRAN_SpecsInternal_Diffuse_Legacy::OutboundUnitSphere( size_t profileindex, size_t heightindex ) const
{
	bool	ok;

	ok = (m_unitsphere != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_SpecsInternal_Diffuse_Legacy::OutboundUnitSphere, the configuration is dirty or unit sphere is NULL");
	}
	return m_unitsphere;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_Diffuse_Legacy::FirstTimeInitializeDiffuseObjects		2010-6-24*/
/** Provides an opportunity to initialize any internal objects before starting the
 *	multiple scatter calculations. 
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsInternal_Diffuse_Legacy::FirstTimeInitializeDiffuseObjects (	)
{
	bool	ok;
	
	ok = true;
	//ok = m_unitsphere->Initialize();
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING," SKTRAN_SpecsInternal_Diffuse_Legacy::FirstTimeInitializeDiffuseObjects, There was an error initializing the unit sphere. Thats not good");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_Diffuse_Legacy::LOSSingleScatterTableFactory		2010-6-22*/
/** **/
/*---------------------------------------------------------------------------*/

const SKTRAN_TableRayLOSFactory* SKTRAN_SpecsInternal_Diffuse_Legacy::LOSSingleScatterTableFactory( ) const
{
	return m_lossinglescatterfactory ;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_Diffuse_Legacy::CreateDiffusePointsTable		2010-5-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsInternal_Diffuse_Legacy::CreateEmptyDiffuseTables( SKTRANSO_TableDiffusePoints**	    diffusepointstable,
																			SKTRANSO_TableSolarTransmission**  solartranstable,
																			SKTRANSO_TableGroundPointDiffuse** groundpointstable,
																			SKTRANSO_InternalEmissionPropertiesTable**	 emissiontable) const
{
	bool	ok;
	SKTRAN_TableDiffusePoints_2D_Height_SZA*	diffusetable;

	diffusetable        = new SKTRAN_TableDiffusePoints_2D_Height_SZA;
	*diffusepointstable = diffusetable;
	*solartranstable    = new SKTRAN_TableSolarTransmission_2D_Height_SZA;
	*groundpointstable  = new SKTRAN_TableGroundPointDiffuse_Colocated(diffusetable);
	*emissiontable      = new SKTRANSO_InternalEmissionPropertiesTable_1D_Height;

	ok =    ( diffusetable       != NULL )
		 && ( *solartranstable   != NULL )
		 && ( *groundpointstable != NULL )
	     && (*emissiontable      != NULL );

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_SpecsUser_Diffuse_Legacy::CreateDiffuseTables, Error allocating memeory");
		if (*emissiontable      != NULL) delete (*emissiontable);
		if (*groundpointstable  != NULL) delete (*groundpointstable);
		if (*solartranstable    != NULL) delete (*solartranstable);
		if (*diffusepointstable != NULL) delete (*diffusepointstable);
		*emissiontable      = NULL;
		*groundpointstable  = NULL;
		*solartranstable    = NULL;
		*diffusepointstable = NULL;
	}
	else
	{
		(*emissiontable)->AddRef();
		(*diffusepointstable)->AddRef();
		(*solartranstable)->AddRef();
		(*groundpointstable)->AddRef();
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::SKTRAN_SpecsUser_Diffuse_Legacy		2010-4-28*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_SpecsUser_Diffuse_Legacy::SKTRAN_SpecsUser_Diffuse_Legacy( const SKTRAN_RayTracingRegionManager* rayregionmanager )
{
//	m_isdirty                             = true;
	m_degreesPerSolarTransSza             = 0.5;					// The number of SZA degrees for each solar transmission entry 
	//m_degreesPerDiffuseSZAProfile         = 20.0;					// The angular separation (degrees) of diffuse profiles in the Solar Zenith Angle  grid;= 0.015;
	m_numDiffuseProfiles							  = 1;
	m_maxdiffusealtitudealongincomingrays = 100000.0;		// By default, max altitude for integrating diffuse rays along source function is at 50 kms.
	m_useUserDefinedSZATransmissiongrid   = false;
//	m_useUserDefinedSZAGroundgrid         = false;
	m_groundRes                           = 5;						// If using "auto incoming zenith" then The number of incoming zenith angles used to cover the ground ((usually coarse to medium)
	m_horizonRes                          = 6;						// If using "auto incoming zenith" thenThe number of incoming zenith angles used to cover the horizon region (usually fine scale)
	m_atmosRes                            = 3;						// If using "auto incoming zenith" thenThe number of incoming zenith angles used to cover the zenith directions looking upward.
	m_unitsphere                          = NULL;
	m_rayregionmanager                    = rayregionmanager;
	m_use_losinternalsinglescattertables  = true;
	m_configuredForMC					  = false;

//	m_gridScatterAngle         = NULL;					// The grid used for caching scattering angle properties.
//	m_solartransmissionSzaGrid = NULL;					// The CosSZA grid used to calculate solar transmission to all elements of all lines of sight
//	m_solartransmissionSlonGrid= NULL;					// The CosSZA grid used to calculate solar transmission to all elements of all lines of sight
//	m_diffuseSzaGrid           = NULL;
//	m_diffuseSlonGrid          = NULL;
//	m_incomingzenithangle      = NULL;					// An array of diffuse zenith angles as a function of heights from m_diffuseradii

//	m_unitsphere               = new SKTRAN_UnitSphereME(169);		// The unit sphere used for outgoing radiances
//	m_unitsphere->AddRef();

	ConfigureDefaults();

}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::~SKTRAN_SpecsUser_Diffuse_Legacy		2010-4-28*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_SpecsUser_Diffuse_Legacy::~SKTRAN_SpecsUser_Diffuse_Legacy()
{
	if (m_unitsphere != NULL) m_unitsphere->Release();
	m_unitsphere = NULL;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureDefaults		2010-4-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureDefaults()
{
	static double g_defaultincomingzenith[] = { 15, 30, 40, 50, 60, 70, 75, 80, 85, 87, 89, 90, 91, 93, 95, 100, 105, 110, 120, 130, 140, 150, 165, 180};
	nx1dArray<double>								diffuseazi;
	size_t											num_incoming_azimuthangles	= 12;
	bool											ok;
	size_t											nzen;
	size_t											idx;
	nx1dArray<double>								zenangles;
	double											lastzen;
	double											angle;

	diffuseazi.Indgen(num_incoming_azimuthangles)   *=  (360.0/num_incoming_azimuthangles);	// Get the azimuth angles same as Sasktran V1											// num_incoming_azimuthangles from ModelParameters.txt


	nzen = N_ELEMENTS(g_defaultincomingzenith);						// Incoming zenith angles
	zenangles.SetSize(nzen);										// Use the non-even grid used in Sasktran Version 1, it works well
	lastzen = 0.0;
	for (idx = 0; idx < nzen; idx++)
	{
		angle = (g_defaultincomingzenith[idx]+lastzen)*(0.5);
		zenangles.At(idx) = angle;
		lastzen = g_defaultincomingzenith[idx];
	}

	ok    =		  ConfigureEvenSpacedShells				 ( 0.0, 1000.0, 100000.0 );
	ok    = ok && ConfigureCosScatteringAngleGrid		 ( 0.5 );
	ok    = ok && ConfigureIncomingAzimuthResolution     ( diffuseazi.UnsafeArrayBasePtr(),       diffuseazi.size()  );
	ok    = ok && ConfigureUserDefinedIncomingZenithAngle( zenangles.UnsafeArrayBasePtr(),        zenangles.size ()  );
	ok    = ok && ConfigureDegreesSZAPerSolarTransmission( 0.5 );
	ok    = ok && ConfigureNumberDiffuseProfiles		 ( 1 );
//	ok    = ok && ConfigureDegreesSZAPerDiffuseProfile   ( 20.0 );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_GridSpecificationsLegacy_V2::ConfigureDefaults, Error configuring the default specifications. This is a problem");
	}
	NXTRACE(( "SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureDefaults, reset m_use_losinternalsinglescattertables back to true\n"));
	m_use_losinternalsinglescattertables = true;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureEvenSpacedShells		2010-2-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureEvenSpacedShells( double minshellheight_meters, double shellwidth, double maxshell)
{
	bool						ok;
	nx1dArray<double>			diffusealts;
	size_t						npts;

	npts                          =  (size_t) ((maxshell - minshellheight_meters)/shellwidth + 1);
	diffusealts.Indgen(npts)	  *=    shellwidth;							// Diffuse Altitude    at 0.0 and then 0.5, 1.5, 2.5, 3.5, 99.5
	diffusealts		             -= 0.5*shellwidth;							// Offset to the 1/2 km point
	diffusealts[0]		          = 0.0;							// Make sure the bottom point is at 0 km
	diffusealts		             += minshellheight_meters;

	ok    = ConfigureDiffuseAltitudeResolution  ( diffusealts.UnsafeArrayBasePtr(),      diffusealts.size() );
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_GridSpecificationsLegacy_V2::ConfigureAltitudeResolution		2008-1-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureDiffuseAltitudeResolution( const double* altitudes, size_t npts )
{
	size_t	i;
	double	lastalt = -6500000.0;
	bool	ok = true;

	m_diffusealtitudes.clear();
	m_diffusealtitudes.reserve( npts );
	for (i = 0; i < npts; i++ )
	{
		ok = ok && (lastalt < altitudes[i] );
//		NXASSERT((ok));
		lastalt = altitudes[i];
		m_diffusealtitudes.push_back(  altitudes[i] );
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_GridSpecificationsLegacy_V2::ConfigureDiffuseAltitudeResolution, The diffuse altitudes are not in ascending order. This will create problems");
	}
	else
	{
		ok  = (m_diffusealtitudes.back() > 990);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_GridSpecificationsLegacy_V2::ConfigureDiffuseAltitudeResolution, The diffuse altitudes appear to be in kilometers, they should be expressed in meters, Please Check");
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureScatteringAngleGrid		2008-1-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureCosScatteringAngleGrid( double degrees_resolution )
{
	m_scatteringresolution_degrees = degrees_resolution;
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureZenithResolution		2008-1-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureIncomingZenithResolutions( size_t groundResolution, size_t horizonResolution, size_t atmosResolution, bool shiftedHorizon )
{
	m_shiftedHorizon = shiftedHorizon;
	m_groundRes =  groundResolution;
	m_horizonRes = horizonResolution;
	m_atmosRes   = atmosResolution;
	m_useUserDefinedincomingzenith = false;
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureUserDefinedZenithAngle		2008-3-12*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureUserDefinedIncomingZenithAngle( double* zenithdegrees, size_t numzen)
{
	size_t	idx;

	m_userdefinedincomingzenith.clear();
	m_userdefinedincomingzenith.reserve(numzen);
	for (idx = 0; idx < numzen; idx++)
	{
		NXASSERT(( (zenithdegrees[idx] >= 0.0) && (zenithdegrees[idx] < 180.1) ));	// Make sure angles are between 0 and pi 
		m_userdefinedincomingzenith.push_back(zenithdegrees[idx]);
	}
	m_useUserDefinedincomingzenith = true;
	return true;

}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureAzimuthResolution		2008-1-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureIncomingAzimuthResolution( const double* azimuths, size_t npts )
{
	size_t i;
	bool	ok = true;
	double	lastval = -99999;

	m_incomingazimuth.clear();
	m_incomingazimuth.reserve( npts );
	for (i = 0; i < npts; i++ )
	{
		ok = ok && ( lastval < azimuths[i] );
		NXASSERT((ok));
		lastval = azimuths[i];
		m_incomingazimuth.push_back(  azimuths[i] );
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_GridSpecificationsLegacy_V2::ConfigureIncomingAzimuthResolution, The incoming azimuths are not in ascending order. This will create problems");
	}
	else
	{
		ok  =    (m_incomingazimuth.front() >  -180)
			  && (m_incomingazimuth.back () >   200)
			  && (m_incomingazimuth.back () <= (360.0));
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_GridSpecificationsLegacy_V2::ConfigureIncomingAzimuthResolution, The incoming azimuths do not appear to be in degrees. Express azimuths in radians in the range 0 to Two Pi");
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureDegreesSZAPerDiffuseProfile		2008-1-14*/
/** **/
/*---------------------------------------------------------------------------*/

//bool SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureDegreesSZAPerDiffuseProfile( double diffusedeltasza_degrees )
//{
//	m_degreesPerDiffuseSZAProfile = diffusedeltasza_degrees;
//	return true;
//}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureDegreesSZAPerDiffuseProfile		2012-06-19*/
/** Replaces ConfigureDegreesSZAPerDiffuseProfile; ensures an odd number of profiles are used so that one occurs at the tangent point. **/
/*---------------------------------------------------------------------------*/
bool SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureNumberDiffuseProfiles(size_t numProfiles)
{
	NXASSERT(0<numProfiles);
	
	if((numProfiles%2) == 1){
		m_numDiffuseProfiles = numProfiles;
	} else{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureNumberDiffuseProfiles, number of profiles should be odd so that one occurs at the tangent point; decreasing numProfiles to %Iu", numProfiles-1);
		m_numDiffuseProfiles = numProfiles-1;
	}
		return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureDegreesSZAPerSolarTransmission		2008-1-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureDegreesSZAPerSolarTransmission( double deltasza_degrees )
{
	m_degreesPerSolarTransSza = deltasza_degrees;
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureMaxJsourceAltitudeAlongDiffuseRays		2009-2-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureMaxJsourceAltitudeAlongDiffuseRays( double heightm)
{
	m_maxdiffusealtitudealongincomingrays = heightm;
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureOutboundUnitSphere		2010-6-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureOutboundUnitSphere(  SKTRAN_UnitSphere_V2* sphere )
{
	if (sphere       != NULL) sphere->AddRef();
	if (m_unitsphere != NULL) m_unitsphere->Release();
	m_unitsphere  = sphere;
	return true;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::UseGlobalSolarTransmissionTableForSingleScatter		2010-6-23*/
/** Flags whether we use the global solar transmission table for line of sight
 *	single scatter or not. The default is false (ie each line of sight ray has its own single scatter table).
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::UseGlobalSolarTransmissionTableForSingleScatter	( bool useglobal)
{
	m_use_losinternalsinglescattertables = !useglobal;
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureUserDefinedTransmissionSZAGrid		2008-3-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureUserDefinedTransmissionSZAGrid( double* cossza, size_t numsza)
{
	size_t	idx;
	bool	ok;
	double	lastval = -2;

	m_szatransmission.clear();
	lastval = -2;
	ok      = true;
	for (idx = 0; idx < numsza; idx++)
	{
		m_szatransmission.push_back( cossza[idx] );
		ok = ok && (cossza[idx] > lastval);
		lastval = cossza[idx];
		NXASSERT(( (cossza[idx] >= -1) && (cossza[idx] <= +1.00) ));
	}
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING, "SKTRAN_GridSpecificationsLegacy_V2::ConfigureUserDefinedTransmissionSZAGrid, User defined grid is not is ascneding order. Please correct. returning to default settings");
	}
	m_useUserDefinedSZATransmissiongrid = ok;
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::MakeIncomingAzimuthGrid		2008-1-11*/
/** 2013-08-28 ndl
 *	Dan Zawada noted that there was no diffuse point on the ground that had downward looking incoming directions.
 *	This created problems when interpolating downward radiances between the ground point and the first
 *	diffuse point. This code fixes that problem.
 *
 *	The first point in the returned array is always a diffuse point on the
 *	ground that is conceptually on the ground with no downward looking directions: it is the
 *	same ground point as defined in previous code. The second point
 *	is now used to represent a diffuse point on the ground which is conceptually
 *	infinitessimally above the ground and has downward looking directions.
 *	This will then work as expected when interpolating radiance in various directions between ground
 *	and the first off ground diffuse point.
 *
 *	Note that this method works as the code in SKTRAN_GridDefBase interpolates points using std::upper_bound
 *	so the code will always index the "off ground" ground point when do ing interpolation.  Got to be careful that we dont change that code.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::MakeDiffuseRadiiGrid( SKTRAN_GridDefDiffuseHeights_V21** userdiffuseheights ) const
{
	bool								ok;
	size_t								outidx;
	size_t								inidx;
	SKTRAN_GridDefDiffuseHeights_V21*	diffuseheights;
	size_t								npts;
	bool								gotfullsphereground;

	diffuseheights = new SKTRAN_GridDefDiffuseHeights_V21;
	ok      = (diffuseheights != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_ERROR, "SKTRAN_SpecsUser_Diffuse_Legacy::MakeDiffuseRadiiGrid, Error allocating memory for azimuth grid");
	}
	else
	{
		diffuseheights->AddRef();
		npts = m_diffusealtitudes.size();
		ok = (npts > 1);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "SKTRAN_SpecsUser_Diffuse_Legacy::MakeDiffuseRadiiGrid, the diffuse profile array must have more than 1 point");
		}
		else
		{
			gotfullsphereground = m_diffusealtitudes.at(0) == m_diffusealtitudes.at(1); // See if we have the full sphere ground point. The first altitude is the ground point hemisphere, the second point is the ground point, full sphere ie infinitessimally above ground
			if (!gotfullsphereground) npts++;
			ok = diffuseheights->AllocateGridArray(npts);
			if (!ok)
			{
				nxLog::Record(NXLOG_ERROR, "SKTRAN_SpecsUser_Diffuse_Legacy::MakeDiffuseRadiiGrid, Error allocating memory for altitude grid");
			}
			else
			{
				diffuseheights->AtVar(0) = m_diffusealtitudes[0];		// This will be the ground point actually on the ground. Its a hemisphere.
				diffuseheights->AtVar(1) = m_diffusealtitudes[0];		// This will be the ground point infinitessimally above the ground. It is a full sphere.
				outidx = 2;										
				inidx  = gotfullsphereground ? 2 : 1;					// If the incoming ray already has the "elevated groud point" then start at 2. otehrwise start at first point really above the ground
				while (inidx < m_diffusealtitudes.size())
				{
					diffuseheights->AtVar(outidx++) = m_diffusealtitudes[inidx++];
				}
	//			ok = (m_diffuseradii->At(0) == m_raytracingshells->At(0));
	//			if (!ok)
	//			{
	//				nxLog::Record(NXLOG_WARNING, "SKTRAN_SpecsUser_Diffuse_Legacy::MakeDiffuseRadiiGrid, The lowest diffuse radii must be the same as the lowest ray tracing shell");
	//			}
			}
		}
	}
	*userdiffuseheights = diffuseheights;
	return ok;
}




/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::MakeDiffuseSzaGrid		2008-1-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::MakeDiffuseSzaGrid( SKTRAN_GridDefCosSZA_V21** userdiffuseSzagrid, SKTRAN_GridDefSLON_V21** userdiffuseslongrid ) const
{
	double						szaDelta;
	//size_t						szaNum;
	size_t						idx;
	size_t						n;
	double						offset;
	bool						ok;
	double						mincossza;
	double						maxcossza;
	double						sza;
	double						cossza;
	double						maxsza;
	double						minsza;
	SKTRAN_GridDefCosSZA_V21*	diffuseSzaGrid;
	SKTRAN_GridDefSLON_V21*		diffuseSLONgrid;

	ok = m_rayregionmanager->GetSZA( &sza, &minsza, &maxsza );
	cossza    = nxmath::cosd( sza );
	mincossza = nxmath::cosd( maxsza);
	maxcossza = nxmath::cosd( minsza);
	NXASSERT(( mincossza <= maxcossza ));
	//szaNum  = (size_t)((maxsza - minsza)/m_degreesPerDiffuseSZAProfile)+ 1;
	//printf("%u,",szaNum);
	if (m_numDiffuseProfiles > 1)
	{
		n         = (m_numDiffuseProfiles -1);
		offset    = n/2.0;							
		szaDelta  = (maxcossza - mincossza)/n;
	}
	else
	{
		szaDelta  = 0.0;
		offset    = 0.0;
	}

	diffuseSzaGrid = new SKTRAN_GridDefCosSZA_V21;
	diffuseSLONgrid = new SKTRAN_GridDefSLON_V21;

	ok = ( (diffuseSzaGrid != NULL) && (diffuseSLONgrid != NULL) );
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING, "SKTRAN_SpecsUser_Diffuse_Legacy::MakeDiffuseSzaGrid, error allocating SKTRAN_GridDefDiffuseProfileCosSZA_V2 object");
	}
	else
	{
		diffuseSzaGrid->AddRef();
		diffuseSLONgrid->AddRef();
		ok =       diffuseSzaGrid ->AllocateGridArray( m_numDiffuseProfiles );
		ok = ok && diffuseSLONgrid->AllocateGridArray( m_numDiffuseProfiles );
		if (!ok)
		{
			nxLog::Record( NXLOG_WARNING, "SKTRAN_SpecsUser_Diffuse_Legacy::MakeDiffuseSzaGrid, error allocaing %Iu points for SZA grid for solar geometry cache", (size_t)m_numDiffuseProfiles );
		}
		else
		{
			for (idx = 0; idx < m_numDiffuseProfiles; idx++)
			{
				diffuseSzaGrid->AtVar(idx) =  cossza + (idx-offset)*szaDelta;
				diffuseSLONgrid->AtVar(idx) = 0.0;
			}
		}
	}
	*userdiffuseSzagrid = diffuseSzaGrid;
	*userdiffuseslongrid = diffuseSLONgrid;
	return ok;
}
/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::MakeIncomingUNitSpheres		2008-1-14*/
/** This makes the incoming unit spheres for the diffuse points. 
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::MakeIncomingUnitSpheres( SKTRAN_GridDefDiffuseHeights_V21* internaldiffuseheights, std::vector< SKTRAN_UnitSphereLatLonGrid*>* inboundunitspheres, std::shared_ptr< const SKTRAN_CoordinateTransform_V2>& coordinatesystem ) const 
{
	size_t									numheights;
	bool										ok;
	size_t										hidx;
	bool										ok1;
	SKTRAN_UnitSphereLatLonGrid*				entry;
	SKTRAN_GridDefDiffuseIncomingZenith_V21*		zengrid     = NULL;
	SKTRAN_GridDefDiffuseIncomingAzimuth_V21*	azimuthgrid = NULL;

	NXTRACE_ONCEONLY(firstime, ("***** 2011-02-09 ***** Check  Addrefs and Releases in SKTRAN_SpecsUser_Diffuse_Legacy::MakeIncomingUnitSpheres\n"));

	numheights            = internaldiffuseheights->NumAltitudes();

	inboundunitspheres->reserve(numheights);
	ok =       (inboundunitspheres->capacity() >= numheights);
	ok = ok && MakeIncomingAzimuthAnglesV2( &azimuthgrid );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_GridSpecificationsLegacy_V2::MakeIncomingZenithGrid, Error allocating memory for Incoming unit sphere grid or azimuth grid");
	}
	else
	{
		for (hidx = 0; hidx < numheights; ++hidx )
		{
			entry   = new SKTRAN_UnitSphereLatLonGrid;
			zengrid = new  SKTRAN_GridDefDiffuseIncomingZenith_V21;
			ok1 = (entry != NULL) && ( zengrid != NULL);
			if (ok1)
			{
				entry->AddRef();
				zengrid->AddRef();

				if (m_useUserDefinedincomingzenith)	ok = ok && SetUserDefinedIncomingZenithGrid( zengrid, (hidx == 0) );
				else
				{
					if (m_shiftedHorizon)
						ok = ok && SetIncomingZenithGridShiftedHorizon( (hidx == 0), internaldiffuseheights->At(hidx), internaldiffuseheights, zengrid, m_rayregionmanager, coordinatesystem.get());
					else
						ok = ok && SetIncomingZenithGrid( (hidx == 0), internaldiffuseheights->At(hidx), internaldiffuseheights, zengrid, m_rayregionmanager, coordinatesystem.get());
				}

				ok = ok && entry->DefineGrid( zengrid, azimuthgrid );
				inboundunitspheres->push_back(entry );

				zengrid->Release();
			}
			else
			{
				nxLog::Record(NXLOG_WARNING, "SKTRAN_GridSpecificationsLegacy_V2::MakeIncomingZenithGrid, Error allocating memory for one element of incoming zenith grid");
				ok = false;
			}
		}
	}
	if (azimuthgrid != NULL) azimuthgrid->Release();
	
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_GridSpecificationsLegacy_V2::MakeIncomingAzimuthGridV2		2008-1-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::MakeIncomingAzimuthAnglesV2( SKTRAN_GridDefDiffuseIncomingAzimuth_V21** userincomingazimuthangle ) const
{
	bool										ok;
	size_t										aziidx;
	size_t										npts;
	SKTRAN_GridDefDiffuseIncomingAzimuth_V21*	incomingazimuthangle;

	NXTRACE_ONCEONLY(firsttime,("**** 2008-12-23 ***** SKTRAN_GridSpecificationsLegacy_V2::MakeIncomingAzimuthAngles, The incoming azimuths should be centered on the cell rather than at the left edge. We should changethe current code.\n"));
	incomingazimuthangle = new SKTRAN_GridDefDiffuseIncomingAzimuth_V21;
	ok      = (incomingazimuthangle != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_ERROR, "SKTRAN_SpecsUser_Diffuse_Legacy::MakeIncomingAzimuthAngles, Error allocating memory for azimuth grid");
	}
	else
	{
		incomingazimuthangle->AddRef();
		npts = m_incomingazimuth.size();
		ok   = incomingazimuthangle->AllocateGridArray(npts);
		if (!ok)
		{
			nxLog::Record(NXLOG_ERROR, "SKTRAN_SpecsUser_Diffuse_Legacy::MakeIncomingAzimuthAngles, Error allocating memory for azimuth buffer");
		}
		else
		{
			for (aziidx = 0; aziidx < npts; ++aziidx)
			{
				incomingazimuthangle->AtVar(aziidx) = m_incomingazimuth[aziidx];
			}
		}
	}
	*userincomingazimuthangle = incomingazimuthangle;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_GridSpecificationsLegacy_V2::ConfigureScatteringAngleGrid		2008-1-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::MakeScatterAngleGrid( SKTRAN_GridDefScatterAngle_V21** userscatteranglegrid ) const
{
	bool							ok;
	SKTRAN_GridDefScatterAngle_V21*	gridScatterAngle;

	gridScatterAngle = new SKTRAN_GridDefScatterAngle_V21;				// create the blank ray tracing  grid
	ok = (gridScatterAngle != NULL);									// make sure it worked
	if (!ok)															// if it did not
	{																	// then 
		nxLog::Record(NXLOG_WARNING,"SKTRAN_SpecsUser_Diffuse_Legacy::MakeScatterAngleGrid, Error allocating memory for scatter angle grid");
	}
	else
	{
		gridScatterAngle->AddRef();
		ok = gridScatterAngle->Configure( m_scatteringresolution_degrees );
		
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_SpecsUser_Diffuse_Legacy::MakeScatterAngleGrid, Error configuring the scatter angle arary");
		}
	}
	*userscatteranglegrid = gridScatterAngle;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::MakeSolarTransmissionSZAGrid		2008-1-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::MakeSolarTransmissionSZAGrid( SKTRAN_GridDefCosSZA_V21** userSZAgrid, SKTRAN_GridDefSLON_V21** userSLONgrid  ) const
{
	double						szaDelta;
	size_t						szaNum;
	size_t						idx;
	bool						ok;
	double						mincossza;
	double						maxcossza;
	double						maxsza;
	double						minsza;
	double						sza;
	std::vector<double>			szatransmission;						//!< The grid of cos(sza) used to determine the placing of the SZA transmission grid.
	SKTRAN_GridDefCosSZA_V21*	solartransmissionSzaGrid = NULL;
	SKTRAN_GridDefSLON_V21*		solartransmissionSlonGrid = NULL;


	ok = m_rayregionmanager->GetSZA( &sza, &minsza, &maxsza );
	if (m_useUserDefinedSZATransmissiongrid)
	{
		szatransmission = m_szatransmission;
	}
	else
	{

		NXTRACE_ONCEONLY(firsttime,("**** 2010-06-23 ***** SKTRAN_SpecsUser_Diffuse_Legacy::MakeSolarTransmissionSZAGrid, Need to review the Solar transmission SZA grid so its wide enough but not too wide\n"));
		if( m_configuredForMC ){
			// configure for Monte Carlo multiscatter (needs solar transmissions for entire illuminated hemisphere)
			minsza = 0.0;
			maxsza = 101.2;	// This is the largest sza where the observer is not in shadow for r_{earth}=6378, r_{atmomax}=6500
		} else{
			minsza = minsza - 10.0;								// we need to add room for all of the diffuse rays that zoom out
			minsza = std::min( minsza, sza-20.0);				// from each of the line of sight rays
			maxsza = maxsza + 10.0;
			maxsza = std::max( maxsza, sza + 20.0);

			if (minsza < 0) minsza = 0.0;
			if (maxsza > 180.0) maxsza = 180.0;
		}

		mincossza = nxmath::cosd( maxsza); 
		maxcossza = nxmath::cosd( minsza);
		szaNum  = (size_t)ceil((maxsza - minsza)/m_degreesPerSolarTransSza)+ 1;
		if ((szaNum & 0x01)  == 0) szaNum++;											// If SzaNum is even, increment to make it odd
		if (szaNum > 1)
		{
			szaDelta  = (maxcossza - mincossza)/(szaNum -1);
		}
		else
		{
			szaDelta  = 0.0;
			mincossza = 0.5*(mincossza + maxcossza);
		}

		szatransmission.reserve( szaNum );
		for (idx = 0; idx < szaNum; idx++)
		{
			szatransmission.push_back( mincossza + idx*szaDelta);
		}
	}

	szaNum = szatransmission.size();
	ok = (szaNum > 0 );
	if (ok)
	{
		solartransmissionSzaGrid  = new SKTRAN_GridDefCosSZA_V21;
		solartransmissionSlonGrid = new SKTRAN_GridDefSLON_V21;
		ok = (solartransmissionSzaGrid != NULL) && (solartransmissionSlonGrid != NULL);
		if (!ok)
		{
			nxLog::Record( NXLOG_WARNING, "SKTRAN_SpecsUser_Diffuse_Legacy::MakeSolarTransmissionSZAGrid, error allocating SKTRAN_GridDefSolarTransCosSZA_V2 object");
		}
		else
		{
			solartransmissionSzaGrid->AddRef();
			solartransmissionSlonGrid->AddRef();
			ok =       solartransmissionSzaGrid->AllocateGridArray( szaNum );
			ok = ok && solartransmissionSlonGrid->AllocateGridArray( szaNum );
			if (!ok)
			{
				nxLog::Record( NXLOG_WARNING, "SKTRAN_SpecsUser_Diffuse_Legacy::MakeSolarTransmissionSZAGrid, error allocaing %Iu points for SZA grid for solar geometry cache", (size_t)szaNum );
			}
			else
			{
				for (idx = 0; idx < szaNum; idx++)
				{
					solartransmissionSzaGrid->AtVar(idx) =  szatransmission[idx];
					solartransmissionSlonGrid->AtVar(idx) =  0.0;
				}
			}
		}
	}
	*userSZAgrid  = solartransmissionSzaGrid;
	*userSLONgrid = solartransmissionSlonGrid;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::MakeOutboundUnitSphere		2010-5-19*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::MakeOutboundUnitSphere( SKTRAN_UnitSphere_V2** userunitsphere ) const
{
	SKTRAN_UnitSphere_V2*					unitsphere;
	bool									ok;

	if ( m_unitsphere != NULL)
	{
		unitsphere = m_unitsphere;
	}
	else
	{
		unitsphere  = new SKTRAN_UnitSphereME(169);		// The unit sphere used for outgoing radiances
	}

	ok = (unitsphere != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,  "SKTRAN_SpecsUser_Diffuse_Legacy::MakeOutboundUnitSphere, Error allocating memory for outbound unit sphere object");
	}
	else
	{
		unitsphere->AddRef();
//		NXTRACE_ONCEONLY(Firsttime,("SKTRAN_SpecsUser_Diffuse_Legacy::MakeOutboundUnitSphere, We should delay making CreatLookuptable until needed as it is slow\n"));
//		ok = unitsphere->CreateLookupTable();  **** THIS IS NOW DONE IN **** SKTRAN_SpecsInternal_Diffuse_Legacy::FirstTimeInitializeDiffuseObjects
	}
	*userunitsphere = unitsphere;
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::MakeMaxDiffuseAltitude		2010-5-19*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::MakeMaxDiffuseAltitude( double* maxJAltitude ) const
{
	*maxJAltitude = m_maxdiffusealtitudealongincomingrays;
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::SetUserDefinedIncomingZenithGrid		2008-3-12*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::SetUserDefinedIncomingZenithGrid( SKTRAN_GridDefDiffuseIncomingZenith_V21* zenithgrid, bool isgroundpoint) const 
{
	bool							ok;
	size_t							numpoints;
	size_t								idx;
	std::vector<double>::const_iterator iter;

	if (isgroundpoint)
	{
		iter       = std::lower_bound( m_userdefinedincomingzenith.begin(), m_userdefinedincomingzenith.end(), 90.0 + 0.00001 );
		numpoints  = (iter - m_userdefinedincomingzenith.begin());
	}
	else
	{
		numpoints = m_userdefinedincomingzenith.size();
	}
	ok  = (numpoints > 0);
	ok  = ok && zenithgrid->AllocateGridArray( numpoints );
	if (ok)
	{
		for (idx = 0; idx < numpoints; idx++)
		{
			zenithgrid->AtVar(idx) = m_userdefinedincomingzenith[idx];
		}
		zenithgrid->SetIsGroundPoint(isgroundpoint);

	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::IsGroundPoint		2010-5-19*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::IsGroundPoint( double h ) const
{
	return ((h-1e-6) < m_diffusealtitudes[0]);
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureIncomingZenithGrid		2008-1-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::SetIncomingZenithGrid( bool											isgroundpoint,
															 double											altitude,
															 const SKTRAN_GridDefDiffuseHeights_V21*		internaldiffuseheights,
															       SKTRAN_GridDefDiffuseIncomingZenith_V21* zenithgrid,
															 const SKTRAN_RayTracingRegionManager*			rayregionmanager,
															 const SKTRAN_CoordinateTransform_V2*			coordinatesystem ) const

{
	double		zen_eighty;
	double		atmosDelta;
	double		horizonDelta;
	double		groundDelta;
	double		horizon;
	size_t		idx;
	size_t		numpoints;
	size_t		zenithCtr;
	size_t		horizonCtr;
	size_t		groundCtr;
	bool		ok;
	size_t		numatmos;
	size_t		numhorizon;
	size_t		numground;
//	bool		isgroundpoint;
	std::vector<double>	zenithbuffer;
	double				earthradius = coordinatesystem->AltitudeToRadius(0.0); // 6360000.0;
	double				groundaltitude;

	NXASSERT((earthradius >6300000.0) && (earthradius < 6400000.0));

	groundaltitude   = internaldiffuseheights->GroundAltitude();
	numatmos         = m_atmosRes;				// Number of zenith points in upward regions above 80 degrees zenith angle
	numhorizon       = m_horizonRes;			// Number of zenith points between 80 degrees zenith and the horizon
	numground        = m_groundRes;				// Number of zenith points below the horizon.

//		isgroundpoint    = IsGroundPoint( altitude );							// See is this is a ground point
	earthradius      = coordinatesystem->AltitudeToRadius( 0.0 );
	zenithgrid->SetIsGroundPoint(isgroundpoint);							//
	zen_eighty       = 80.0;												// begin horizon zeniths at 80 degrees 

	if (isgroundpoint )														// If this is for a point that hits the ground
	{																		// then only get half the number of points
		numhorizon   = numhorizon/2;										// For some reason halves the number of points for ground points
		numground    = 0;													// There are no points below the ground
		horizon      = 90.0;
		atmosDelta   = zen_eighty/numatmos;									// Atmosphere region between zenith angle 0 and 80 degrees.
		horizonDelta = ( horizon  - zen_eighty )/numhorizon;				// Horizon is SZA between 80 and the horizon.
		groundDelta  = 0.0;
	}
	else
	{
		horizon      = 90.0 + 180.0*acos( (groundaltitude + earthradius)/(altitude + earthradius) )/nxmath::Pi;	// Find the horizon.
		atmosDelta   = zen_eighty/numatmos;									// Atmosphere region between zenith angle 0 and 80 degrees.
		horizonDelta = ( horizon	- zen_eighty )/numhorizon;				// Horizon is SZA between 80 and the horizon.
		groundDelta  = ( 180.0		- horizon    )/numground;				// Ground is SZA between horizon and 180 degrees.
	}

	numpoints   = numatmos + numhorizon + numground + 1;
	zenithbuffer.assign(numpoints, 0.0);
	ok          = zenithgrid->AllocateGridArray( numpoints-1 );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_GridSpecificationsLegacy_V2::ConfigureIncomingZenithGrid, Error allocating memory for %Iu elements", (size_t)numpoints );
	}
	else
	{
		idx = 0;

		for( zenithCtr = 0; zenithCtr < numatmos; ++zenithCtr )
		{
			zenithbuffer[idx++] = zenithCtr * atmosDelta;
		}

		for( horizonCtr = 0; horizonCtr < numhorizon+1; ++horizonCtr )
		{
			zenithbuffer[idx++]  = zen_eighty  + horizonCtr * horizonDelta;
		}

		if (numground > 0)
		{
			for( groundCtr = numground; groundCtr > 0; --groundCtr )
			{
				zenithbuffer[idx++] = 180 - (groundCtr-1)* groundDelta;
			}
		}
		NXASSERT(( idx == numpoints ));
		for (idx=0; idx < numpoints-1; ++idx )
		{
			zenithgrid->AtVar(idx) = 0.5*(zenithbuffer[idx] + zenithbuffer[idx+1]);
		}
	}
	
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::ConfigureIncomingZenithGridShiftedHorizon		2013-11-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::SetIncomingZenithGridShiftedHorizon(	bool											isgroundpoint,
																			double											altitude,
																			const SKTRAN_GridDefDiffuseHeights_V21*			internaldiffuseheights,
																			SKTRAN_GridDefDiffuseIncomingZenith_V21*		zenithgrid,
																			const SKTRAN_RayTracingRegionManager*			rayregionmanager,
																			const SKTRAN_CoordinateTransform_V2*			coordinatesystem ) const

{
	double		atmosDelta;
	double		horizonDelta;
	double		groundDelta;
	double		horizon;
	double		startHorizon;
	double		endHorizon;
	double		horizonSize = 20;
	size_t		idx;
	size_t		numpoints;
	size_t		zenithCtr;
	size_t		horizonCtr;
	size_t		groundCtr;
	bool		ok;
	size_t		numatmos;
	size_t		numhorizon;
	size_t		numground;
//	bool		isgroundpoint;
	std::vector<double>	zenithbuffer;
	double				earthradius = 6360000.0;
	double				groundaltitude;

	
	groundaltitude   = internaldiffuseheights->GroundAltitude();
	numatmos         = m_atmosRes;				// Number of zenith points in upward regions above 80 degrees zenith angle
	numhorizon       = m_horizonRes;			// Number of zenith points between 80 degrees zenith and the horizon
	numground        = m_groundRes;				// Number of zenith points below the horizon.

//		isgroundpoint    = IsGroundPoint( altitude );							// See is this is a ground point
	earthradius      = coordinatesystem->AltitudeToRadius( 0.0 );
	zenithgrid->SetIsGroundPoint(isgroundpoint);							//

	if (isgroundpoint )														// If this is for a point that hits the ground
	{																		// then only get half the number of points
		numhorizon   = numhorizon/2;										// For some reason halves the number of points for ground points
		numground    = 0;													// There are no points below the ground
		horizon      = 90.0;
		startHorizon = horizon-horizonSize;
		endHorizon   = horizon+horizonSize;
		atmosDelta   = startHorizon/numatmos;									// Atmosphere region between zenith angle 0 and 80 degrees.
		horizonDelta = ( horizon  - startHorizon )/numhorizon;				// Horizon is SZA between 80 and the horizon.
		groundDelta  = 0.0;
	}
	else
	{
		horizon      = 90.0 + 180.0*acos( (groundaltitude + earthradius)/(altitude + earthradius) )/nxmath::Pi;	// Find the horizon.
		startHorizon = horizon-horizonSize;
		endHorizon   = horizon+horizonSize;
		atmosDelta   = startHorizon/numatmos;									// Atmosphere region between zenith angle 0 and 80 degrees.
		horizonDelta = ( endHorizon	- startHorizon )/numhorizon;				// Horizon is SZA horizonSize degrees around the horizon.
		groundDelta  = ( 180.0		- endHorizon    )/numground;				// Ground is SZA between horizon and 180 degrees.
	}

	numpoints   = numatmos + numhorizon + numground + 1;
	zenithbuffer.assign(numpoints, 0.0);
	ok          = zenithgrid->AllocateGridArray( numpoints-1 );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_GridSpecificationsLegacy_V2::ConfigureIncomingZenithGrid, Error allocating memory for %Iu elements", (size_t)numpoints );
	}
	else
	{
		idx = 0;

		for( zenithCtr = 0; zenithCtr < numatmos; ++zenithCtr )
		{
			zenithbuffer[idx++] = zenithCtr * atmosDelta;
		}

		for( horizonCtr = 0; horizonCtr < numhorizon+1; ++horizonCtr )
		{
			zenithbuffer[idx++]  = startHorizon  + horizonCtr * horizonDelta;
		}

		if (numground > 0)
		{
			for( groundCtr = numground; groundCtr > 0; --groundCtr )
			{
				zenithbuffer[idx++] = 180 - (groundCtr-1)* groundDelta;
			}
		}
		NXASSERT(( idx == numpoints ));
		for (idx=0; idx < numpoints-1; ++idx )
		{
			zenithgrid->AtVar(idx) = 0.5*(zenithbuffer[idx] + zenithbuffer[idx+1]);
		}
	}
	
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::::MakLOSInternalScatterFactory		2010-6-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::MakLOSInternalScatterFactory( SKTRAN_TableRayLOSFactory_Legacy** lossinglescatterfactory ) const
{
	SKTRAN_TableRayLOSFactory_Legacy*	factory;
	bool								ok;

	if (m_use_losinternalsinglescattertables)
	{
		factory = new SKTRAN_TableRayLOSFactory_Legacy;
		ok = (factory != NULL);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_SpecsUser_Diffuse_Legacy::::MakLOSInternalScatterFactory, Error allocating memory for LOS single scatter factory");
		}
		else
		{
			factory->AddRef();
		}
	}
	else
	{
		factory = NULL;
		ok      = true;
	}
	*lossinglescatterfactory = factory;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::IsProperlyDefined		2010-5-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::IsProperlyDefined() const
{
	NXTRACE_ONCEONLY(firsttime,("**** SKTRAN_SpecsUser_Diffuse_Legacy::IsProperlyDefined, needs a proper implementation\n"));
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy::CreateInternalSpecs		2010-5-19*/
/** This function is the instant when we transfer the users diffuse settings
 *	into the engines internal diffuse settings. After this there is no reliance
 *	upon the users diffuse settings. This will happen as part of the call to
 *	engine.ConfigureModel 
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Diffuse_Legacy::CreateInternalSpecs(       SKTRAN_SpecsInternal_Diffuse_V21**						engineinternalspecs,
														   std::shared_ptr< const SKTRAN_CoordinateTransform_V2>&       coordinatesystem
														 ) const
{

	SKTRAN_SpecsInternal_Diffuse_Legacy*				internalspecs;
	bool													ok;

	internalspecs = new SKTRAN_SpecsInternal_Diffuse_Legacy;
	ok = (internalspecs != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_SpecsUser_Diffuse_Legacy::CreateInternalSpecs, Error allocating memory");
	}
	else
	{
		internalspecs->AddRef();
		ok =       MakeDiffuseRadiiGrid				( &internalspecs->m_diffuseradii         );
		ok = ok && MakeDiffuseSzaGrid				( &internalspecs->m_diffuseSzaGrid, &internalspecs->m_diffuseSlonGrid );
//		ok = ok && MakeIncomingZenithAngles			( &internalspecs->m_incomingzenithangle, coordinatesystem );
//		ok = ok && MakeIncomingAzimuthAngles		( &internalspecs->m_incomingazimuthangle );
		ok = ok && MakeIncomingUnitSpheres			( internalspecs->m_diffuseradii, &internalspecs->m_incomingunitspheres, coordinatesystem );
		ok = ok && MakeScatterAngleGrid				( &internalspecs->m_gridScatterAngle     );
		ok = ok && MakeSolarTransmissionSZAGrid		( &internalspecs->m_solartransmissionSzaGrid, &internalspecs->m_solartransmissionSlonGrid);
		ok = ok && MakeOutboundUnitSphere			( &internalspecs->m_unitsphere );
		ok = ok && MakeMaxDiffuseAltitude			( &internalspecs->m_maxdiffusealtitudealongincomingrays );
		ok = ok && MakLOSInternalScatterFactory		( &internalspecs->m_lossinglescatterfactory );

		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_SpecsUser_Diffuse_Legacy::CreateInternalSpecs, Error making new specifications. This is not going to be good");
		}
	}
	*engineinternalspecs = internalspecs;
	return ok;
}
