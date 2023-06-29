//#include "sktran_hr_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_User_Diffuse		2013-08-23*/
/** @ingroup HR_SPECS 
 *  Class which provides user options for the management of diffuse profiles.
 *  This currently includes
 *    - Number of diffuse profiles
 *    - Their placement
 *    - Resolution on incoming/outgoing spheres
 * 
 *  Currently the default options are,
 *    - 1 Diffuse profile placed at the tangent point
 *       - 100 Points, spaced at 500m, 1500m, etc,
 *    - 6 Zenith angles before horizon area, 8 in horizon area, and 10 after horizon
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_HR_Specs_User_Diffuse
{
	private:
		size_t						m_numprofiles;						//!< Number of diffuse profiles in the look plane
		size_t						m_numbeforehoriz;					//!< Number of zenith angles before the horizon area
		size_t						m_numhoriz;							//!< Number of zenith angles in the horizon area
		size_t						m_numafterhoriz;					//!< Number of zenith angles after the horizon area
		size_t						m_numazi;							//!< Number of azimuth angles
		size_t						m_numoutgoing;						//!< Number of rays on the outgoing unit sphere
		size_t						m_numofflook;						//!< number of planes off the look direction to play diffuse profiles
		double						m_angleofflook;						//!< Maximum angle off of the look plane to place diffuse profiles
		double						m_heightres;						//!< Height resolution to use on diffuse profiles in m
		double						m_maxdiffuseheight;					//!< This option is not normally used as diffuse profiles should span the entire altitude range. It is used for backward compatibility with old SASKTRAN code.
//		double						m_surfaceheight;					//!< Minimum height for diffuse profiles in meters
		//bool						m_forcelinearplacement;				//!< if true then diffuse profiles are placed linearly along LOS instead of in SZA
		bool						m_forcev21incomingsphere;
		std::vector<double>			m_manualdiffuseheights;
		std::vector<size_t>			m_diagnosticscatorders;				//!< Orders of scatter for which diagnostic data are collected
		std::vector<size_t>			m_diagnosticdiffprofs;				//!< Diffuse profiles for which diagnostic data are collected
		std::vector<nxVector>		m_manualdiffuselocations;
		std::vector<double>			m_manualdiffuselatlons;
		std::vector<double>			m_manualdiffuseszas;
		std::vector<double>			m_manualdiffuselospositions;
		nxVector					m_diffuseplanereference;
		nxVector					m_diffuseplanenormal;
		std::vector<double>			m_diffuseplaneangles;
		double						m_horizonsize;

		SKTRAN_HR_DiffuseIncomingSphereType		m_incomingspheretype;
		SKTRAN_HR_DiffuseProfilePlacementType	m_diffuseplacement;
        SKTRAN_HR_DiffuseMatrixStorageMethod    m_diffuseMatrixStorageMethod;

	private:
		void						ConfigureDefaults();
	public:
		SKTRAN_HR_Specs_User_Diffuse() { ConfigureDefaults(); }
		~SKTRAN_HR_Specs_User_Diffuse() {};
		/**
		 *  Sets the total number of diffuse profiles in the look direction plane
		 *  Note the total number of diffuse profiles will be m_numprofiles * m_numofflook
		 **/
		bool						SetNumProfiles( size_t numprofiles ) { m_numprofiles = numprofiles; return true; }
		size_t						GetNumProfiles() const { return m_numprofiles; }
		void						SetNumBeforeHoriz( size_t num ) { m_numbeforehoriz = num; }
		size_t						GetNumBeforeHoriz() const { return m_numbeforehoriz; }
		void						SetNumHoriz( size_t num ) { m_numhoriz = num; }
		size_t						GetNumHoriz() const { return m_numhoriz; }
		void						SetNumAfterHoriz( size_t num ) { m_numafterhoriz = num; }
		size_t						GetNumAfterHoriz() const { return m_numafterhoriz; }
		void						SetNumAzi( size_t num ) { m_numazi = num; }
		size_t						GetNumAzi() const { return m_numazi; }
		void						SetNumOutgoing( size_t num ) { m_numoutgoing = num; }
		size_t						GetNumOutgoing() const { return m_numoutgoing; }
		void						SetNumoffLook( size_t num) { m_numofflook = num; }
		size_t						GetNumOffLook() const { return m_numofflook; }
		void						SetAngleOffLook( double angle ) { m_angleofflook = angle; }
		double						GetAngleOffLook() const { return m_angleofflook; }
		void						SetHeightRes( double heightres ) { m_heightres = heightres; }
		double						GetHeightRes() const { return m_heightres; }
		void						SetMaxDiffuseHeight( double maxheight ) { m_maxdiffuseheight = maxheight; }
//		void						SetSurfaceHeight( double groundheight){ m_surfaceheight = groundheight;}
//		double						GetSurfaceHeight() const { return m_surfaceheight;}
		double						GetMaxDiffuseHeight() const { return m_maxdiffuseheight; }
		void						SetForceLinearPlacement( bool force ) { if( force)
																			{
																				m_diffuseplacement = SKTRAN_HR_DiffuseProfilePlacement_LinearLOS;
																			}
																			else
																			{
																				m_diffuseplacement = SKTRAN_HR_DiffuseProfilePlacement_LinearSZAForceTP;
																			}}
		void						SetForceV21IncomingSphere( bool force ) { m_forcev21incomingsphere = force; }
		void						SetHorizonSize( double size ) { m_horizonsize = size; }
		void						SetManualDiffuseLocations( const std::vector<nxVector>& locations ) { m_manualdiffuselocations = locations; }
		void						SetManualDiffuseLatLons( const std::vector<double>& latlons ) { m_manualdiffuselatlons = latlons; }
		void						SetManualDiffuseSZAs( const std::vector<double>& szas ) { m_manualdiffuseszas = szas; }
		void						SetManualDiffuseLOSPositions( const std::vector<double>& positions ) { m_manualdiffuselospositions = positions; }
		void						SetManualDiffusePlaneReference( const nxVector& reference ) { m_diffuseplanereference = reference; }
		void						SetManualDiffusePlaneNormal( const nxVector& normal ) { m_diffuseplanenormal = normal; }
		void						SetManualDiffusePlaneAngles( const std::vector<double>& angles ) { m_diffuseplaneangles= angles; }
		void						SetManualDiffuseHeights( const std::vector<double>& heights ) { m_manualdiffuseheights = heights; }
		const std::vector<double>&	GetManualDiffuseHeights() const { return m_manualdiffuseheights; }
		void						SetDiagnosticScatterOrders( const std::vector<size_t>& orders ) { m_diagnosticscatorders = orders; }
		void						SetDiagnosticDiffuseProfiles( const std::vector<size_t>& profiles ) { m_diagnosticdiffprofs = profiles; }
		void						SetDiffusePlacementType( SKTRAN_HR_DiffuseProfilePlacementType type ) { m_diffuseplacement = type; }
        void                        SetScatterMatrixStorageMethod( SKTRAN_HR_DiffuseMatrixStorageMethod method ){ m_diffuseMatrixStorageMethod = method;}
        SKTRAN_HR_DiffuseMatrixStorageMethod GetScatteringMatrixStorageMethod ( ) const { return m_diffuseMatrixStorageMethod;}

		void						SetDiffuseIncomingType( SKTRAN_HR_DiffuseIncomingSphereType type ) { m_incomingspheretype = type; }
		friend class SKTRAN_HR_Specs_Internal_Diffuse;
		friend class SKTRAN_HR_Specs_Internal_Core;
};
