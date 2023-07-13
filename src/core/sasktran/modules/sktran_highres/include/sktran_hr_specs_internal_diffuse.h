//#include "sktran_hr_internals.h"

class SKTRAN_HR_Diffuse_Point;
class SKTRAN_HR_Diffuse_Second_Order_Source;

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_Diffuse		2014-10-30*/
/** Factory class for creating the diffuse points used for the diffuse table.
 *  Note that the diffuse table is created in SKTRAN_HR_Specs_Internal_Core
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_HR_Specs_Internal_Diffuse
{
	private:
		std::shared_ptr<const SKTRAN_CoordinateTransform_V2>	m_coords;
		size_t													m_numabove;			//!< Number of points to use above the horizon
		size_t													m_numhoriz;			//!< Number of points to use in the horizon
		size_t													m_numbelow;			//!< Number of points to use below the horizon
		size_t													m_numazi;			//!< Number of azimuth resolutions to use
		//SKTRAN_UnitSphereME*									m_outgoingsphere;
		SKTRAN_HR_OutgoingSphereObject_Base*					m_outgoingsphereobj;
		size_t													m_numoutgoing;
		bool													m_uselegacyv21sphere;
		double													m_horizonsize;
		SKTRAN_HR_DiffuseIncomingSphereType						m_incomingspheretype;
        SKTRAN_HR_DiffuseMatrixStorageMethod                    m_diffuseMatrixStorageMethod;

	private:
		bool							MakeDefaultIncomingAziGrid( SKTRAN_GridDefDiffuseIncomingAzimuth_V21& azigrid, const double& altitude );
		bool							MakeDefaultIncomingZenGrid( SKTRAN_GridDefDiffuseIncomingZenith_V21&  zengrid, const double& altitude, const bool& isgroundpoint );
		bool							MakeSasktran21IncomingZenGrid( SKTRAN_GridDefDiffuseIncomingZenith_V21& zengrid, const double& altitude, const bool& isgroundpoint );
		bool							MakeSasktran21HorizonShiftZenGrid( SKTRAN_GridDefDiffuseIncomingZenith_V21& zengrid, const double& altitude, const bool& isgroundpoint );
		bool							ReleaseResources();
		bool							ConfigureDefaults();
	public:
										SKTRAN_HR_Specs_Internal_Diffuse();
		virtual						   ~SKTRAN_HR_Specs_Internal_Diffuse();
		bool							ConfigureIncomingUnitSphere		( const SKTRAN_UnitSphere_V2** unitsphere, const double& altitude, bool isground );
		bool							ConfigureOutgoingUnitSphere		();
		bool							MakeDiffusePoint				( SKTRAN_HR_Diffuse_Point& diffusepoint, const HELIODETIC_POINT& location, bool isground );
		bool							MakeSecondOrderSource			( SKTRAN_HR_Diffuse_Second_Order_Source** source );
		bool							Initialize						( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords );
		bool							Configure						( const SKTRAN_HR_Specs_User& specs );
        SKTRAN_HR_DiffuseMatrixStorageMethod GetDiffuseMatrixStorageMethod ( ) { return m_diffuseMatrixStorageMethod;}
};
