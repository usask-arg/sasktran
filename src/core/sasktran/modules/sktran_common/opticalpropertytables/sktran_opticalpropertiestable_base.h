//#pragma once
//
//#include <modules/sktran_common/include/sktran_common_internals.h>
//#include "sktran_specifications_base.h"


/*-----------------------------------------------------------------------------
*					SKTRAN_TableOpticalProperties_Base		2007-11-20*/
/** @ingroup optprop
*	A class used to cache the optical properties of the atmosphere for
*	a range of altitudes and scattering angles.
**/
/*---------------------------------------------------------------------------*/
class SKTRAN_TableOpticalProperties_Base : public nxUnknown
{
	protected:
		std::shared_ptr< const SKTRAN_CoordinateTransform_V2>	m_coordinates;
        std::unique_ptr< SKTRAN_PolarizationProperties_Base >	m_scatprops;

	public:
		enum OPTIONSENUM {  OPTION_IS_HORIZONTALLY_UNIFORM = 1};				//!< True if the optical properties are only a function of altitude (and time which is always ignored)


	public:
		 											SKTRAN_TableOpticalProperties_Base		();
		virtual 								   ~SKTRAN_TableOpticalProperties_Base		();
		bool										SetCoords								( std::shared_ptr< const SKTRAN_CoordinateTransform_V2>& coordinates );
		void                                        SetPolarizationProperties               ( std::unique_ptr< SKTRAN_PolarizationProperties_Base >& polprops );
		const SKTRAN_CoordinateTransform_V2*		CoordinatesPtr							() const					{ return m_coordinates.get();}
        std::shared_ptr< const SKTRAN_CoordinateTransform_V2> Coordinates                   ( ) const { return m_coordinates; }
		static bool									GroundBRDFAngles						( const HELIODETIC_POINT& point, const HELIODETIC_UNITVECTOR& incomingray, const HELIODETIC_UNITVECTOR& outgoingray, double* mu_in, double* mu_out, double* cosdphi );
		bool										Get_AlbedoForDeprecatedLegacyCode		( const HELIODETIC_POINT& point, double* brdf ) const;

	public:
		virtual bool								IsOptionTrue							( SKTRAN_TableOpticalProperties_Base::OPTIONSENUM options) const = 0;
		virtual bool								ConfigureOptical						( double wavelen, SKTRAN_AtmosphericOpticalState_V21& opticalstate ) = 0;

		virtual double								TotalExtinctionPerCM					( const HELIODETIC_POINT& point ) const = 0;
		virtual bool								GetBRDF									( const HELIODETIC_POINT& point, double mu_in, double mu_out, double cosdphi, double* brdf ) const = 0;
		virtual bool								GetBRDFGeodetic							( const GEODETIC_INSTANT& point, double mu_in, double mu_out, double cosdphi, double* brdf ) const = 0;
		virtual bool								GetScatteringCoefficientCM2				( const HELIODETIC_POINT& point, double cosangle,  SKTRAN_PhaseMatrixScalar* phasematrix ) const = 0;
		virtual double								ScatteringExtinctionPerCM				( const HELIODETIC_POINT& point ) const = 0;
		virtual bool								GetEffectiveExtinctionPerCMWithHeight1	( const HELIODETIC_POINT& startpoint, const HELIODETIC_POINT& endpoint, double* sigma0, double* sigma1 ) const = 0;

		virtual bool								GetEffectiveExtinctionPerCMWithHeight1	( const SKTRAN_RayStorage_Base* r, size_t startPtIndex, double* sigma0, double* sigma1 ) const { return false; }
	public:
		virtual bool								GetLinearExtinctionPerCMVector			( const std::vector< HELIODETIC_POINT>& quadpoints, std::vector<double>& sigmak, std::vector<double>& sigmaf, size_t numpoints ) const { nxLog::Record(NXLOG_WARNING, "SKTRAN_TableOpticalProperties_Base::GetLinearExtinctionPerCMVector is not implemented"); return false; };
		virtual bool                                GetScatteringMatrixCM2                  ( const HELIODETIC_POINT& point, double cosAngle, SKTRAN_ScatMat_MIMSNC&  pmatrix   ) const {return false;} // Have to give dummy implementation for legacy
		virtual bool                                GetResultOfUnpolarizedScatterCM2        ( const HELIODETIC_POINT& point, double cosAngle, SKTRAN_Stokes_NC& stokesvec ) const {return false;}

        virtual bool                                CreateInterpolationForPoint             ( const HELIODETIC_POINT& point, SKTRAN_TableOpticalProperties_PointCache& interpolator ) const = 0; // Returns a function that gives the scattering matrix as a function of scattering angle 

	public:
		virtual double								TotalExtinctionPerCM					(double wavelength, const HELIODETIC_POINT& point) const = 0;
		virtual bool								GetBRDF									(double wavelength, const HELIODETIC_POINT& point, double mu_in, double mu_out, double cosdphi, double* brdf) const = 0;
		virtual bool								GetBRDFGeodetic							(double wavelength, const GEODETIC_INSTANT& point, double mu_in, double mu_out, double cosdphi, double* brdf) const = 0;
		virtual bool								GetEffectiveExtinctionPerCMWithHeight1	(double wavelength, const HELIODETIC_POINT& startpoint, const HELIODETIC_POINT& endpoint, double* sigma0, double* sigma1) const = 0;
		virtual bool								GetEffectiveExtinctionPerCMWithHeight1	(double wavelength, const SKTRAN_RayStorage_Base* r, size_t startPtIndex, double* sigma0, double* sigma1) const { return false;}
		virtual bool								GetScatteringCoefficientCM2				(double wavelength, const HELIODETIC_POINT& point, double cosangle, SKTRAN_PhaseMatrixScalar* scatcoeff) const = 0;
		virtual double								ScatteringExtinctionPerCM				(double wavelength, const HELIODETIC_POINT& point) const = 0;
		virtual bool								GetScatteringMatrixCM2					(double wavelength, const HELIODETIC_POINT& point, double cosAngle, SKTRAN_ScatMat_MIMSNC&  pmatrix) const { return false;}
		virtual bool								GetResultOfUnpolarizedScatterCM2		(double wavelength, const HELIODETIC_POINT& point, double cosAngle, SKTRAN_Stokes_NC& stokesvec) const { return false;}
		virtual bool								CreateInterpolationForPoint				(double wavelength, const HELIODETIC_POINT& point, SKTRAN_TableOpticalProperties_PointCache& interpolator) const = 0;

};


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableEmission_Base		 2015- 3- 2*/
/** @ingroup emission
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_TableEmission_Base
{
	protected:
		std::shared_ptr< const SKTRAN_CoordinateTransform_V2>	m_coordinates;

	public:
													SKTRAN_TableEmission_Base		(){};
		virtual 								   ~SKTRAN_TableEmission_Base		(){};
		bool										SetCoords								(std::shared_ptr< const SKTRAN_CoordinateTransform_V2> coordinates){ m_coordinates = coordinates; return true;}
		const SKTRAN_CoordinateTransform_V2*		CoordinatesPtr							() const					{ return m_coordinates.get();}

	public:
		virtual double								GetIsotropicRadianceInAtmosphere		( const HELIODETIC_POINT& point ) const = 0;
		virtual double								GetIsotropicGroundRadiance				( const HELIODETIC_POINT& point ) const = 0;
		virtual bool								IsEmpty									() const = 0;
		virtual void								SetEmpty								( bool value) = 0;

};
