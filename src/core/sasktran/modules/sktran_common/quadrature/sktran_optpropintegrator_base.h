//#pragma once

//#include <iostream>
//#include "sktran_common_internals.h"

class SKTRAN_SolarTransmission_Base;
class SKTRAN_Source_Term;

/*-----------------------------------------------------------------------------
 *					SKTRAN_OpticalPropertiesIntegrator_Base		2014-2-7*/
/** @ingroup odintegrate
 */
/*---------------------------------------------------------------------------*/

class SKTRAN_OpticalPropertiesIntegrator_Base : public nxUnknown
{
	protected:
		const SKTRAN_TableOpticalProperties_Base*		m_opticalprops;

	protected:
		virtual double									OpticalDepthOfCell							( const SKTRAN_RayOptical_Base* ray, size_t cellidx ) const = 0;

		// given a scatter point defined by start + scatterDistance * look, move the scatter point back slightly if it is right on the edge of the atmosphere
		virtual bool									RoundScatterPosition						( const HELIODETIC_POINT& start, const HELIODETIC_UNITVECTOR& look, const double& scatterDistance, const SKTRAN_CoordinateTransform_V2* coords, HELIODETIC_VECTOR & scatterVector, HELIODETIC_POINT& scatterPosition ) const;

	public:
														SKTRAN_OpticalPropertiesIntegrator_Base ( );
		virtual										   ~SKTRAN_OpticalPropertiesIntegrator_Base ( );
		virtual void									ReleaseResources                        ( );

		virtual bool									CalculateRayScalarTransmission					( SKTRAN_RayOptical_Base*			ray,
																									  double*							transmisison,
																									  bool								totaltransmissiononly,
																									  bool								usecachedtransmission
																									) const = 0;
		virtual bool                                    CalculateRayScalarTransmission_withMinContainer	( SKTRAN_RayOptical_Base* ray, double* transmission, bool totaltransmissiononly, bool usecachedtransmission  ) const = 0;
		virtual bool									CalculateRayScalarTransmissionVector				( SKTRAN_RayOptical_Base* ray, double* transmission, bool totaltransmissiononly, bool usecachedtransmission, std::vector<double>* sigmak = NULL, std::vector<double>* sigmaf = NULL ) const = 0;

		virtual bool									SingleScatterRadiance						( SKTRAN_RayOptical_Base* ray, double& radiance, const SKTRAN_Source_Term& solartable ) const = 0;

		// Find a random scatter position along a ray, and P( scatterPoint | thisRay ) for the domain of P being distance along the ray
		virtual bool									FindNextScatterPosition						( double rand, double resolution, SKTRAN_TableOpticalProperties_Base* cachedOpticalProperties, SKTRAN_RayOptical_Base const* originalRay, HELIODETIC_VECTOR& scatterVector, double& scatterProb, double& randomOpticalDepth, double userImposedDistance ) const = 0;
		virtual bool									CalculatePartialOpticalDepth				( SKTRAN_RayOptical_Base const* ray, size_t cellindex, double fraction, double& opticaldepth ) const = 0;

		bool											SetOpticalProps								( const SKTRAN_TableOpticalProperties_Base* optprop );
		const SKTRAN_TableOpticalProperties_Base*		GetOpticalProps								( ) const { return m_opticalprops; }

};

