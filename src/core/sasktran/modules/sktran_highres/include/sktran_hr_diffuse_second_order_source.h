//#include "sktran_hr_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Second_Order_Source	2014-10-30*/
/**  Allows for analytical calculation of the second order source term.
 *   Whenever a source is request, a diffuse point is placed at the location
 *   of the request, and the second order source is calculated.
 *
 *   Used primarily only for testing purposes
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_HR_Diffuse_Second_Order_Source : public SKTRAN_Source_Term
{
	private:
		const SKTRAN_UnitSphere_V2*							m_incomingsphere;
		std::shared_ptr< const SKTRAN_RayFactory_Base>		m_rayfactory;
		const SKTRAN_OpticalPropertiesIntegrator_Base*		m_integrator;
		const SKTRAN_TableOpticalProperties_Base*			m_opticaltable;
		const SKTRAN_SolarTransmission_Base*				m_solartable;
        const SKTRAN_EmissionTable_Base*                    m_emissiontable;

	private:
		bool												RotateIncomingRay			( const nxVector& inray, const HELIODETIC_POINT& pt, HELIODETIC_UNITVECTOR& outlook ) const;
		bool												RadianceFromLook			( const HELIODETIC_UNITVECTOR& look, const HELIODETIC_POINT& pt, double& radiance ) const;
		const SKTRAN_RayFactory_Base*						RayFactoryPtr				() const { return m_rayfactory.get();}

	public:
															SKTRAN_HR_Diffuse_Second_Order_Source();
														   ~SKTRAN_HR_Diffuse_Second_Order_Source();

		virtual bool										SourceTermAtPoint					( const SKTRAN_SourceTermQueryObject_Base& qobj,
																								  double*							source ) const override;

		virtual bool										SourceTermAtPoint					( const SKTRAN_SourceTermQueryObject_Base& qobj,
																								  SKTRAN_Stokes_NC*					source ) const override;

		virtual bool										SetUnitSphere						( const SKTRAN_UnitSphere_V2& unitsphere );

		virtual bool										SetOptical							( std::shared_ptr< const SKTRAN_RayFactory_Base>	rayfactory,
																								  const SKTRAN_OpticalPropertiesIntegrator_Base&	integrator,
																								  const SKTRAN_TableOpticalProperties_Base&			opticaltable,
																								  const SKTRAN_SolarTransmission_Base&				solartable, 
                                                                                                  const SKTRAN_EmissionTable_Base&                  emissions );
		// These are unused by HR engine 
        
		virtual bool										GroundSourceAtPoint					( const SKTRAN_SourceTermQueryObject_Base& qobj,
																								  double*							source   ) const override { skRTStokesVector::SetToZero( *source ); return true;}

		virtual	bool										GroundSourceAtPoint					( const SKTRAN_SourceTermQueryObject_Base& qobj,
																								  SKTRAN_Stokes_NC*					source   ) const override { skRTStokesVector::SetToZero( *source ); return true;}
		

};
