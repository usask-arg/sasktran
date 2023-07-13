//#pragma once

//#include "sktran_common_internals.h"




class SKTRAN_RayFactory_Base;
/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_Base		2013-05-24*/
/** @ingroup solartransmission
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_SolarTransmission_Base : public SKTRAN_Source_Term
{
	private:
        std::unique_ptr<SKTRAN_Sun_Base>                    m_sun;
		SKTRAN_OpticalPropertiesIntegrator_Base*			m_integrator;
		std::weak_ptr< const SKTRAN_RayFactory_Base>		m_rayfactory;

	protected:
		virtual bool										FillTable_ClassSpecific             ( )          = 0;

		bool												EstimateNormalizedPolarizationVector( const HELIODETIC_POINT&      pt,
																								  const HELIODETIC_UNITVECTOR& look,
																								  SKTRAN_Stokes_NC&            vec,
																								  const HELIODETIC_BASIS&      basis ) const;

		bool												EstimateNormalizedPolarizationVector( const double&				   wavelength,
																								  const HELIODETIC_POINT&      pt,
																								  const HELIODETIC_UNITVECTOR& look,
																								  SKTRAN_Stokes_NC&            vec,
																								  const HELIODETIC_BASIS&      basis ) const;

        virtual double                                      CosAngleToSource                    ( const HELIODETIC_UNITVECTOR& look, const HELIODETIC_POINT* location) const { return m_sun->CosAngleToSun(look);}
		

	private:
		void												ReleaseResources                    ( );	// These should all be private


	public:
															SKTRAN_SolarTransmission_Base		();
		virtual											   ~SKTRAN_SolarTransmission_Base		();
		const SKTRAN_RayFactory_Base*						RayFactory							() const { return m_rayfactory.lock().get();}
		const SKTRAN_OpticalPropertiesIntegrator_Base*		Integrator							() const { return m_integrator;}
		const SKTRAN_Sun_Base*								Sun									() const { return m_sun.get();}
        bool												FillTable							( );
        bool												SetSun                              ( std::unique_ptr<SKTRAN_Sun_Base>& sun );
        bool												ConfigureOptical					( std::weak_ptr< const SKTRAN_RayFactory_Base>	rayfactory, SKTRAN_OpticalPropertiesIntegrator_Base* integrator);

		// Need to bring overloaded methods into scope 
		using SKTRAN_Source_Term::SourceTermAtPoint;
		using SKTRAN_Source_Term::GroundSourceAtPoint;
		using SKTRAN_Source_Term::MonteCarlo_SingleScatteredRadianceAtPoint;
		using SKTRAN_Source_Term::MonteCarlo_GroundScatteredRadianceAtPoint;
	
	public:
        virtual bool										MakeThreadSafeFor            ( size_t numThreads ) = 0;
        virtual bool										SourceTermAtPoint			 ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC* source   ) const override;
		virtual bool										GroundSourceAtPoint			 ( const SKTRAN_SourceTermQueryObject_Base& outboundray, double* source ) const override;
		virtual bool										GroundSourceAtPoint			 ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC* source   ) const override;
		virtual bool										MonteCarlo_SingleScatteredRadianceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC& radiance ) const override;
		virtual bool										MonteCarlo_GroundScatteredRadianceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC& radiance ) const override;
		
		virtual bool										SourceTermAtPoint( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC* source) const override;
		virtual bool										GroundSourceAtPoint(const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& outboundray, double* source) const override;
		virtual bool										GroundSourceAtPoint(const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC* source) const override;
		virtual bool										MonteCarlo_SingleScatteredRadianceAtPoint ( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC& radiance) const override;
		virtual bool										MonteCarlo_GroundScatteredRadianceAtPoint ( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC& radiance) const override;

        // These should be removed from the interface eventually 
        virtual bool  TransmissionAtPoint            ( const HELIODETIC_POINT&   point, double&  transmission) const = 0;
        virtual bool  TransmissionAtVector           ( const HELIODETIC_VECTOR&  vec,   double&  transmission) const = 0;
		virtual bool  TransmissionAtPoint			 ( const double& wavelength, const HELIODETIC_POINT&   point, double&  transmission) const { return false; }
		virtual bool  TransmissionAtVector			 ( const double& wavelength, const HELIODETIC_VECTOR&  vec,   double&  transmission) const { return false; }

};

/*-----------------------------------------------------------------------------
*					SKTRAN_SolarTransmission_DoNothing		2016-10-21*/
/** \ingroup common
 *  A dummy solar transmission table that does nothing. Its purpose is to allow
 *  access to the base class' access methods during RT calculations where the 
 *  solar transmission is assumed to be zero, e.g. for photochemical emissions.
**/
/*---------------------------------------------------------------------------*/
class SKTRAN_SolarTransmission_DoNothing : public SKTRAN_SolarTransmission_Base
{
    protected:
        virtual bool  FillTable_ClassSpecific        ( ) { return true;}

    public:
                      SKTRAN_SolarTransmission_DoNothing  ( ) { ;}
                     ~SKTRAN_SolarTransmission_DoNothing  ( ) { ;}
        virtual bool  MakeThreadSafeFor                   ( size_t numThreads ){ return true;}
        virtual bool  SourceTermAtPoint                   ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_StokesScalar* source   ) const override { *source        = 0.0;  return true;}
		virtual bool  GroundSourceAtPoint                 ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_StokesScalar* source   ) const override { *source        = 0.0;  return true;}
		virtual bool  MonteCarlo_SingleScatteredRadianceAtPoint      ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_StokesScalar& radiance ) const override { radiance      = 0.0;  return true;}
		virtual bool  MonteCarlo_GroundScatteredRadianceAtPoint      ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_StokesScalar& radiance ) const override { radiance      = 0.0;  return true;}
        virtual bool  SourceTermAtPoint                              ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC*    source   ) const override { source  ->SetTo(0.0); return true;} 
        virtual bool  GroundSourceAtPoint                            ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC*    source   ) const override { source  ->SetTo(0.0); return true;} 
        virtual bool  MonteCarlo_SingleScatteredRadianceAtPoint      ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC&    radiance ) const override { radiance .SetTo(0.0); return true;}
        virtual bool  MonteCarlo_GroundScatteredRadianceAtPoint      ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC&    radiance ) const override { radiance .SetTo(0.0); return true;}

		virtual bool  SourceTermAtPoint								 ( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_StokesScalar* source   ) const override { *source       = 0.0;  return true;}
		virtual bool  GroundSourceAtPoint							 ( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_StokesScalar* source   ) const override { *source       = 0.0;  return true;}
		virtual bool  MonteCarlo_SingleScatteredRadianceAtPoint      ( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_StokesScalar& radiance ) const override { radiance      = 0.0;  return true;}
		virtual bool  MonteCarlo_GroundScatteredRadianceAtPoint      ( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_StokesScalar& radiance ) const override { radiance      = 0.0;  return true;}
        virtual bool  SourceTermAtPoint                              ( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC*    source   ) const override { source  ->SetTo(0.0); return true;}
        virtual bool  GroundSourceAtPoint                            ( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC*    source   ) const override { source  ->SetTo(0.0); return true;}
        virtual bool  MonteCarlo_SingleScatteredRadianceAtPoint      ( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC&    radiance ) const override { radiance .SetTo(0.0); return true;}
        virtual bool  MonteCarlo_GroundScatteredRadianceAtPoint      ( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC&    radiance ) const override { radiance .SetTo(0.0); return true;}

        // These should be removed from the interface eventually 
        virtual bool  TransmissionAtPoint                 ( const HELIODETIC_POINT&   point, double&  transmission) const override { transmission=0.0; return true;}
        virtual bool  TransmissionAtVector                ( const HELIODETIC_VECTOR&  vec,   double&  transmission) const override { transmission=0.0; return true;} 

};


