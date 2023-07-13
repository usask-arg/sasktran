
/*-----------------------------------------------------------------------------
 *					SKTRAN_SourceTermIntegrator_Base		 2014- 11- 6*/
/** @ingroup srcintegrate
 *	Base class for evaluating the integrals of source terms along a ray
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_SourceTermIntegrator_Base : public nxUnknown
{
	protected:
		mutable int  m_fidx;
		const SKTRAN_TableOpticalProperties_Base*		m_opticalprops;

	private:
		void                                            ReleaseResources( );

	public:
														SKTRAN_SourceTermIntegrator_Base            ( );
													   ~SKTRAN_SourceTermIntegrator_Base            ( );

		bool											SetOpticalProps								( const SKTRAN_TableOpticalProperties_Base* optprop );
		const SKTRAN_TableOpticalProperties_Base*		GetOpticalProps								( ) const { return m_opticalprops; }
		virtual bool									IntegrateSourceTerm							( const SKTRAN_RayOptical_Base* ray, double& radiance,           const std::vector<SKTRAN_Source_Term*>& sources ) const = 0;
		virtual bool                                    IntegrateSourceTerm                         ( const SKTRAN_RayOptical_Base* ray, SKTRAN_Stokes_NC& radiance, const std::vector<SKTRAN_Source_Term*>& sources ) const = 0;
};