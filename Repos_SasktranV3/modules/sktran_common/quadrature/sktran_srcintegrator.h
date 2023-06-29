

/*-----------------------------------------------------------------------------
 *					SKTRAN_OpticalPropertiesIntegrator_Base		2014-2-7*/
/** @ingroup srcintegrate
 *	Integrate the contribution of source terms along a ray whose transmission has already been calculated. 
 *  Assume that sources for each integration segment are well-represented by the source at the center of 
 *  those segments. 
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_SourceTermIntegrator_Order0 : public SKTRAN_SourceTermIntegrator_Base
{
	template<typename radtype> class Integration_Impl{ // Why is this not just a template function? 
		public:
			bool IntegrateSourceTerm	( const SKTRAN_RayOptical_Base* ray, radtype& radiance, const std::vector<SKTRAN_Source_Term*>& sources, const SKTRAN_TableOpticalProperties_Base* optprops ) const;
	};

	public:
		virtual bool  IntegrateSourceTerm ( const SKTRAN_RayOptical_Base* ray, double& radiance,           const std::vector<SKTRAN_Source_Term*>& sources ) const override;
		virtual bool  IntegrateSourceTerm ( const SKTRAN_RayOptical_Base* ray, SKTRAN_Stokes_NC& radiance, const std::vector<SKTRAN_Source_Term*>& sources ) const override;
	
};
									
/*-----------------------------------------------------------------------------
 *					SKTRAN_OpticalPropertiesIntegrator_Base		2014-2-7*/
/** @ingroup srcintegrate
 *	Integrate the contribution of source terms along a ray whose transmission has already been calculated. 
 *  Assume that sources for each integration segment are well-represented by the source at the center of 
 *  those segments. 
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_SourceTermIntegrator_Order2 : public SKTRAN_SourceTermIntegrator_Base
{
	private:
		template<typename radtype> bool     GetQuadraticCoeff         ( const radtype& sourcestart, const radtype& sourcemid, const radtype& sourceend, double ds, radtype& b, radtype& c ) const;
		template<typename radtype> radtype  GetRadianceContrib        ( const radtype& sourcestart, const radtype& sourcemid, const radtype& sourceend, double ds, double sigmastart, double sigmaend ) const;
		template<typename radtype> bool     IntegrateSourceTerm_Impl  ( const SKTRAN_RayOptical_Base* ray, radtype& radiance, const std::vector<SKTRAN_Source_Term*>& sources ) const;

	public:
		virtual bool									IntegrateSourceTerm							( const SKTRAN_RayOptical_Base* ray, double& radiance,           const std::vector<SKTRAN_Source_Term*>& sources ) const override;
		virtual bool                                    IntegrateSourceTerm                         ( const SKTRAN_RayOptical_Base* ray, SKTRAN_Stokes_NC& radiance, const std::vector<SKTRAN_Source_Term*>& sources ) const override;
};
				