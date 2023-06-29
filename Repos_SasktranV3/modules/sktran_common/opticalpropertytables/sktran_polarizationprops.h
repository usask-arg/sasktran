


/*-----------------------------------------------------------------------------
 *					SKTRAN_PolarizationProperties_NoPolarization		 2016- 7- 12*/
/** @ingroup polarize
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_PolarizationProperties_NoPolarization : public SKTRAN_PolarizationProperties_Base
{
    private:
    std::vector< double > m_phaseFunction;

	public:
                       SKTRAN_PolarizationProperties_NoPolarization ( );
		virtual       ~SKTRAN_PolarizationProperties_NoPolarization ( );
		virtual bool   Allocate                         ( size_t numcachepoints ) override;
		virtual bool   StorePolarizationPropsCM2        ( size_t index, const skRTPhaseMatrix& pmatrix, SKTRAN_AtmosphericOpticalState_V21& state ) override;
        virtual bool   GetPhaseFunctionCM2              ( const SKTRAN_GridIndex gridindex[], const double gridweight[], size_t numels, double& phaseFunction ) const override;
		virtual bool   GetScatteringMatrixCM2           ( const SKTRAN_GridIndex gridindex[], const double gridweight[], size_t numels, SKTRAN_ScatMat_MIMSNC& matrix ) const override;
		virtual bool   GetResultOfUnpolarizedScatterCM2 ( const SKTRAN_GridIndex gridindex[], const double gridweight[], size_t numels, SKTRAN_Stokes_NC& stokesvec ) const override;

        virtual double PhaseMatrixAccess                  ( size_t idx ) const override { return m_phaseFunction[idx]; }
        virtual size_t NumElements                        ( )            const override { return m_phaseFunction.size(); }
		virtual bool NeedsFullMueller					  ( )			 const override { return false; }

};


/*-----------------------------------------------------------------------------
 *					SKTRAN_PolarizationProperties_Polarized		 2016- 7- 12*/
/** @ingroup polarize
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_PolarizationProperties_Polarized : public SKTRAN_PolarizationProperties_Base
{
	public:
		typedef SKTRAN_ScatMat_MIMSNC ScatMatType; // Should be templated

	private:
		std::vector<ScatMatType>   m_muellerMatrices; 

	public:
                       SKTRAN_PolarizationProperties_Polarized ( );
		virtual       ~SKTRAN_PolarizationProperties_Polarized ( );
		virtual bool   Allocate	                               ( size_t numcachepoints ) override;
		virtual bool   StorePolarizationPropsCM2               ( size_t index, const skRTPhaseMatrix& pmatrix, SKTRAN_AtmosphericOpticalState_V21& state ) override;
        virtual bool   GetPhaseFunctionCM2                     ( const SKTRAN_GridIndex gridindex[], const double gridweight[], size_t numels, double& phaseFunction ) const override;
		virtual bool   GetScatteringMatrixCM2                  ( const SKTRAN_GridIndex gridindex[], const double gridweight[], size_t numels, SKTRAN_ScatMat_MIMSNC& matrix ) const override;
		virtual bool   GetResultOfUnpolarizedScatterCM2        ( const SKTRAN_GridIndex gridindex[], const double gridweight[], size_t numels, SKTRAN_Stokes_NC& stokesvec ) const override;

        virtual double PhaseMatrixAccess                  ( size_t idx ) const override { return m_muellerMatrices[idx].p11();}
        virtual size_t NumElements                        ( )            const override { return m_muellerMatrices.size();}
		virtual bool NeedsFullMueller                     ( )			 const override { return true; }
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_PolarizationProperties_Polarized_Eddington		 2016- 7- 12*/
/** SKTRAN_PolarizationProperties_Polarized_Eddington
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_PolarizationProperties_Polarized_Eddington : public SKTRAN_PolarizationProperties_Polarized
{
	private:
		std::vector<double>				m_eddingtonExtinction;

	public:
                       SKTRAN_PolarizationProperties_Polarized_Eddington      ( );
		virtual       ~SKTRAN_PolarizationProperties_Polarized_Eddington      ( );
		virtual bool   Allocate	                 ( size_t numcachepoints ) override;
		virtual bool   StorePolarizationPropsCM2 ( size_t index, const skRTPhaseMatrix& pmatrix, SKTRAN_AtmosphericOpticalState_V21& state ) override;

};

