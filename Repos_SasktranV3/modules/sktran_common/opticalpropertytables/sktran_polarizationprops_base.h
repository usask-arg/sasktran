/*-----------------------------------------------------------------------------
 *					SKTRAN_PolarizationProperties_Base		2014-05-27*/
/** @ingroup polarize
 *	Cache the normalized Mueller matrices on the same grid as is done for optical
 *  property tables. 
 **/
/*---------------------------------------------------------------------------*/
class SKTRAN_PolarizationProperties_Base 
{
	public:
                       SKTRAN_PolarizationProperties_Base ( ) {}
		virtual       ~SKTRAN_PolarizationProperties_Base ( ) {}
		virtual bool   Allocate	                          ( size_t numcachepoints ) = 0;
		virtual bool   StorePolarizationPropsCM2          ( size_t index, const skRTPhaseMatrix& pmatrix, SKTRAN_AtmosphericOpticalState_V21& state ) = 0;
        virtual bool   GetPhaseFunctionCM2                ( const SKTRAN_GridIndex gridindex[], const double gridweight[], size_t numels, double& phaseFunction ) const = 0;
		virtual bool   GetScatteringMatrixCM2             ( const SKTRAN_GridIndex gridindex[], const double gridweight[], size_t numels, SKTRAN_ScatMat_MIMSNC& matrix ) const = 0;
		virtual bool   GetResultOfUnpolarizedScatterCM2   ( const SKTRAN_GridIndex gridindex[], const double gridweight[], size_t numels, SKTRAN_Stokes_NC& stokesvec ) const = 0;

        virtual double PhaseMatrixAccess                  ( size_t idx ) const = 0;
        virtual size_t NumElements                        ( )            const = 0;
		virtual bool NeedsFullMueller                     ( )			 const = 0;
};

