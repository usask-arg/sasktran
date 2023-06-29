
class SKTRAN_Engine_Base;

enum ENUM_SKTRAN_DIAGNOSTIC_STAGE { SKTRAN_DIAGNOSE_CONFIGUREOPTICAL,			//!< diagnostic invoked after tables are initialized with incoling solat irradiance
									SKTRAN_DIAGNOSE_FIRSTORDERINCOMING,			//!< diagnostic invoked after table sare initialized with first order incoming irradiance
									SKTRAN_DIAGNOSE_SCATTERINCOMING,			//!< diagnostic invoked after incoming radys of table have been scattered to output
									SKTRAN_DIAGNOSE_CALCULATEINCOMING,
									SKTRAN_DIAGNOSE_FINISH
								  };

typedef  bool (* SKTRAN_DiagnosticFunction)( enum ENUM_SKTRAN_DIAGNOSTIC_STAGE stage, double wavelengthnm, size_t orderofscatter, bool processingok, SKTRAN_Engine_Base* engine );



/*-----------------------------------------------------------------------------
 *					SKTRAN_DiagnosticInterface		2009-1-15*/
/** @ingroup diagnostics
 *	A class used to invoke user defined diagnostic functions at various 
 *	stages of processing. Really useful for extracting tables
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_DiagnosticInterface
{
	private:

		std::list<SKTRAN_DiagnosticFunction>	m_diagnosticfunctions;		//!< A list of diagnostic functions called at the end of processing each wavelength
		SKTRAN_Engine_Base*					m_engineinterface;
		double									m_wavelennm;

	public:
												SKTRAN_DiagnosticInterface	();
											   ~SKTRAN_DiagnosticInterface	();
		void									Initialize					(SKTRAN_Engine_Base*	engine); 
		void									SetWavelength				( double m_wavelennm );
		bool									AddDiagnosticFunction		(SKTRAN_DiagnosticFunction f);
		bool									ClearDiagnosticFunctions	();
		bool									Diagnose					( ENUM_SKTRAN_DIAGNOSTIC_STAGE stage, size_t orderofscatter, bool processingok );
};

