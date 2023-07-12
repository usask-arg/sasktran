//#include "sktran_hr_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_Source_Term		2014-10-30*/
/** Creates a source term from a diffuse table
 *
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_HR_Diffuse_Source : public SKTRAN_Source_Term
{
	private:
		const SKTRAN_HR_Diffuse_Table_CPU*		m_diffusetable;
	public:
												SKTRAN_HR_Diffuse_Source() {};
		virtual								   ~SKTRAN_HR_Diffuse_Source() {};

		
		bool							SetDiffuseTable							( const SKTRAN_HR_Diffuse_Table_CPU& diffusetable ) { m_diffusetable = &diffusetable; return true; }
		virtual double					IntegrateScalarIncoming					( const SKTRAN_SourceTermQueryObject_Base& qobj,
																				  std::function<double(const nxVector&, const nxVector&)> func ) const 
																				{ return m_diffusetable->IntegrateScalarIncoming( qobj, func); }

    virtual SKTRAN_Stokes_NC					IntegrateVectorIncoming			( const SKTRAN_SourceTermQueryObject_Base& qobj,
                                                                                   std::function<SKTRAN_ScatMat_MIMSNC(const nxVector&, const nxVector&)> func ) const
                                                                                { return m_diffusetable->IntegrateVectorIncoming( qobj, func); }

		virtual bool					SourceTermAtPoint						( const SKTRAN_SourceTermQueryObject_Base& qobj, double* source             ) const override;
		virtual bool					SourceTermAtPoint						( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC* source   ) const override;
		virtual	bool					GroundSourceAtPoint						( const SKTRAN_SourceTermQueryObject_Base& qobj, double*           source   ) const override;
		virtual bool					GroundSourceAtPoint						( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC* source   ) const override; 
};

