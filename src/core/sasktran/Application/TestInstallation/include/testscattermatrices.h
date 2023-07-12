
class TestScatterMatrices{

	private:
		const double m_tolerance;

	private:
		template<typename SCATMAT1, typename SCATMAT2>
		bool Compare_Impl( const SCATMAT1& p, const SCATMAT2& s) const;
		bool Compare ( const skRTPhaseMatrix& p,               const SKTRAN_ScatMat_MIMSNC& s)         const { return Compare_Impl(p,s);}
		bool Compare ( const skRTPhaseMatrix& p,               const SKTRAN_PhaseMat_MIMSNC& s) const { return Compare_Impl(p,s);}
		bool Compare ( const skRTPhaseMatrix& p,               const SKTRAN_ScatMat_Rot& s)            const { return Compare_Impl(p,s);}
		bool Compare ( const SKTRAN_ScatMat_MIMSNC& p,         const SKTRAN_ScatMat_MIMSNC& s)         const { return Compare_Impl(p,s);}
		bool Compare ( const SKTRAN_PhaseMat_MIMSNC& p, const SKTRAN_PhaseMat_MIMSNC& s) const { return Compare_Impl(p,s);}
		bool Compare ( const SKTRAN_ScatMat_Rot& p,            const SKTRAN_ScatMat_Rot& s)            const { return Compare_Impl(p,s);}
		bool Compare ( const SKTRAN_PhaseMat_MIMSNC& p, const SKTRAN_ScatMat_MIMSNC& s)         const { return Compare_Impl(p,s);}

		bool Compare ( double d1, double d2 ) const;
		bool Compare ( const SKTRAN_Stokes_NC& v1, const SKTRAN_Stokes_NC& v2) const;

		bool TestRONC           ( ) const;
		bool TestRot            ( ) const;
		bool TestRONC_Composite ( ) const;

	public:
		TestScatterMatrices () : m_tolerance(1e-10){};
		bool RunTest() const;
};