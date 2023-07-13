#pragma once

class EtaCalculator_Base
{
	protected:
		class etaStore
		{
			public:
				double cosEta_refToScatt, sinEta_refToScatt;
				double cosEta_scattToRef, sinEta_scattToRef;
		
										etaStore(){;} // Don't initialize doubles to zero
		};
		const double m_dtol = 1.0; //-1e-10;
		const double m_denomtol = 1e-10;

	public: 
		virtual			   ~EtaCalculator_Base  ( ) {;}
		virtual void		SetPoint            ( const SKTRAN_HR_Diffuse_Point& point ) = 0;
		virtual void		UpdateOutgoingIndex ( const SKTRAN_HR_Diffuse_Point& point, size_t outidx ) = 0;
        virtual void        SetOutgoingDirection (const HELIODETIC_UNITVECTOR& out) {};
		virtual void		ScattToRef          ( SKTRAN_Stokes_NC& stokes ) const = 0;
		virtual void		RefToScatt          ( SKTRAN_Stokes_NC& stokes ) const = 0;
		virtual void		CalculateEtas       ( const SKTRAN_HR_Diffuse_Point& point, size_t inidx  ) = 0;
		virtual void		ScattMatToPhaseMat  ( const SKTRAN_ScatMat_MIMSNC& s, SKTRAN_PhaseMat_MIMSNC& p ) const = 0; // Could template over the first argument to allow other scat mat symmetries 
};


class EtaCalculator_DoRotation : public EtaCalculator_Base
{
	private:
	struct inPlaneVec { double x, y; };

	HELIODETIC_UNITVECTOR      m_outray;

	double                     m_invSinOutToSun;
	double                     m_nol2;
	double                     m_sinOutToSun;
	inPlaneVec                 m_outgoing_inPlane; 
	size_t                     m_outidx;
	std::vector< inPlaneVec >  m_incoming_inPlane; 
	std::vector<double>        m_sinInToSun;
	std::vector<double>        m_invSinInToSun;

	/*double                     m_sinOutToSun;*/
	/*HELIODETIC_UNITVECTOR      m_qPtVecs[3];*/
	//double                     m_invSinOutToSun;
	//double                     m_nol2;
	//inPlaneVec                 m_outgoing_inPlane; 
	//size_t                     m_outidx;

	double m_crts, m_srts, m_cstr, m_sstr;

	HELIODETIC_UNITVECTOR      m_qPtVecs[3];

	private:
		inline size_t sub2ind ( size_t inidx, size_t outidx, size_t numin ) const { return (outidx*numin) + inidx; }

	public: 
		virtual void SetPoint            ( const SKTRAN_HR_Diffuse_Point& point ) override;
		virtual void UpdateOutgoingIndex ( const SKTRAN_HR_Diffuse_Point& point, size_t outidx ) override;
        virtual void SetOutgoingDirection (const HELIODETIC_UNITVECTOR& out) override;

        virtual void ScattToRef          ( SKTRAN_Stokes_NC& stokes ) const override;
		virtual void RefToScatt          ( SKTRAN_Stokes_NC& stokes ) const override;
		virtual void CalculateEtas       ( const SKTRAN_HR_Diffuse_Point& point, size_t inidx ) override;
		virtual void ScattMatToPhaseMat  ( const SKTRAN_ScatMat_MIMSNC& s, SKTRAN_PhaseMat_MIMSNC& p ) const override;
};

class EtaCalculator_NoRotation : public EtaCalculator_Base
{
	public: 
		virtual void SetPoint            ( const SKTRAN_HR_Diffuse_Point& point ) override;
		virtual void UpdateOutgoingIndex ( const SKTRAN_HR_Diffuse_Point& point, size_t outidx ) override;
		virtual void ScattToRef          ( SKTRAN_Stokes_NC& stokes ) const override;
		virtual void RefToScatt          ( SKTRAN_Stokes_NC& stokes ) const override;
		virtual void CalculateEtas       ( const SKTRAN_HR_Diffuse_Point& point, size_t inidx ) override;
		virtual void ScattMatToPhaseMat  ( const SKTRAN_ScatMat_MIMSNC& s, SKTRAN_PhaseMat_MIMSNC& p ) const override;
};    


