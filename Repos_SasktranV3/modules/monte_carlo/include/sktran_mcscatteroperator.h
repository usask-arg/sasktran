#include "sktran_montecarlo_internals.h"

class SKTRAN_MCScatterOperator_Base
{

	protected:
		std::vector< const SKTRAN_Source_Term* >        m_sources;
		SKTRAN_HPFOSet*                                 m_hpfos;
		const SKTRAN_TableOpticalProperties_Base*       m_optprops;  
		SKTRAN_TableOpticalProperties_MCBase const*     m_mcoptprops;

	protected:
		void ChangePhotonBasis_atmoScatter (const double& cosTheta, SKTRAN_MCPhoton_Base* photon, const SKTRAN_RNG& rng, int order=2) const;
		void ChangePhotonBasis_groundScatter(const HELIODETIC_UNITVECTOR& groundScatterVector, const SKTRAN_RNG& rng, SKTRAN_MCPhoton_Base* photon, int order=2) const;
		
    protected:
        virtual void RotatePolarizedPhoton(SKTRAN_MCPhoton_Base* photon, double cosPhi, double sinPhi, int order=2) const = 0;  

	private:
		void ReleaseResources( );

public:
			     SKTRAN_MCScatterOperator_Base(  );
			     SKTRAN_MCScatterOperator_Base( const SKTRAN_MCScatterOperator_Base& other );
	virtual     ~SKTRAN_MCScatterOperator_Base(  );

	        bool SetCoordinateSystem  ( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>&   coords );
	virtual	bool SetOpticalProperties ( const SKTRAN_TableOpticalProperties_Base* optProp );
			
	        void ClearSourceTerms        ( );
		    void AddSourceTerm           ( const SKTRAN_Source_Term* s);
	virtual bool SetHPFlipOpSymmetry     ( SKTRAN_HPFOSet* hpfos );
			SKTRAN_HPFOSet * const GetHPFlipOp  ( ) { return m_hpfos;}
	virtual bool AcceptScatterPoint      ( const SKTRAN_MCPhoton_Base* mcphoton );
	virtual bool RandomScatter           ( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* photon, const SKTRAN_RNG& rng, int order ) const = 0;
	virtual bool CollectSourceTerms      ( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* scatterRay, const SKTRAN_CoordinateTransform_V2& coords, size_t threadid ) const = 0;
	        bool DeclareNewRay           ( );
			bool InputScatterOrder       ( int order );
	const std::vector< const SKTRAN_Source_Term* >& SourceTermVector ( ) { return m_sources;}

	virtual SKTRAN_MCScatterOperator_Base* ProduceCopy ( ) const = 0;
};

class SKTRAN_MCScatterOperatorContainer{
	private:
		std::shared_ptr<SKTRAN_MCScatterOperator_Base> m_op;

	public:
		SKTRAN_MCScatterOperatorContainer ( );
		SKTRAN_MCScatterOperatorContainer ( std::shared_ptr<SKTRAN_MCScatterOperator_Base>& op );
		SKTRAN_MCScatterOperatorContainer ( const SKTRAN_MCScatterOperatorContainer& other );
	   ~SKTRAN_MCScatterOperatorContainer ( );

             SKTRAN_MCScatterOperator_Base* Op( )       { return m_op.get();}
       const SKTRAN_MCScatterOperator_Base* Op( ) const { return m_op.get();}
};

class SKTRAN_MCScatterOperator_Scalar : public SKTRAN_MCScatterOperator_Base
{
	private:
		void         ReleaseResources        (  );
		virtual bool randomGroundScatter     ( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* photon, const SKTRAN_RNG& rng ) const;
		virtual bool randomAtmoScatter       ( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* photon, const SKTRAN_RNG& rng, int order=0 ) const;

    protected:
		virtual void RotatePolarizedPhoton(SKTRAN_MCPhoton_Base* photon, double cosPhi, double sinPhi, int order=2) const override {return;}

	public:
		SKTRAN_MCScatterOperator_Scalar      (  );
		SKTRAN_MCScatterOperator_Scalar      ( const SKTRAN_MCScatterOperator_Scalar& other );
		virtual ~SKTRAN_MCScatterOperator_Scalar     (  );

		virtual bool RandomScatter           ( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* photon, const SKTRAN_RNG& rng, int order ) const override;
		virtual bool CollectSourceTerms      ( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* scatterRay, const SKTRAN_CoordinateTransform_V2& coords, size_t threadid ) const override;

		virtual SKTRAN_MCScatterOperator_Scalar* ProduceCopy ( ) const override { return new SKTRAN_MCScatterOperator_Scalar(*this);}
};

class SKTRAN_MCScatterOperator_Polarized : public SKTRAN_MCScatterOperator_Base
{

	protected:
		void         ReleaseResources        (  );
     	bool         randomGroundScatter     ( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* photon, const SKTRAN_RNG& rng, int order=2 ) const;
		bool         randomAtmoScatter       ( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* photon, const SKTRAN_RNG& rng, int order=2 ) const;

	protected:
		virtual void RotatePolarizedPhoton(SKTRAN_MCPhoton_Base* photon, double cosPhi, double sinPhi, int order=2) const override;

	public:
		SKTRAN_MCScatterOperator_Polarized   (  );
		SKTRAN_MCScatterOperator_Polarized   ( const SKTRAN_MCScatterOperator_Polarized& other );
		virtual ~SKTRAN_MCScatterOperator_Polarized  (  );

		virtual bool RandomScatter           ( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* photon, const SKTRAN_RNG& rng, int order ) const override;
		virtual bool CollectSourceTerms      ( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* scatterRay, const SKTRAN_CoordinateTransform_V2& coords, size_t threadid ) const override;

		virtual SKTRAN_MCScatterOperator_Polarized* ProduceCopy ( ) const override { return new SKTRAN_MCScatterOperator_Polarized(*this);}

		virtual bool AcceptScatterPoint      ( const SKTRAN_MCPhoton_Base* mcphoton ) override;

		virtual bool SetHPFlipOpSymmetry     ( SKTRAN_HPFOSet* hpfos ) override;

};


class SKTRAN_MCScatterOperator_PseudoPolarized : public SKTRAN_MCScatterOperator_Polarized
{
	private:
		mutable SKTRAN_Stokes_NC m_foovec;

	private:
     	bool         randomGroundScatter     ( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* photon, const SKTRAN_RNG& rng, int order=2 ) const;
		bool         randomAtmoScatter       ( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* photon, const SKTRAN_RNG& rng, int order=2 ) const;

	protected:
		virtual void RotatePolarizedPhoton(SKTRAN_MCPhoton_Base* photon, double cosPhi, double sinPhi, int order=2) const override;

	public:
		virtual bool RandomScatter           ( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* photon, const SKTRAN_RNG& rng, int order ) const;
		virtual bool CollectSourceTerms      ( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* scatterRay, const SKTRAN_CoordinateTransform_V2& coords, size_t threadid ) const;

		virtual SKTRAN_MCScatterOperator_PseudoPolarized* ProduceCopy ( ) const override { return new SKTRAN_MCScatterOperator_PseudoPolarized(*this);}


};

class SKTRAN_MCScatterOperator_ScalarInelastic : public SKTRAN_MCScatterOperator_Scalar
{
	private:
		double m_lastwavelength;
		double m_minwavelength;
		double m_maxwavelength;
		void         ReleaseResources();
		virtual bool randomAtmoScatterInelastic(const HELIODETIC_POINT& scatterPoint, const double randNum, const SKTRAN_RNG& rng, SKTRAN_MCPhoton_Base* photon, int order = 0) const;

	public:
		SKTRAN_MCScatterOperator_ScalarInelastic();
		SKTRAN_MCScatterOperator_ScalarInelastic(const SKTRAN_MCScatterOperator_ScalarInelastic& other);
		virtual ~SKTRAN_MCScatterOperator_ScalarInelastic();

		virtual	bool SetOpticalProperties(const SKTRAN_TableOpticalProperties_Base* optProp) override;

		virtual bool RandomScatter(const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* photon, const SKTRAN_RNG& rng, int order) const override;

		virtual SKTRAN_MCScatterOperator_ScalarInelastic* ProduceCopy() const override { return new SKTRAN_MCScatterOperator_ScalarInelastic(*this); }

		virtual bool AcceptScatterPoint(const SKTRAN_MCPhoton_Base* mcphoton) override;

};