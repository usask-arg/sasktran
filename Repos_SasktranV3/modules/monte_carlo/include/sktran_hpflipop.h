#include "sktran_montecarlo_internals.h" 

class SKTRAN_HVFlipOp
{
	private:
		double m_store[9]; 

	public:
		SKTRAN_HVFlipOp                   ( );
		SKTRAN_HVFlipOp                   ( const SKTRAN_HVFlipOp& other   );
		SKTRAN_HVFlipOp        Define    ( const SKTRAN_MCPhoton_Base* photon  ); // Define the Householder transform across the plane containing the observer and its scatter point
		HELIODETIC_UNITVECTOR  operator* ( const HELIODETIC_UNITVECTOR& v ) const;
		HELIODETIC_VECTOR      operator* ( const HELIODETIC_VECTOR& v     ) const;
		SKTRAN_MCBasis         operator* ( const SKTRAN_MCBasis& b        ) const;
		SKTRAN_HVFlipOp        LeftApply ( const SKTRAN_HVFlipOp& hpfo    ); // Left multiply myself by #hpfo
};

class SKTRAN_HPFOSet : public nxUnknown
{
	protected:
		std::vector<SKTRAN_HVFlipOp>   m_hvfos;
		std::shared_ptr<const SKTRAN_CoordinateTransform_V2> m_coords;
		int m_order;
		
	public:
		              SKTRAN_HPFOSet        ( );
		virtual      ~SKTRAN_HPFOSet        ( );
		        bool  SetCoords             ( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords );
		virtual bool  ProducePointCache     ( const HELIODETIC_POINT& hp, const HELIODETIC_UNITVECTOR& look, std::vector<HELIODETIC_POINT>& pc, std::vector< HELIODETIC_UNITVECTOR >& lc ) const = 0;
		virtual bool  ApplyFlipsToBasis     ( const SKTRAN_MCBasis& basis, std::vector<SKTRAN_MCBasis>& flippedBasis ) const = 0; 
		virtual bool  InputOrder            ( int order );
		int GetOrder ( ) const { return m_order;}
		virtual bool  InputSourcesVariance  ( double variance ) = 0;
		virtual bool  AddHPFO               ( const SKTRAN_HVFlipOp& hpfo ) = 0;
		virtual bool  AddHPFO               ( const SKTRAN_MCPhoton_Base* mcphoton ) = 0;
		virtual bool  DeclareNewRay         ( ); // Invalidate any existing flip operators
		virtual bool  CopyInto              ( SKTRAN_HPFOSet& dest ) const;
		virtual bool  SetMaxNumberFlipOrders ( int maxOrder ) = 0;

		virtual SKTRAN_HPFOSet* ProduceCopy ( ) const = 0;
};

class SKTRAN_HPFOSet_NoSymmetry : public SKTRAN_HPFOSet
{
	public:
		SKTRAN_HPFOSet_NoSymmetry           ( ){ }
		SKTRAN_HPFOSet_NoSymmetry           ( const SKTRAN_HPFOSet_NoSymmetry& other ) {other.CopyInto(*this);}
		virtual ~SKTRAN_HPFOSet_NoSymmetry  ( ){ }
		virtual bool ProducePointCache      ( const HELIODETIC_POINT& hp, const HELIODETIC_UNITVECTOR& look, std::vector<HELIODETIC_POINT>& pc, std::vector< HELIODETIC_UNITVECTOR >& lc ) const override;
		virtual bool ApplyFlipsToBasis      ( const SKTRAN_MCBasis& basis, std::vector<SKTRAN_MCBasis>& flippedBasis ) const override; 
		virtual bool InputSourcesVariance   ( double variance ) override;
		virtual bool AddHPFO                ( const SKTRAN_HVFlipOp& hpfo ) override;
		virtual bool AddHPFO                ( const SKTRAN_MCPhoton_Base* mcphoton ) override;
		
		virtual bool SetMaxNumberFlipOrders ( int maxOrder ) override;

		virtual SKTRAN_HPFOSet_NoSymmetry* ProduceCopy ( ) const { return new SKTRAN_HPFOSet_NoSymmetry(*this);}
};


class SKTRAN_HPFOSet_HorizSymmetry : public SKTRAN_HPFOSet
{
    private:
		int    m_order_stopPushing;
		int    m_order_stopPushing_max;
		int    m_order_clear;
		int    m_numRunsSinceCheck_stopPushing;
		int    m_numRunsSinceCheck_clear;
		int    m_numRunsBetweenChecks;
		double m_lastReportedVariance;
		double m_runningvariance_stopPushing;
		double m_runningvariance_clear;
		double m_highThreshold;
		double m_lowThreshold;
		std::vector<double> m_variance;

	private:
		bool AdjustLimits( );

	public:
		SKTRAN_HPFOSet_HorizSymmetry          ( );
		SKTRAN_HPFOSet_HorizSymmetry          ( const SKTRAN_HPFOSet_HorizSymmetry& other );
		virtual ~SKTRAN_HPFOSet_HorizSymmetry ( ){}
		virtual bool ProducePointCache        ( const HELIODETIC_POINT& hp, const HELIODETIC_UNITVECTOR& look, std::vector<HELIODETIC_POINT>& pc, std::vector< HELIODETIC_UNITVECTOR >& lc ) const override;
		virtual bool ApplyFlipsToBasis        ( const SKTRAN_MCBasis& basis, std::vector<SKTRAN_MCBasis>& flippedBasis ) const override; 
		virtual bool InputSourcesVariance     ( double variance ) override;
		virtual bool AddHPFO                  ( const SKTRAN_HVFlipOp& hpfo ) override;
		virtual bool AddHPFO                  ( const SKTRAN_MCPhoton_Base* mcphoton ) override;
		virtual bool DeclareNewRay            ( ) override; // Invalidate any existing flip operators

		virtual bool SetMaxNumberFlipOrders ( int maxOrder ) override;
		
		virtual SKTRAN_HPFOSet_HorizSymmetry* ProduceCopy ( ) const { return new SKTRAN_HPFOSet_HorizSymmetry(*this);}
};

