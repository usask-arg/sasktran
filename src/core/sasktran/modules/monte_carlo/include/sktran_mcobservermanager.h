#include "sktran_montecarlo_internals.h" 

// Manage observer LOS
// This should be deprecated -- observer should be described outside of engine classes
class SKTRAN_MCObserverManager
{
	private:
		SKTRAN_LineOfSightArray_V21 m_los;
		std::vector<double> m_pdfs;
		std::vector<double> m_cdfs;

		size_t Sub2Ind ( size_t startIdx, size_t endIdx) const {return startIdx*m_los.NumRays() + endIdx;}


	public:
		SKTRAN_MCObserverManager  ( );
		~SKTRAN_MCObserverManager ( );

		bool   SetRays                       ( const SKTRAN_LineOfSightArray_V21& rays );
		void   CacheProbabilityDistributions ( );
		size_t DrawRandomLOS                 ( double rand ) const;
		size_t DrawProposalLOSMutation       ( size_t currentIndex, double rand ) const;
		double TransitionProb                ( size_t startIdx, size_t endIdx ) const;
};
