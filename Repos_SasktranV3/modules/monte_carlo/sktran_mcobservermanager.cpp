#include "include/sktran_montecarlo_internals.h"


SKTRAN_MCObserverManager::SKTRAN_MCObserverManager( )
{
}

SKTRAN_MCObserverManager::~SKTRAN_MCObserverManager( )
{
}

bool SKTRAN_MCObserverManager::SetRays( const SKTRAN_LineOfSightArray_V21& rays )
{
	return m_los.DeepCopy( rays );
}


void SKTRAN_MCObserverManager::CacheProbabilityDistributions( )
{
	//bool ok = true;
	//double norm, thisp;

	//m_cdfs.resize(m_los.NumRays() * m_los.NumRays() );
	//m_pdfs.resize(m_los.NumRays() * m_los.NumRays() );

	//// For now just do this terrible thing:
	//for(int fromidx=0; fromidx<m_los.NumRays(); ++fromidx){
	//	norm = 0.0;
	//	for(int endidx=0; endidx<m_los.NumRays(); ++endidx){
	//		thisp = m_los.Entry(fromidx)->Look() & m_los.Entry(endidx)->Look();
	//		m_pdfs[Sub2Ind(fromidx,endidx)] = thisp;
	//		norm += thisp;
	//	}
	//	for(int endidx=0; endidx<m_los.NumRays(); ++endidx){
	//		m_pdfs[Sub2Ind(fromidx,endidx)] /= norm;
	//	}
	//	m_cdfs[Sub2Ind(fromidx,0)] = m_pdfs[Sub2Ind(fromidx,0)];
	//	for(int endidx=1; endidx<m_los.NumRays()-1; ++endidx){
	//		m_cdfs[Sub2Ind(fromidx,endidx)] = m_cdfs[Sub2Ind(fromidx,endidx-1)] + m_pdfs[Sub2Ind(fromidx,endidx)];
	//	}
	//	m_cdfs[Sub2Ind(fromidx,m_los.NumRays()-1)] = 1.0;
	//}

	return;
}

double SKTRAN_MCObserverManager::TransitionProb( size_t startIdx, size_t endIdx ) const
{
	//return m_pdfs[Sub2Ind(startIdx,endIdx)];
	return  1.0 / m_los.NumRays();
}


size_t SKTRAN_MCObserverManager::DrawRandomLOS ( double rand ) const
{
	double resultd = m_los.NumRays() * rand;
	size_t result  = min(size_t(std::floor( resultd )), m_los.NumRays()-1);
	return result;
}

size_t SKTRAN_MCObserverManager::DrawProposalLOSMutation( size_t currentIndex, double rand ) const
{
	//auto it = std::lower_bound(m_cdfs.begin()+Sub2Ind(currentIndex,0), m_cdfs.begin()+Sub2Ind(currentIndex,0)+m_los.NumRays(), rand);
	//return std::distance(m_cdfs.begin()+Sub2Ind(currentIndex,0), it);
	return DrawRandomLOS( rand );
}