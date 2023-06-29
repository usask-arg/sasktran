#include "include/sktran_montecarlo_internals.h"

SKTRAN_Sun_RandomDisc::SKTRAN_Sun_RandomDisc( )
{
	nxLog::Record(NXLOG_INFO,"SKTRAN_Sun_RandomDisc::Constructor, We need to check out that this class works properly after moving the thread stroage around. Might have to adjust how we initialize the random number generators");
	m_currentSun    = nullptr;
	m_rngs          = nullptr;
	m_sineApexAngle = 0.0;
	m_currentSun1.SetThreadStorageInitializer( std::bind( &SKTRAN_Sun_RandomDisc::InitializeThreadEntry, this, std::placeholders::_1) );
}

SKTRAN_Sun_RandomDisc::~SKTRAN_Sun_RandomDisc( )
{
	ReleaseResources();
}

bool SKTRAN_Sun_RandomDisc::InitializeSunVector( HELIODETIC_UNITVECTOR* threadentry) const
{
	threadentry->SetCoords(0.0, 0.0, 1.0);
	return true;
}

bool SKTRAN_Sun_RandomDisc::Initialize( double sineApexAngle, const std::vector<SKTRAN_RNG>& randGens, const size_t numThreads )
{
	bool ok = true;

	m_currentSun1.Clear();


	ok = ok && 0<=sineApexAngle;
	ok = ok && 0<numThreads;
	if(ok)
	{
		m_sineApexAngle = sineApexAngle;
		// ###########################################
		// UNSAFE -- fix memory leak before committing
		// ###########################################
		nxLog::Record( NXLOG_ERROR, "SKTRAN_Sun_RandomDisc::Initialize, This class was broken when threadid handling was changed." );
//		m_rngs          = randGens;			// ndl: Would need a rng container to do reference counting
		m_rngs = &randGens;
//		m_currentSun    = new HELIODETIC_UNITVECTOR[numThreads];
		m_currentSun = new HELIODETIC_UNITVECTOR[numThreads];
	}
//	ok = ok && NULL!=m_currentSun;
	if(ok)
	{
		for(size_t sunidx=0; sunidx<numThreads; sunidx++)
		{ 
			InitializeSunVector( &m_currentSun[sunidx] );
		}
	}
	else
	{
		nxLog::Record(NXLOG_ERROR, "SKTRAN_Sun_RandomDisc::Initialize, Could not initialize random sun.");
	}

	return ok;
}

void SKTRAN_Sun_RandomDisc::ReleaseResources( )
{
	// #################
	// Hack to make disc sun work 
	// ########################
	delete[] m_currentSun; m_currentSun = NULL;
	m_rngs = nullptr;				// Would be better if we did some reference counting, but just trust engine for now
	m_sineApexAngle = -1.0;
}

		
double SKTRAN_Sun_RandomDisc::GetSineApexAngle( ) const
{
	return m_sineApexAngle;
}

bool SKTRAN_Sun_RandomDisc::InitializeThreadEntry( SKTRAN_Sun_RandomDisc_ThreadEntry* threaddata ) const
{
	bool	ok = true;
	//ok = InitializeSunVector( &threaddata->m_sun);
	return ok;
}


void SKTRAN_Sun_RandomDisc::UpdateSun() const
{
	double r, theta, norm;

	//m_currentSun1.LookupUpThreadData( &threaddata );
	//r     = std::sqrt( threaddata->m_rng()) * m_sineApexAngle;	// Select a random circle on the disc
	//theta = threaddata->m_rng() * nxmath::TWOPI;		// Select a random point on that circle
	//norm  = std::sqrt(1.0 + r*r);						// Normalization factor 
	//threaddata->m_sun.SetCoords(r*std::cos(theta)/norm, r*std::sin(theta)/norm, 1.0/norm);
	r = std::sqrt( (*m_rngs)[omp_get_thread_num()]() * m_sineApexAngle );
	theta = (*m_rngs)[omp_get_thread_num()]() * nxmath::TWOPI;
	norm = std::sqrt(1.0 + r*r);
	double x = r*std::cos(theta)/norm;
	double y = r*std::sin(theta)/norm;
	double z = 1.0/norm;
	norm = std::sqrt( x*x + y*y + z*z );
	m_currentSun[omp_get_thread_num()].SetCoords( x/norm, y/norm, z/norm );

}

const HELIODETIC_UNITVECTOR& SKTRAN_Sun_RandomDisc::GetSunUnit( ) const
{

	//m_currentSun1.LookupUpThreadData( &threaddata );
	//return threaddata->m_sun;

	return m_currentSun[omp_get_thread_num()];
}

void SKTRAN_Sun_RandomDisc::SunUnitVector( HELIODETIC_UNITVECTOR* h) const
{

	//m_currentSun1.LookupUpThreadData( &threaddata );
	//h = threaddata->m_sun;

	*h = GetSunUnit();
}

/*
void SKTRAN_Sun_RandomDisc::SunUnitVector( HELIODETIC_VECTOR&     h ) const
{
	//SKTRAN_Sun_RandomDisc_ThreadEntry*	threaddata;

	//m_currentSun1.LookupUpThreadData( &threaddata );
	//h.SetCoords(threaddata->m_sun);

	HELIODETIC_UNITVECTOR hu;
	SunUnitVector( hu );
	h.SetCoords( hu );
}
*/

double SKTRAN_Sun_RandomDisc::CosAngleToSun( const HELIODETIC_UNITVECTOR& h) const
{
	//SKTRAN_Sun_RandomDisc_ThreadEntry*	threaddata;
	//m_currentSun1.LookupUpThreadData( &threaddata );
	//return h.X()*threaddata->m_sun.X() + h.Y()*threaddata->m_sun.Y() + h.Z()*threaddata->m_sun.Z();

	HELIODETIC_UNITVECTOR hu;
	SunUnitVector( &hu );
	return h.X()*hu.X() + h.Y()*hu.Y() + h.Z()*hu.Z();
}

//double SKTRAN_Sun_RandomDisc::CosAngleToSun( const HELIODETIC_VECTOR& h ) const
//{
//	return ComponentToSun(h)/h.Magnitude();
//}

double SKTRAN_Sun_RandomDisc::ComponentToSun( const HELIODETIC_VECTOR& h) const
{
	//SKTRAN_Sun_RandomDisc_ThreadEntry*	threaddata;
	//m_currentSun1.LookupUpThreadData( &threaddata );
	//return h.X()*threaddata->m_sun.X() + h.Y()*threaddata->m_sun.Y() + h.Z()*threaddata->m_sun.Z();

	HELIODETIC_UNITVECTOR hu;
	SunUnitVector( &hu );
	return h.X()*hu.X() + h.Y()*hu.Y() + h.Z()*hu.Z();
}
