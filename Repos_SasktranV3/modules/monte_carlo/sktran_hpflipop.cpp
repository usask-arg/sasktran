#include "include/sktran_montecarlo_internals.h"

SKTRAN_HVFlipOp::SKTRAN_HVFlipOp( )
{
	for(double* s=m_store; s<(m_store+9); ++s){
		*s = -9999.9;
	}
}

SKTRAN_HVFlipOp::SKTRAN_HVFlipOp( const SKTRAN_HVFlipOp& other )
{
	std::copy(other.m_store, other.m_store+9, m_store);
	//for(int a=0; a<9; ++a){
	//	m_store[a] = other.m_store[a];
	//}
}

SKTRAN_HVFlipOp SKTRAN_HVFlipOp::Define( const SKTRAN_MCPhoton_Base* photon )
{
	nxVector v1( photon->photonOptical()->LookVector().X(),  photon->photonOptical()->LookVector().Y(),  photon->photonOptical()->LookVector().Z()  );
	nxVector v2( photon->photonOptical()->GetObserver().X(), photon->photonOptical()->GetObserver().Y(), photon->photonOptical()->GetObserver().Z() );

    nxVector u = v1.Cross(v2);
	if(1e-10 < u.Magnitude()) u = u.UnitVector();

	// Householder transform
	m_store[0] = 1.0 - 2.0*u.X()*u.X();
	m_store[1] =     - 2.0*u.X()*u.Y();
	m_store[2] =     - 2.0*u.X()*u.Z();
	m_store[3] =     - 2.0*u.Y()*u.X();
	m_store[4] = 1.0 - 2.0*u.Y()*u.Y();
	m_store[5] =     - 2.0*u.Y()*u.Z();
	m_store[6] =     - 2.0*u.Z()*u.X();
	m_store[7] =     - 2.0*u.Z()*u.Y();
	m_store[8] = 1.0 - 2.0*u.Z()*u.Z();

	return *this;
}

HELIODETIC_VECTOR SKTRAN_HVFlipOp::operator* ( const HELIODETIC_VECTOR& v ) const
{
	// Do matrix multiplication
	double x,y,z;
	HELIODETIC_VECTOR result;

	x = m_store[0]*v.X() + m_store[1]*v.Y() + m_store[2]*v.Z();
	y = m_store[3]*v.X() + m_store[4]*v.Y() + m_store[5]*v.Z();
	z = m_store[6]*v.X() + m_store[7]*v.Y() + m_store[8]*v.Z();
	
	result.SetCoords(x, y, z);

	return result;
}

HELIODETIC_UNITVECTOR SKTRAN_HVFlipOp::operator* ( const HELIODETIC_UNITVECTOR& v ) const
{
	// Do matrix multiplication
	double x,y,z;
	HELIODETIC_UNITVECTOR result;

	x = m_store[0]*v.X() + m_store[1]*v.Y() + m_store[2]*v.Z();
	y = m_store[3]*v.X() + m_store[4]*v.Y() + m_store[5]*v.Z();
	z = m_store[6]*v.X() + m_store[7]*v.Y() + m_store[8]*v.Z();
	
	result.SetCoords(x, y, z);

	return result;
}


SKTRAN_MCBasis SKTRAN_HVFlipOp::operator* ( const SKTRAN_MCBasis& b ) const
{
	SKTRAN_MCBasis result;
	
	result.x = (*this)*b.x;
	result.y = (*this)*b.y;
	result.z.SetCoords(result.x.Y()*result.y.Z() - result.x.Z()*result.y.Y(),
		               result.y.X()*result.x.Z() - result.y.Z()*result.x.X(), 
					   result.x.X()*result.y.Y() - result.x.Y()*result.y.X());

	return result;
}

SKTRAN_HVFlipOp SKTRAN_HVFlipOp::LeftApply(const SKTRAN_HVFlipOp& hpfo)
{
	// Left multiply myself by #hpfo
	double r1,r2,r3;
	
	r1 = hpfo.m_store[0]*m_store[0] + hpfo.m_store[1]*m_store[3] + hpfo.m_store[2]*m_store[6];
	r2 = hpfo.m_store[3]*m_store[0] + hpfo.m_store[4]*m_store[3] + hpfo.m_store[5]*m_store[6];
	r3 = hpfo.m_store[6]*m_store[0] + hpfo.m_store[7]*m_store[3] + hpfo.m_store[8]*m_store[6];
	m_store[0] = r1; m_store[3] = r2; m_store[6] = r3;
	r1 = hpfo.m_store[0]*m_store[1] + hpfo.m_store[1]*m_store[4] + hpfo.m_store[2]*m_store[7];
	r2 = hpfo.m_store[3]*m_store[1] + hpfo.m_store[4]*m_store[4] + hpfo.m_store[5]*m_store[7];
	r3 = hpfo.m_store[6]*m_store[1] + hpfo.m_store[7]*m_store[4] + hpfo.m_store[8]*m_store[7];
	m_store[1] = r1; m_store[4] = r2; m_store[7] = r3;
	r1 = hpfo.m_store[0]*m_store[2] + hpfo.m_store[1]*m_store[5] + hpfo.m_store[2]*m_store[8];
	r2 = hpfo.m_store[3]*m_store[2] + hpfo.m_store[4]*m_store[5] + hpfo.m_store[5]*m_store[8];
	r3 = hpfo.m_store[6]*m_store[2] + hpfo.m_store[7]*m_store[5] + hpfo.m_store[8]*m_store[8];
	m_store[2] = r1; m_store[5] = r2; m_store[8] = r3;

	return *this;
};

SKTRAN_HPFOSet::SKTRAN_HPFOSet( )
{
	m_coords = NULL;
}

SKTRAN_HPFOSet::~SKTRAN_HPFOSet( )
{
}

bool SKTRAN_HPFOSet::DeclareNewRay( )
{
	m_hvfos.clear(); 
	return true;
}

bool SKTRAN_HPFOSet::InputOrder( int order )
{
	m_order = order;
	return true;
}

bool SKTRAN_HPFOSet::CopyInto( SKTRAN_HPFOSet& dest ) const
{
	bool ok = true;

	dest.m_coords = m_coords;
	std::copy( m_hvfos.begin(), m_hvfos.end(), dest.m_hvfos.begin() );

	return ok;
}

bool SKTRAN_HPFOSet::SetCoords( std::shared_ptr< const SKTRAN_CoordinateTransform_V2>& coords )
{
	bool ok = true;
	m_coords = coords;

	return ok;
}

bool SKTRAN_HPFOSet_NoSymmetry::ProducePointCache( const HELIODETIC_POINT& hp, const HELIODETIC_UNITVECTOR& look, std::vector<HELIODETIC_POINT>& pc, std::vector< HELIODETIC_UNITVECTOR >& lc ) const
{
	bool ok = true;
	
	pc.clear    ( );
	pc.reserve  ( 1 );
	pc.push_back( hp );

	lc.clear    ( );
	lc.reserve  ( 1 );
	lc.push_back( look );
	return ok;
}

bool SKTRAN_HPFOSet_NoSymmetry::ApplyFlipsToBasis ( const SKTRAN_MCBasis& basis, std::vector<SKTRAN_MCBasis>& flippedBasis ) const
{
	bool ok = true;

	flippedBasis.clear( );
	flippedBasis.reserve( 1 );
	flippedBasis.push_back( basis );

	return ok;
}
		
bool SKTRAN_HPFOSet_NoSymmetry::InputSourcesVariance ( double )
{
	return true;
}


bool SKTRAN_HPFOSet_NoSymmetry::AddHPFO(const SKTRAN_HVFlipOp& hpfo)
{
	bool ok = true;

	return true;
}


bool SKTRAN_HPFOSet_NoSymmetry::AddHPFO(const SKTRAN_MCPhoton_Base* mcphoton)
{
	bool ok = true;

	return ok;
}

bool SKTRAN_HPFOSet_NoSymmetry::SetMaxNumberFlipOrders ( int maxOrder )
{
	bool ok = true;

	return ok;
}

SKTRAN_HPFOSet_HorizSymmetry::SKTRAN_HPFOSet_HorizSymmetry( )
{
	m_order_stopPushing             = 1;
	m_order_stopPushing_max         = 10;
	m_order_clear                   = 5;
	m_numRunsSinceCheck_stopPushing = 0;
	m_numRunsSinceCheck_clear       = 0;
	m_numRunsBetweenChecks          = 64;
	m_lastReportedVariance          = 0.0;
	m_runningvariance_stopPushing   = 0.0;
	m_runningvariance_clear         = 0.0;
	m_highThreshold                 = 0.05*0.05;
	m_lowThreshold                  = 0.01*0.01;
}


SKTRAN_HPFOSet_HorizSymmetry::SKTRAN_HPFOSet_HorizSymmetry( const SKTRAN_HPFOSet_HorizSymmetry& other )
{
	m_order_stopPushing             = other.m_order_stopPushing;
	m_order_stopPushing_max         = other.m_order_stopPushing_max;
	m_order_clear                   = other.m_order_clear;
	m_numRunsSinceCheck_stopPushing = other.m_numRunsSinceCheck_stopPushing;
	m_numRunsSinceCheck_clear       = other.m_numRunsSinceCheck_clear;
	m_numRunsBetweenChecks          = other.m_numRunsBetweenChecks;
	m_lastReportedVariance          = other.m_lastReportedVariance;
	m_runningvariance_stopPushing   = other.m_runningvariance_stopPushing;
	m_runningvariance_clear         = other.m_runningvariance_clear;
	m_highThreshold                 = other.m_highThreshold;
	m_lowThreshold                  = other.m_lowThreshold;
	other.CopyInto(*this);
	m_variance.resize(m_order_clear-m_order_stopPushing);
	for( auto it=m_variance.begin(); it<m_variance.end(); ++it){
		*it = 0.0;
	}

}

bool SKTRAN_HPFOSet_HorizSymmetry::ProducePointCache( const HELIODETIC_POINT& hp, const HELIODETIC_UNITVECTOR& look, std::vector<HELIODETIC_POINT>& pc, std::vector< HELIODETIC_UNITVECTOR >& lc ) const
{
	bool ok = true;

	HELIODETIC_UNITVECTOR u      = hp.Vector().UnitVector();
	double                radius = hp.Radius();
	HELIODETIC_UNITVECTOR tempU;
	HELIODETIC_POINT      tempP;

	pc.clear();
	lc.clear();
	pc.reserve(m_hvfos.size()+1);
	lc.reserve(m_hvfos.size()+1);

	pc.push_back( hp );
	lc.push_back( look );
	for( std::vector<SKTRAN_HVFlipOp>::const_iterator hvfo = m_hvfos.begin(); hvfo<m_hvfos.end(); ++hvfo){
		tempU =  (*hvfo)*u ;
		tempP.Initialize( tempU, radius, m_coords.get() );
		pc.push_back( tempP );
		tempU = (*hvfo)*look;
		lc.push_back(tempU);
	}

	return ok;
}

bool SKTRAN_HPFOSet_HorizSymmetry::ApplyFlipsToBasis ( const SKTRAN_MCBasis& basis, std::vector<SKTRAN_MCBasis>& flippedBasis ) const
{
	bool ok = true;

	flippedBasis.reserve( m_hvfos.size()+1 );
	flippedBasis.push_back( basis );
	for( std::vector<SKTRAN_HVFlipOp>::const_iterator hvfo = m_hvfos.begin(); hvfo<m_hvfos.end(); ++hvfo ){
		flippedBasis.push_back( (*hvfo)*basis );
	}
	return ok;
}

bool SKTRAN_HPFOSet_HorizSymmetry::DeclareNewRay( )
{
	++m_numRunsSinceCheck_stopPushing;
	return SKTRAN_HPFOSet::DeclareNewRay( );
}
bool SKTRAN_HPFOSet_HorizSymmetry::AdjustLimits( )
{
	bool ok = true;
	
	double aveOp  = 1.0 / m_numRunsSinceCheck_stopPushing;

	for( auto it=m_variance.begin(); it<m_variance.end(); ++it ){
		*it = *it * aveOp;
	}
	double minVar = 0.01*0.01*pow(2.0,m_order_stopPushing-1); // Minimum variance for which it is beneficial to do flips
	double maxVar = 0.05*0.01*pow(2.0,m_order_stopPushing-1); // Minimum variance for which we should try doing more flips
	if( maxVar < m_variance.front() && m_order_stopPushing<m_order_stopPushing_max ) ++m_order_stopPushing;
	else if(1<m_order_stopPushing && m_variance.front() < minVar) --m_order_stopPushing;
	if( maxVar < m_variance.back() || m_order_clear<m_order_stopPushing ) ++m_order_clear;
	else if( m_order_stopPushing<m_order_clear && m_variance.back()<minVar ) --m_order_clear;

	m_variance.resize(m_order_clear-m_order_stopPushing + (m_order_clear==m_order_stopPushing?1:0) );
	for( auto it=m_variance.begin(); it<m_variance.end(); ++it){
		*it = 0.0;
	}

	return ok;
}

bool SKTRAN_HPFOSet_HorizSymmetry::InputSourcesVariance( double variance )
{
	bool ok = true;

	if( ( m_order_stopPushing<m_order && m_order<=m_order_clear+1 ) ){
		m_variance[m_order-m_order_stopPushing-1] += variance;
		++m_numRunsSinceCheck_stopPushing;;
	} else if( m_order_stopPushing==m_order && m_order_stopPushing==m_order_clear ){
		m_variance[m_order-m_order_stopPushing  ] += variance;
		++m_numRunsSinceCheck_stopPushing;
	}

	if( m_numRunsBetweenChecks < m_numRunsSinceCheck_stopPushing ){
		ok = ok && AdjustLimits( );
		m_numRunsSinceCheck_stopPushing = 0;
	}
	return true;
}


bool SKTRAN_HPFOSet_HorizSymmetry::AddHPFO(const SKTRAN_HVFlipOp& hpfo)
{
	bool ok = true;

	if( m_order <= m_order_stopPushing ){
		m_hvfos.reserve(2*m_hvfos.size()+1);
		m_hvfos.push_back(hpfo);
		auto opidx = m_hvfos.size()-1;
		for( std::vector<SKTRAN_HVFlipOp>::const_iterator myhpfo = m_hvfos.cbegin(); 0<opidx; ++myhpfo, --opidx ){
			//m_hvfos.emplace_back(*myhpfo);
			//m_hvfos.back().LeftApply(hpfo);
			m_hvfos.emplace_back(hpfo);
			m_hvfos.back().LeftApply(*myhpfo);
		}
	} else if( m_order_clear == m_order ){
		m_hvfos.clear();
	}

	return ok;
}

bool SKTRAN_HPFOSet_HorizSymmetry::AddHPFO(const SKTRAN_MCPhoton_Base* mcphoton)
{
	bool ok = true;

	SKTRAN_HVFlipOp hpfo;

	return ok && AddHPFO ( hpfo.Define(mcphoton) );
}


bool SKTRAN_HPFOSet_HorizSymmetry::SetMaxNumberFlipOrders ( int maxOrder )
{
	bool ok = true;

	ok = ok && 0<maxOrder;
	if(ok) m_order_stopPushing_max = maxOrder;
	if(m_order_stopPushing_max < m_order_stopPushing ) m_order_stopPushing = m_order_stopPushing_max;

	return ok;
}