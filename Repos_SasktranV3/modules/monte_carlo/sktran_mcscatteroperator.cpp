#include "include/sktran_montecarlo_internals.h"

SKTRAN_MCScatterOperatorContainer::SKTRAN_MCScatterOperatorContainer( )
{
	m_op = nullptr;
	nxLog::Record(NXLOG_ERROR, "SKTRAN_MCScatterOperatorContainer, Cannot call default constructor. This should be =delete but that isn't supported in VS yet.");
}

SKTRAN_MCScatterOperatorContainer::SKTRAN_MCScatterOperatorContainer( const SKTRAN_MCScatterOperatorContainer& other )
{
	m_op = std::shared_ptr<SKTRAN_MCScatterOperator_Base>( other.Op()->ProduceCopy( ) );
}


SKTRAN_MCScatterOperatorContainer::SKTRAN_MCScatterOperatorContainer( std::shared_ptr<SKTRAN_MCScatterOperator_Base>& op )
{
	m_op = op;
}

SKTRAN_MCScatterOperatorContainer::~SKTRAN_MCScatterOperatorContainer( )
{
}

SKTRAN_MCScatterOperator_Base::SKTRAN_MCScatterOperator_Base( )
{
	m_hpfos      = nullptr;
	m_optprops   = nullptr;
	m_mcoptprops = nullptr;
}

SKTRAN_MCScatterOperator_Base::SKTRAN_MCScatterOperator_Base( const SKTRAN_MCScatterOperator_Base& other )
{
	m_hpfos = other.m_hpfos->ProduceCopy();
	m_hpfos->AddRef();

	m_sources.reserve( other.m_sources.size() );
	for( auto s = other.m_sources.begin(); s<other.m_sources.end(); ++s){
		m_sources.push_back(*s);
		(*s)->AddRef();
	}

	m_optprops   = other.m_optprops;
	if(nullptr!=m_optprops) m_optprops->AddRef();
	m_mcoptprops = other.m_mcoptprops;
}

SKTRAN_MCScatterOperator_Base::~SKTRAN_MCScatterOperator_Base( )
{
	ClearSourceTerms( );
	if(NULL!=m_hpfos) m_hpfos->Release( ); m_hpfos = NULL;
}

bool SKTRAN_MCScatterOperator_Base::SetHPFlipOpSymmetry( SKTRAN_HPFOSet* hpfos )
{
	bool ok = true;

	if(nullptr!=  hpfos)   hpfos->AddRef();
	if(nullptr!=m_hpfos) m_hpfos->Release();
	m_hpfos = hpfos;

	return ok;
}
bool SKTRAN_MCScatterOperator_Base::SetOpticalProperties ( const SKTRAN_TableOpticalProperties_Base* optProp )
{
	bool ok = true;

	ok = ok && nullptr!=optProp;
	if(ok) optProp->AddRef();
	ReleaseResources( );
	m_optprops = optProp;
	m_mcoptprops = dynamic_cast< const SKTRAN_TableOpticalProperties_MCBase* >( optProp );

	ok = ok && nullptr!=optProp;
	if(!ok) nxLog::Record(NXLOG_ERROR,"SKTRAN_MCScatterOperator_Base::SetOpticalProperties, Couldn't use opt props for mc scatter operator.");

	return ok;

}

void SKTRAN_MCScatterOperator_Base::ReleaseResources(  )
{
	if(nullptr!=m_optprops) m_optprops->Release();
	m_optprops = nullptr;
	m_mcoptprops = nullptr;
}


bool SKTRAN_MCScatterOperator_Base::AcceptScatterPoint( const SKTRAN_MCPhoton_Base* mcphoton )
{
	bool ok = true;
	
	m_hpfos->AddHPFO( mcphoton );

	return ok;
}

bool SKTRAN_MCScatterOperator_Base::DeclareNewRay( )
{
	return m_hpfos->DeclareNewRay( );
}

bool SKTRAN_MCScatterOperator_Base::InputScatterOrder( int order )
{
	return m_hpfos->InputOrder( order );
}

/*-----------------------------------------------------------------------------
 *			SKTRAN_Engine_MC_V21::GenerateRandomLookVector_groundScatter	2011-07-12*/
/** 
 *	A look direction #newLook is constructed assuming Lambertian scatter
 *	at some point $r\times\hat{r}$, where $\hat{r}=#groundScatterVector$ is
 *	normal to the reflecting surface. 
**/
/*-----------------------------------------------------------------------------*/
void SKTRAN_MCScatterOperator_Base::ChangePhotonBasis_groundScatter(const HELIODETIC_UNITVECTOR& groundScatterVector, const SKTRAN_RNG& rng, SKTRAN_MCPhoton_Base* photon, int order) const
{

	double					cosTheta, sinTheta;	// Theta is chosen from a Lambertian distribution
	double					cosPhi, sinPhi;		// Phi is chosen from a uniform distribution [0, 2pi)
	HELIODETIC_VECTOR       temp;

	// Determine angle between #newLook and #groundScatterVector
	cosTheta = sqrt( rng() );		// Lambertian scatter: P_{cos(\theta)} = sqrt(x), x \in [0 1]
	sinTheta = sqrt( 1.0 - cosTheta*cosTheta );

	// Assume scatter is symmetric in rotation about #groundScatterVector
	double phi = 2.0 * nxmath::Pi * rng();
	cosPhi = std::cos( phi );
	sinPhi = std::sin( phi );
	
	double	norm;
	double	x,y,z;
	x = -groundScatterVector.X();
	y = -groundScatterVector.Y();
	z = -groundScatterVector.Z();
	if( abs(x)<0.99 ){
		norm = sqrt(y*y + z*z);
		photon->GetBasisVar().x.SetCoords(cosTheta*x + sinTheta*(			-sinPhi*norm), 
						  cosTheta*y + sinTheta*( cosPhi*z + sinPhi*x*y)/norm, 
						  cosTheta*z + sinTheta*(-cosPhi*y + sinPhi*x*z)/norm);
	} else{
		norm = sqrt(x*x + z*z);
		photon->GetBasisVar().x.SetCoords(cosTheta*x + sinTheta*(-cosPhi*z + sinPhi*x*y)/norm,
						  cosTheta*y + sinTheta*(			-sinPhi*norm), 
						  cosTheta*z + sinTheta*( cosPhi*x + sinPhi*y*z)/norm);
	}

	// update second component of basis
	if( abs(photon->GetBasis().x.X()) < 0.99 ){
		// [1] == cross(x,[0])
		temp.SetCoords(0.0, photon->GetBasis().x.Z(), -photon->GetBasis().x.Y());
	} else{
		// [1] == cross(y,[0])
		temp.SetCoords( -photon->GetBasis().x.Z(), 0.0, photon->GetBasis().x.X());
	}
	norm = temp.Magnitude();
	photon->GetBasisVar().y.SetCoords(temp.X() / norm, temp.Y() / norm, temp.Z() / norm); 
	photon->GetBasisVar().z.SetCoords(photon->GetBasis().x.Y()*photon->GetBasis().y.Z() - photon->GetBasis().x.Z()*photon->GetBasis().y.Y(),
		                             photon->GetBasis().x.Z()*photon->GetBasis().y.X() - photon->GetBasis().x.X()*photon->GetBasis().y.Z(),
                                     photon->GetBasis().x.X()*photon->GetBasis().y.Y() - photon->GetBasis().x.Y()*photon->GetBasis().y.X() );

	// Phi rotates CCW around b0 -- eta rotates CW around b0
	RotatePolarizedPhoton(photon, cosPhi, -sinPhi, order); 

	return;
}

/*-----------------------------------------------------------------------------
 *			SKTRAN_MCScatterOperator_Base::GenerateRandomLookVector		2011-07-12*/
/** 
 *	A look direction #newLook is constructed using the phase function at
 *	#scatterPoint as a probability density function. Assumes scattering pdf
 *	is uniform in rotation about #prevLook (uniform distribution in azimuth).
**/
/*-----------------------------------------------------------------------------*/
void SKTRAN_MCScatterOperator_Base::ChangePhotonBasis_atmoScatter( const double& cosTheta, SKTRAN_MCPhoton_Base* photon, const SKTRAN_RNG& rng, int order ) const
{
	double tx1,ty1,tz1;
	double tx2,ty2,tz2;
	double tx3,ty3,tz3;
	double phi;
	double sinTheta;
	double cosPhi, sinPhi;
	double norm;
	const SKTRAN_MCBasis& basis = photon->GetBasis();

    sinTheta = std::sqrt(1 - cosTheta*cosTheta);

    // In a time forward sense, the photon was moving in the basis which we'll have at the end of this function, then scattered
    // into some temporary basis, then that basis was rotated through some angle phi into the original basis we see above. 
    // We undo these steps:

	// Find the angle through which basis.y and .z were rotated
	phi = 2.0*nxmath::Pi * rng();
	cosPhi = std::cos( phi );	// Assume distribution of this angle is uniform
	sinPhi = std::sin( phi );

	// t1 was the direction of basis.y before rotation through phi
	tx1 = cosPhi*basis.y.X() - sinPhi*basis.z.X();
	ty1 = cosPhi*basis.y.Y() - sinPhi*basis.z.Y();
	tz1 = cosPhi*basis.y.Z() - sinPhi*basis.z.Z();
	norm = std::sqrt(tx1*tx1 + ty1*ty1 + tz1*tz1); // Machine rounding can be on the order of 1e-6
	tx1 /= norm;
	ty1 /= norm;
	tz1 /= norm;

    // t2 was the direction of basis.z before rotation through phi
    tx2 = cosPhi*basis.z.X() + sinPhi*basis.y.X();
	ty2 = cosPhi*basis.z.Y() + sinPhi*basis.y.Y();
	tz2 = cosPhi*basis.z.Z() + sinPhi*basis.y.Z();
	norm = std::sqrt(tx2*tx2 + ty2*ty2 + tz2*tz2); // Machine rounding can be on the order of 1e-6
	tx2 /= norm;
	ty2 /= norm;
	tz2 /= norm;

    // Apply the rotation to the components of the electric field
	RotatePolarizedPhoton(photon, cosPhi, sinPhi, order); 

	// Now the scattering through theta happens

    // t3 was the direction of basis.y before scattering through theta
    tx3 = cosTheta*tx1 + sinTheta*basis.x.X();
    ty3 = cosTheta*ty1 + sinTheta*basis.x.Y();
    tz3 = cosTheta*tz1 + sinTheta*basis.x.Z();
	norm = std::sqrt(tx3*tx3 + ty3*ty3 + tz3*tz3); // Machine rounding can be on the order of 1e-6
	tx3 /= norm;
	ty3 /= norm;
	tz3 /= norm;

    // t1 was the direction of basis.x before scattering through theta
    tx1 = cosTheta*basis.x.X() - sinTheta*tx1;
    ty1 = cosTheta*basis.x.Y() - sinTheta*ty1;
    tz1 = cosTheta*basis.x.Z() - sinTheta*tz1;
	norm = std::sqrt(tx1*tx1 + ty1*ty1 + tz1*tz1); // Machine rounding can be on the order of 1e-6
	tx1 /= norm;
	ty1 /= norm;
	tz1 /= norm;

    // Set the coordinates
	photon->GetBasisVar().x.SetCoords( tx1, ty1, tz1 );
	photon->GetBasisVar().y.SetCoords( tx3, ty3, tz3 );
	photon->GetBasisVar().z.SetCoords( tx2, ty2, tz2 );

	return;
}

bool SKTRAN_MCScatterOperator_Base::SetCoordinateSystem( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords )
{
	return m_hpfos->SetCoords( coords );
}

void SKTRAN_MCScatterOperator_Base::ClearSourceTerms ( )
{
	for( auto source = m_sources.begin(); source<m_sources.end(); ++source ){
		(*source)->Release();
	}
	m_sources.clear();
}

void SKTRAN_MCScatterOperator_Base::AddSourceTerm ( const SKTRAN_Source_Term* s)
{
	s->AddRef();
	m_sources.push_back( s );
}


SKTRAN_MCScatterOperator_Scalar::SKTRAN_MCScatterOperator_Scalar(  )
{
	m_optprops = nullptr;
}


SKTRAN_MCScatterOperator_Scalar::SKTRAN_MCScatterOperator_Scalar( const SKTRAN_MCScatterOperator_Scalar& other ) : SKTRAN_MCScatterOperator_Base( other )
{
}

SKTRAN_MCScatterOperator_Scalar::~SKTRAN_MCScatterOperator_Scalar(  )
{
	ReleaseResources();
}

void SKTRAN_MCScatterOperator_Scalar::ReleaseResources(  )
{
	if(nullptr!=m_optprops) m_optprops->Release();
	m_optprops = nullptr;
}

bool SKTRAN_MCScatterOperator_Scalar::randomGroundScatter ( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* photon, const SKTRAN_RNG& rng ) const
{
	ChangePhotonBasis_groundScatter ( scatterPoint.Vector().UnitVector(), rng, photon );
	photon->photonOptical()->MoveObserver       ( scatterPoint.Vector(), HELIODETIC_UNITVECTOR(photon->GetBasis().x).Negate() );
	
	return true;
}

bool SKTRAN_MCScatterOperator_Scalar::RandomScatter   ( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* photon, const SKTRAN_RNG& rng, int order ) const
{
	bool ok = true;

	if( photon->m_isGroundScatter ){
		ok = ok && randomGroundScatter( scatterPoint, photon, rng );
	} else{
		ok = ok && randomAtmoScatter  ( scatterPoint, photon, rng, order );
	}

	return ok;
}

bool SKTRAN_MCScatterOperator_Scalar::randomAtmoScatter   ( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* photon, const SKTRAN_RNG& rng, int ) const
{
	bool ok = true;

	double cosTheta, primaryPhase, secondaryPhase;

	//ok = ok && m_mcoptprops->GetCosScatteringAngle ( scatterPoint, rng(), cosTheta, nullptr ); // Don't need the phase matrix
	ok = ok && m_mcoptprops->GetCosScatteringAngle( photon->photonOptical()->GetWavelength(), scatterPoint, rng(), cosTheta, nullptr); // Don't need the phase matrix

	// for an elastic scatter, we must include the ratio of the elastic cross sections (as well as the ratio of the phase functions, which is normally just 1)
	ok = ok && m_optprops->GetScatteringCoefficientCM2(photon->CurrentWavelength(), scatterPoint, cosTheta, &primaryPhase);
	//primaryPhase /= m_optprops->ScatteringExtinctionPerCM(photon->CurrentWavelength(), scatterPoint);
	for (size_t idx = 0; idx < photon->m_numWavelengths; idx++)
	{
		if (idx != photon->m_primaryWavelengthIndex)
		{
			ok = ok && m_optprops->GetScatteringCoefficientCM2(photon->CurrentWavelengths()[idx], scatterPoint, cosTheta, &secondaryPhase);
			//secondaryPhase /= m_optprops->ScatteringExtinctionPerCM(photon->CurrentWavelengths()[idx], scatterPoint);
			photon->ScatterFactors()[idx] = secondaryPhase / primaryPhase;
		}
	}

	if(ok){
		ChangePhotonBasis_atmoScatter          ( cosTheta, photon, rng );
		photon->photonOptical()->MoveObserver   ( scatterPoint.Vector(), HELIODETIC_UNITVECTOR(photon->GetBasis().x).Negate() );
	}
			
	return ok;
}

/* Copied over from the pv2 branch: Add soon */
//
//bool SKTRAN_MCScatterOperator_Polarized::ChooseAzimuthAngle ( const skRTPhaseMatrix& pmatrix, const skRTStokesVector& vec, const SKTRAN_RNG& rng, double& aziAngle ) const
//{
//	bool ok = true;
//
//	double area_sin   = sqrt(pow(pmatrix.At(1,2),2.0)+pow(pmatrix.At(1,3),2.0)) * sqrt(pow(vec.At(2),2.0)+pow(vec.At(3),2.0)); // Area under part of P(azi) that varies as (1+sin(2*azi+phaseShift))
//	double area_const = pmatrix.At(1,1)*vec.At(1) + pmatrix.At(1,4)*vec.At(4) - area_sin; // Area under part of P(azi) that is constant
//	double target = nxmath::TwoPi * rng();
//	double result;
//	const double tolerance = nxmath::TwoPi / 3000.0;
//
//	// I_result = [p11*I + p14*V - sqrt(P12^2+p13^2)*sqrt(Q^2+U^2)]
//	//          + sqrt(P12^2+p13^2)*sqrt(Q^2+U^2) * [1 + sin(2*azi + atan2(Q,U) + atan2(P13,P12) )]
//	double r = rng();
//	//if(r < area_const / (area_sin+area_const) ){
//	if(true){
//		// Azimuth is randomly distributed
//		result = target;
//	} else{
//		// Azimuth is distributed as 1+sin(2*azi+phaseShift)
//		// CDF(azi) = azi + sin(azi)*sin(azi+phaseShift)
//		double phaseShift = atan2(vec.At(2),vec.At(3)) + atan2(pmatrix.At(1,3),pmatrix.At(1,2));
//		double step = nxmath::Pi;
//		result = 0.0;
//		double cdfval = 0.0;
//		double oldcdfval;
//		int ctr = 0;
//		const int maxctr = 1000;
//		while( (tolerance < (target - cdfval)) && ((++ctr)<maxctr) ){ // Do "binary" search to find correct azimuth
//			oldcdfval = cdfval;
//			cdfval = (result+step) + sin(result+step)*sin(result+phaseShift+step);
//			if( cdfval < target ){
//				result = result + step;
//			} else{
//				cdfval = oldcdfval;
//			}
//			step = 0.5 * step;
//		}
//	}
//
//	aziAngle = result;
//
//	return ok;
//}



bool SKTRAN_MCScatterOperator_Scalar::CollectSourceTerms( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* scatterRay, const SKTRAN_CoordinateTransform_V2& coords, size_t threadid ) const
{
	bool ok = true;

	double sourcesRadiance = 0.0;
	double sourcesRadiance2 = 0.0;
	double sourcesVariance = 0.0;

	std::vector< HELIODETIC_POINT > pc;
	std::vector< HELIODETIC_UNITVECTOR > lc;
	ok = ok && m_hpfos->ProducePointCache( scatterPoint, scatterRay->photonOptical()->LookVector(), pc, lc );

    SKTRAN_SourceTermQueryObject_ModifiablePolarized qobj( scatterPoint, scatterRay->photonOptical()->LookVector() );
	
	scatterRay->m_sourcesRadiance = 0.0;
	scatterRay->m_sourcesRadiance2 = 0.0;
	for (auto&& ps : scatterRay->photonSources()) ps.SetVector(0.0);
	for (auto source = m_sources.begin(); source < m_sources.end(); ++source)
	{
		auto l = lc.cbegin(); // Can't put two storage types in the FOR initialization, apparently
		for (auto p = pc.cbegin(); p < pc.end(); ++p, ++l)
		{
			//printf("%1.10e, %1.10e, %1.10e\n"  , p->Vector().X(), p->Vector().Y(), p->Vector().Z());
			//printf("%1.10e, %1.10e, %1.10e\n\n", l->X(), l->Y(), l->Z());
			qobj.UpdateQuery(*p, *l);
			if (scatterRay->m_isGroundScatter) {
				ok = ok && (*source)->MonteCarlo_GroundScatteredRadianceAtPoint(qobj, scatterRay); 
			}
			else {
				ok = ok && (*source)->MonteCarlo_SingleScatteredRadianceAtPoint(qobj, scatterRay);
			}
			sourcesRadiance += scatterRay->photonSource().GetRecentContribSca() / pc.size();
			sourcesRadiance2 += scatterRay->photonSource().GetRecentContribSca() * scatterRay->photonSource().GetRecentContribSca() / pc.size();
		}
	}

	bool elasticRaman = false, done = false;
	while (!done)
	{
		auto ps = scatterRay->photonSources(elasticRaman).cbegin();
		auto sw = scatterRay->ScatterWeights(elasticRaman).cbegin();
		for (auto pr = scatterRay->photonRadiances(elasticRaman).begin(); pr != scatterRay->photonRadiances(elasticRaman).end(); pr++, ps++, sw++)
		{
			pr->AddToVector(ps->GetScalar() * (*sw / pc.size()));
		}
		if (elasticRaman) done = true;
		else elasticRaman = true;
	}
	
	// just use variance from the primary wavelength
	sourcesVariance = 1e-10 < scatterRay->m_sourcesRadiance ? (scatterRay->m_sourcesRadiance2 - scatterRay->m_sourcesRadiance * scatterRay->m_sourcesRadiance) / scatterRay->m_sourcesRadiance : 0.0;
	if (0.0 < sourcesVariance) ok = ok && m_hpfos->InputSourcesVariance(sourcesVariance);	

	return ok;
}


/****************
 * Polarized
 ****************/

SKTRAN_MCScatterOperator_Polarized::SKTRAN_MCScatterOperator_Polarized(  )
{
	m_optprops = nullptr;
}

SKTRAN_MCScatterOperator_Polarized::SKTRAN_MCScatterOperator_Polarized ( const SKTRAN_MCScatterOperator_Polarized& other ) : SKTRAN_MCScatterOperator_Base( other )
{
	m_optprops = other.m_optprops;
	if(nullptr!=m_optprops) m_optprops->AddRef();
}

SKTRAN_MCScatterOperator_Polarized::~SKTRAN_MCScatterOperator_Polarized(  )
{
	ReleaseResources();
}

void SKTRAN_MCScatterOperator_Polarized::ReleaseResources(  )
{
}


bool SKTRAN_MCScatterOperator_Polarized::RandomScatter   ( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* photon, const SKTRAN_RNG& rng, int order ) const
{
	bool ok = true;

	// If all particles are randomly oriented, then all we care about is scattering angle -- reflections are fine across the photon-plane

	if( photon->m_isGroundScatter ){
		ok = ok && randomGroundScatter( scatterPoint, photon, rng, order );
	} else{
		ok = ok && randomAtmoScatter  ( scatterPoint, photon, rng, order );
	}

	return ok;
}

bool SKTRAN_MCScatterOperator_Polarized::randomGroundScatter ( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* photon, const SKTRAN_RNG& rng, int order ) const
{
	SKTRAN_ScatMat_MIMSNC       depolarizeMatrix;
	
	depolarizeMatrix.SetTo(0.0);
	depolarizeMatrix.AssignAt(1,1,1.0);
	photon->AddPhaseOpInPath(depolarizeMatrix);

	// This assumes the polarization basis doesn't have to be changed since we're about to depolarize the photon (rot matrix is 1-2-1 block diagonal)
	ChangePhotonBasis_groundScatter ( scatterPoint.Vector().UnitVector(), rng, photon, order );
	photon->photonOptical()->MoveObserver       ( scatterPoint.Vector(), HELIODETIC_UNITVECTOR(photon->GetBasis().x).Negate() );
	
	return true;
}

bool SKTRAN_MCScatterOperator_Polarized::randomAtmoScatter   ( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* photon, const SKTRAN_RNG& rng, int order ) const
{
	bool ok = true;

	SKTRAN_ScatMat_MIMSNC pmatrix;
	double cosTheta = 0.0;

	ok = ok && m_mcoptprops->GetCosScatteringAngle ( scatterPoint, rng(), cosTheta, &pmatrix );

	if(ok){
		if( 0.0<pmatrix.At(1,1) ){		
		
			ChangePhotonBasis_atmoScatter           ( cosTheta, photon, rng, order ); // Change photon look direction and polarization basis
			photon->photonOptical()->MoveObserver    ( scatterPoint.Vector(), HELIODETIC_UNITVECTOR(photon->GetBasis().x).Negate() ); 

			photon->AddPhaseOpInPath(pmatrix);   // Apply phase matrix for this scatter
		} else{
			// Chose an absorbing point to scatter on 
			pmatrix.SetTo(0.0);
		}
	} else{
		nxLog::Record(NXLOG_ERROR, "SKTRAN_MCScatterOperator_Polarized::randomAtmoScatter, Couldn't get scattering angle, returned %1.16e", cosTheta);
	}

	return ok;
}


bool SKTRAN_MCScatterOperator_PseudoPolarized::RandomScatter   ( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* photon, const SKTRAN_RNG& rng, int order ) const
{
	bool ok = true;

	// If all particles are randomly oriented, then all we care about is scattering angle -- reflections are fine across the photon-plane

	if( photon->m_isGroundScatter ){
		ok = ok && randomGroundScatter( scatterPoint, photon, rng, order );
	} else{
		ok = ok && randomAtmoScatter  ( scatterPoint, photon, rng, order );
	}

	return ok;
}

bool SKTRAN_MCScatterOperator_PseudoPolarized::randomGroundScatter ( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* photon, const SKTRAN_RNG& rng, int order ) const
{
	SKTRAN_ScatMat_MIMSNC       depolarizeMatrix;

	// This assumes the polarization basis doesn't have to be changed since we're about to depolarize the photon (rot matrix is 1-2-1 block diagonal)
	ChangePhotonBasis_groundScatter ( scatterPoint.Vector().UnitVector(), rng, photon, order );
	photon->photonOptical()->MoveObserver       ( scatterPoint.Vector(), photon->GetBasis().x );
	
	
	if(2==order){ //  Hack to get fake vectors
		printf("pseudovector ground scatter off los is not implemented\n");
		//depolarizeMatrix.SetTo(0.0);
		//depolarizeMatrix.At(1,1) = 1.0;
		//photon->AddPhaseOpInPath(depolarizeMatrix);
	}

	return true;
}

bool SKTRAN_MCScatterOperator_PseudoPolarized::randomAtmoScatter   ( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* photon, const SKTRAN_RNG& rng, int order ) const
{
	bool ok = true;

	SKTRAN_ScatMat_MIMSNC pmatrix;
	SKTRAN_Stokes_NC      svec;
	double                cosTheta = 0.0;

	ok = ok && m_mcoptprops->GetCosScatteringAngle ( scatterPoint, rng(), cosTheta, &pmatrix );
	ok = ok && 0.0<pmatrix.At(1,1);

	if(ok){
		if(2==order){ // Hack to get fake vectors
			ok = ok && m_optprops->GetResultOfUnpolarizedScatterCM2( scatterPoint, cosTheta, svec ); 
            svec.Normalize();
			if(ok) m_foovec = svec;
		}
		
		// This part should happen all the time
		//pmatrix *= 1.0/pmatrix.At(1,1);		// Need to get rid of 1/kscatt normalization for MC
		//photon->AddPhaseOpInPath(pmatrix);   // Apply phase matrix for this scatter
		
		ChangePhotonBasis_atmoScatter           ( cosTheta, photon, rng, order ); // Change photon look direction and polarization basis
		photon->photonOptical()->MoveObserver    ( scatterPoint.Vector(), HELIODETIC_UNITVECTOR(photon->GetBasis().x).Negate() ); 
	}
			
	return ok;
}



bool SKTRAN_MCScatterOperator_Polarized::CollectSourceTerms( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* scatterRay, const SKTRAN_CoordinateTransform_V2& coords, size_t threadid ) const
{
	bool ok = true;

	SKTRAN_Stokes_NC svec, radVec, radVec2;
	double eta = 0.0;
	double cosToSource = 0.0;
	SKTRAN_Stokes_NC vectemp; skRTStokesVector::SetToZero(vectemp);
	
	SKTRAN_PhaseMatrixScalar sourcesRadiance  = 0.0;
	SKTRAN_PhaseMatrixScalar sourcesRadiance2 = 0.0;
	double sourcesVariance = 0.0;
	radVec .SetTo ( 0.0 );
	radVec2.SetTo ( 0.0 );

	std::vector< HELIODETIC_POINT > pc;
	std::vector< HELIODETIC_UNITVECTOR > lc;
	std::vector< SKTRAN_MCBasis > bc;
	ok = ok && m_hpfos->ProducePointCache( scatterPoint, scatterRay->photonOptical()->LookVector(), pc, lc );
	ok = ok && m_hpfos->ApplyFlipsToBasis( scatterRay->GetBasis(), bc );
    
    SKTRAN_SourceTermQueryObject_ModifiablePolarized qobj( scatterPoint, scatterRay->photonOptical()->LookVector() );
	for( auto source = m_sources.begin(); source<m_sources.end(); ++source )
	{ 
		auto l = lc.cbegin(); // Can't put two storage types in the FOR initialization
		auto b = bc.cbegin();
		for( auto p = pc.cbegin(); p<pc.end(); ++p, ++l, ++b )
		{
            qobj.UpdateQuery( *p, *l );
            qobj.UpdateBasis( *b );
			if( scatterRay->m_isGroundScatter ){
				ok = ok && (*source)->MonteCarlo_GroundScatteredRadianceAtPoint ( qobj, svec );
			} else{
				ok = ok && (*source)->MonteCarlo_SingleScatteredRadianceAtPoint ( qobj, svec );
			}	

			sourcesRadiance  += svec.I() / pc.size();
			sourcesRadiance2 += svec.I()*svec.I() / pc.size();
			radVec           += svec*(1.0/pc.size());
			radVec2          += svec*(svec.I()/pc.size() );
		}
	}

	sourcesVariance = 1e-10<sourcesRadiance?(sourcesRadiance2 - sourcesRadiance*sourcesRadiance) / sourcesRadiance : 0.0;
	if( 0.0<sourcesVariance ) ok = ok && m_hpfos->InputSourcesVariance( sourcesVariance );
	
	scatterRay->photonRadiance().SetDebugVector ( radVec );
	scatterRay->LeftApplyPathTo(radVec);
	scatterRay->photonRadiance().AddToVector( radVec * scatterRay->ScatterWeight() );

	return ok;
}
bool SKTRAN_MCScatterOperator_PseudoPolarized::CollectSourceTerms( const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* scatterRay, const SKTRAN_CoordinateTransform_V2& coords, size_t threadid ) const
{
	bool ok = true;

	SKTRAN_PhaseMatrixScalar radianceTemp;
	SKTRAN_Stokes_NC         svec, radVec, radVec2;
	double eta = 0.0;
	double cosToSource = 0.0;

	SKTRAN_PhaseMatrixScalar sourcesRadiance  = 0.0;
	SKTRAN_PhaseMatrixScalar sourcesRadiance2 = 0.0;
	double sourcesVariance = 0.0;
	radVec .SetTo ( 0.0 );
	radVec2.SetTo ( 0.0 );

	std::vector< HELIODETIC_POINT > pc;
	std::vector< HELIODETIC_UNITVECTOR > lc;
	std::vector< SKTRAN_MCBasis > bc;
	ok = ok && m_hpfos->ProducePointCache( scatterPoint, scatterRay->photonOptical()->LookVector(), pc, lc );
	ok = ok && m_hpfos->ApplyFlipsToBasis( scatterRay->GetBasis(), bc );
    
    SKTRAN_SourceTermQueryObject_ModifiablePolarized qobj( scatterPoint, scatterRay->photonOptical()->LookVector() );
	for( auto source = m_sources.begin(); source<m_sources.end(); ++source )
	{ 
		auto l=lc.cbegin(); // Can't put two storage types in the FOR initialization
		auto b = bc.cbegin();
		for( auto p = pc.cbegin(); p<pc.end(); ++p, ++l, ++b )
		{
            qobj.UpdateQuery( scatterPoint, scatterRay->photonOptical()->LookVector() );
            qobj.UpdateBasis( *b );
			if( scatterRay->m_isGroundScatter ){
				ok = ok && (*source)->MonteCarlo_GroundScatteredRadianceAtPoint( qobj, radianceTemp );
			} else{
				ok = ok && (*source)->MonteCarlo_SingleScatteredRadianceAtPoint( qobj, radianceTemp );
			}	
			if(1<m_hpfos->GetOrder( )){
				svec = m_foovec;
			} else{
                qobj.UpdateQuery( scatterPoint, scatterRay->photonOptical()->LookVector() );
                qobj.UpdateBasis( *b );
				if( scatterRay->m_isGroundScatter ){
					ok = ok && (*source)->MonteCarlo_GroundScatteredRadianceAtPoint( qobj, svec );
				} else{
					ok = ok && (*source)->MonteCarlo_SingleScatteredRadianceAtPoint( qobj, svec );
				}
				svec *= 0.0 < svec.I() ? 1.0/svec.I() : 1.0;
			}

			sourcesRadiance  += radianceTemp / pc.size();
			sourcesRadiance2 += radianceTemp*radianceTemp / pc.size();
			radVec           += svec*(radianceTemp/pc.size());
			radVec2          += svec*(radianceTemp*radianceTemp/pc.size() );
		}
		// Get contribution for unpolarized source going through the scattering path back to observer
		//sourcesRadiance += radianceTemp;	// Assume all sources are unpolarized
	}

	sourcesVariance = 1e-10<sourcesRadiance?(sourcesRadiance2 - sourcesRadiance*sourcesRadiance) / sourcesRadiance : 0.0;
	if( 0.0<sourcesVariance ) ok = ok && m_hpfos->InputSourcesVariance( sourcesVariance );
	
	scatterRay->LeftApplyPathTo(radVec);
	scatterRay->photonRadiance().AddToVector( radVec       * scatterRay->ScatterWeight());

	return ok;
}

void SKTRAN_MCScatterOperator_Polarized::RotatePolarizedPhoton(SKTRAN_MCPhoton_Base* photon, double cosEta, double sinEta, int) const
{
	photon->AddRotateOpInPath( SKTRAN_ScatMat_Rot( cosEta, sinEta ) );	// Rotate the photon into the observer basis before measuring it

}


void SKTRAN_MCScatterOperator_PseudoPolarized::RotatePolarizedPhoton(SKTRAN_MCPhoton_Base* photon, double cosEta, double sinEta, int order) const
{
	//SKRTFLOAT cosTwoEta = cosEta*cosEta - sinEta*sinEta;
	//SKRTFLOAT sinTwoEta = 2.0*cosEta*sinEta;
	//skRTPhaseMatrix R;

	//if(2==order){
	//	R.SetTo(0.0);
	//	R.At(1,1) =  (SKRTFLOAT)1.0;
	//	R.At(2,2) =  cosTwoEta;
	//	R.At(2,3) =  sinTwoEta;
	//	R.At(3,2) = -sinTwoEta;
	//	R.At(3,3) =  cosTwoEta;
	//	R.At(4,4) =  (SKRTFLOAT)1.0;
	//
	//	photon->AddPhaseOpInPath(R);	// Rotate the photon into the observer basis before measuring it
	//}
}

bool SKTRAN_MCScatterOperator_Polarized::AcceptScatterPoint( const SKTRAN_MCPhoton_Base* mcphoton )
{
	bool ok = true;
	
	ok = ok && SKTRAN_MCScatterOperator_Base::AcceptScatterPoint( mcphoton );

	return ok;
}

bool SKTRAN_MCScatterOperator_Polarized::SetHPFlipOpSymmetry( SKTRAN_HPFOSet* hpfos )
{
	bool ok = true;

	ok = ok && SKTRAN_MCScatterOperator_Base::SetHPFlipOpSymmetry( hpfos );
	ok = ok && m_hpfos->SetMaxNumberFlipOrders( 1 );

	return ok;
}



SKTRAN_MCScatterOperator_ScalarInelastic::SKTRAN_MCScatterOperator_ScalarInelastic()
{
	m_optprops = nullptr;
	m_lastwavelength = 0.0;
}


SKTRAN_MCScatterOperator_ScalarInelastic::SKTRAN_MCScatterOperator_ScalarInelastic(const SKTRAN_MCScatterOperator_ScalarInelastic& other) : SKTRAN_MCScatterOperator_Scalar(other)
{
	m_minwavelength = other.m_minwavelength;
	m_maxwavelength = other.m_maxwavelength;
	m_lastwavelength = other.m_lastwavelength;
}

SKTRAN_MCScatterOperator_ScalarInelastic::~SKTRAN_MCScatterOperator_ScalarInelastic()
{
	ReleaseResources();
}

bool SKTRAN_MCScatterOperator_ScalarInelastic::SetOpticalProperties(const SKTRAN_TableOpticalProperties_Base * optProp)
{
	bool ok = true;
	ok = ok && SKTRAN_MCScatterOperator_Base::SetOpticalProperties(optProp);
	ok = ok && m_mcoptprops->InelasticProperties()->GetWavelengthRange(m_minwavelength, m_maxwavelength);
	return ok;
}

void SKTRAN_MCScatterOperator_ScalarInelastic::ReleaseResources()
{
	if (nullptr != m_optprops) m_optprops->Release();
	m_optprops = nullptr;
}


bool SKTRAN_MCScatterOperator_ScalarInelastic::randomAtmoScatterInelastic(const HELIODETIC_POINT& scatterPoint, const double randNum, const SKTRAN_RNG& rng, SKTRAN_MCPhoton_Base* photon, int order) const
{
	bool ok = true;

	double cosTheta;

	ok = ok && m_mcoptprops->InelasticProperties()->GetCosScatteringAngle(photon->CurrentWavelength(), scatterPoint, rng(), cosTheta);

	ok = ok && m_mcoptprops->InelasticProperties()->GetIncomingWavelength(photon->CurrentWavelengths(), photon->m_primaryWavelengthIndex, scatterPoint, randNum, photon->CurrentWavelengths(), photon->ScatterFactors());

	if (ok) {
		ChangePhotonBasis_atmoScatter(cosTheta, photon, rng);
		photon->photonOptical()->MoveObserver(scatterPoint.Vector(), HELIODETIC_UNITVECTOR(photon->GetBasis().x).Negate());
	}

	return ok;
}



bool SKTRAN_MCScatterOperator_ScalarInelastic::RandomScatter(const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* photon, const SKTRAN_RNG& rng, int order) const
{
	bool ok = true;
	if (photon->m_isGroundScatter)
	{
		ok = ok && SKTRAN_MCScatterOperator_Scalar::RandomScatter(scatterPoint, photon, rng, order);
	}
	else
	{
		bool elastic = true;
		double r = 0.0;

		if (photon->m_manualScatter)
		{
			elastic = photon->m_elasticScatter;
			r = photon->m_randNum;
		}
		else
		{
			double wavelength = photon->photonOptical()->GetWavelength();
			double kelastic = m_optprops->ScatteringExtinctionPerCM(wavelength, scatterPoint);
			double kinelastic = m_mcoptprops->InelasticProperties()->InelasticExtinctionPerCM(wavelength, scatterPoint);
			double k = photon->m_randNum * (kelastic + kinelastic);
			elastic = k > kinelastic;
			if (!elastic) r = k / kinelastic;
		}

		if (elastic)
			ok = ok && SKTRAN_MCScatterOperator_Scalar::RandomScatter(scatterPoint, photon, rng, order);
		else
			ok = ok && randomAtmoScatterInelastic(scatterPoint, r, rng, photon, order);
	}
	return ok;
}

bool SKTRAN_MCScatterOperator_ScalarInelastic::AcceptScatterPoint(const SKTRAN_MCPhoton_Base * mcphoton)
{
	bool ok = true;

	ok = ok && SKTRAN_MCScatterOperator_Base::AcceptScatterPoint(mcphoton);

	if (m_lastwavelength != mcphoton->CurrentWavelength())
	{
		double minwavelength = *std::min_element(mcphoton->CurrentWavelengths().begin(), mcphoton->CurrentWavelengths().end());
		double maxwavelength = *std::max_element(mcphoton->CurrentWavelengths().begin(), mcphoton->CurrentWavelengths().end());
		if (minwavelength < m_minwavelength)
		{
			nxLog::Record(NXLOG_WARNING, "SKTRAN_MCScatterOperator_ScalarInelastic::AcceptScatterPoint, New wavelength %.1f nm is outside the defined range %.1f to %.1f; optical properties are being extrapolated.", minwavelength, m_minwavelength, m_maxwavelength);
		}
		if (maxwavelength > m_maxwavelength)
		{
			nxLog::Record(NXLOG_WARNING, "SKTRAN_MCScatterOperator_ScalarInelastic::AcceptScatterPoint, New wavelength %.1f nm is outside the defined range %.1f to %.1f; optical properties are being extrapolated.", maxwavelength, m_minwavelength, m_maxwavelength);
		}
		m_lastwavelength = mcphoton->CurrentWavelength();
	}

	return ok;
}
