#include "include/sktran_montecarlo_internals.h"

SKTRAN_PhotonLog_AveKernel::SKTRAN_PhotonLog_AveKernel( )
{
	m_altTolerance = 1.0;

	m_minAlt        = 0.0;
	m_numAlts       = 200;
	m_deltaAlt      = 500.0;
	
	m_minCosAngle   = -1.0;
	m_numCosAngles  = 201;
	m_deltaCosAngle = 0.01; 
	m_grid.resize(0);

}


SKTRAN_PhotonLog_AveKernel::~SKTRAN_PhotonLog_AveKernel ( )
{

}


bool SKTRAN_PhotonLog_AveKernel::ConfigureKernel ( double minAlt, size_t numAlts, double deltaAlt, double minCosAngle, size_t numCosAngles, double deltaCosAngle, size_t numThreads )
{
	bool ok = true;

	ok = ok && 0.0<deltaAlt && 0.0<deltaCosAngle && 0<numThreads;

	if(ok){
		m_minAlt        = minAlt;
		m_numAlts       = numAlts;
		m_deltaAlt      = deltaAlt;
		m_minCosAngle   = minCosAngle;
		m_numCosAngles  = numCosAngles;
		m_deltaCosAngle = deltaCosAngle;
		
		m_grid.resize( numThreads );

		for( std::vector< std::vector<double> >::iterator tidx=m_grid.begin(); tidx<m_grid.end(); tidx++)
		{
			tidx->resize( (numAlts+1) * numCosAngles ); // Number of altitudes, plus ground scatter
		}

		WipeKernel( );

	} else{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_AveragingKernel::ConfigureKernel, Requires 0.0<delta, 0<numThreads.");
	}
	
	if(!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_AveragingKernel::ConfigureKernel, Could not configure kernel.");
	}

	return ok;

}


bool SKTRAN_PhotonLog_AveKernel::ConfigureObserverGeometry( const SKTRAN_RayOptical_Base& observerRay, double earthRadius )
{
	bool ok = true;

	const SKTRAN_RayOptical_Straight* straightOpt = reinterpret_cast<const SKTRAN_RayOptical_Straight*>(&observerRay);
	
	HELIODETIC_UNITVECTOR look  = straightOpt->LookVector();
//	HELIODETIC_VECTOR     tangentPoint = straightOpt->GeometryRay()->GetObserver() + HELIODETIC_VECTOR(straightOpt->GeometryRay()->LookVector(),straightOpt->GeometryRay()->DistanceToTangentPoint(0.0));
	HELIODETIC_VECTOR     tangentPoint = straightOpt->GetObserver() + HELIODETIC_VECTOR(straightOpt->LookVector(), std::sqrt(std::pow(straightOpt->GetObserver().Magnitude(),2) - std::pow(straightOpt->GetObserver() & straightOpt->LookVector(), 2)) ); 
		
	m_earthRadius = earthRadius;
	m_lookUnit.SetCoords  ( look.X(), look.Y(), look.Z() );
	m_upUnit.SetCoords    ( tangentPoint.UnitVector().X(), tangentPoint.UnitVector().Y(), tangentPoint.UnitVector().Z() );
	m_rightUnit.SetCoords ( m_lookUnit.Y()*m_upUnit.Z() - m_lookUnit.Z()*m_upUnit.Y(), 
		                    m_lookUnit.Z()*m_upUnit.X() - m_lookUnit.X()*m_upUnit.Z(), 
		                    m_lookUnit.X()*m_upUnit.Y() - m_lookUnit.Y()*m_upUnit.X() ); // look cross up
	
	ok = ok && 0.999 < m_lookUnit.Magnitude()   && 1.001 > m_lookUnit.Magnitude() &&
		       0.999 < m_upUnit.Magnitude()     && 1.001 > m_upUnit.Magnitude()   &&
		       0.999 < m_rightUnit.Magnitude()  && 1.001 > m_rightUnit.Magnitude();
	if( !ok )
	{
		nxLog::Record(NXLOG_ERROR, "SKTRAN_AveragingKernel::ConfigureObserverGeometry, Could not create orthogonal system -- be careful using curved rays.");
	}

	return ok;
}


void SKTRAN_PhotonLog_AveKernel::WipeKernel( )
{
	for(std::vector<std::vector<double> >::iterator tidx=m_grid.begin(); tidx<m_grid.end(); tidx++)
	{
		for(std::vector<double>::iterator gidx=tidx->begin(); gidx<tidx->end(); gidx++)
		{
			*gidx = 0.0;
		}
	}
}


bool SKTRAN_PhotonLog_AveKernel::FindGridWeights ( HELIODETIC_VECTOR vec, bool groundScatter, size_t* indices, double* weights )
{
	bool ok = true;

	double cosAlongLook;
	double lowAltWeight;
	double lowAngWeight;
	double planeScaleFactor =  sqrt(1.0-pow(m_rightUnit & vec.UnitVector(),2));
	size_t lowAng;
	
	// Get altitude weights
	double alt = vec.Magnitude() - m_minAlt - m_earthRadius;
	size_t lowAlt = (size_t) floor(alt / m_deltaAlt);
	if(groundScatter)
	{
		lowAlt = 0;
		lowAltWeight = 1.0;
	} else{
		if(lowAlt < m_numAlts && 0.0<alt)
		{
			lowAltWeight = 1 - ( (alt-lowAlt*m_deltaAlt) / m_deltaAlt );
		} else if(0.0 < (alt+m_altTolerance)){ // Rounding error caused scatter just below the Earth
			lowAlt = 0;
			lowAltWeight = 1.0;
		} else if( m_numAlts*m_deltaAlt + m_altTolerance < vec.Magnitude() ){ // Rounding error caused scatter just above the atmosphere
			lowAlt = m_numAlts-2;
			lowAltWeight = 0.0;
		} else{
			ok = false;
			nxLog::Record(NXLOG_ERROR, "SKTRAN_AveragingKernel::FindGridWeights, Altitude out of range.");
		}
		++lowAlt;
	}

	// Get angular component
	if(1e-16<planeScaleFactor){
		planeScaleFactor = 1.0/planeScaleFactor;
		cosAlongLook = (m_lookUnit&vec.UnitVector())*planeScaleFactor - m_minCosAngle;
		lowAng = (size_t)floor(cosAlongLook / m_deltaCosAngle);
		if(lowAng < m_numCosAngles){
			lowAngWeight = 1 - ( (cosAlongLook-lowAng*m_deltaCosAngle) / m_deltaCosAngle );
		} else{
			ok = false;
			nxLog::Record(NXLOG_ERROR, "SKTRAN_AveragingKernel::FindGridWeights, Angle out of range.");
		}
	} else{
		ok = false;
		nxLog::Record(NXLOG_ERROR,"SKTRAN_AveragingKernel::FindGridWeights, Grid position is not unique, ignoring entry.");
	}


	if(ok){
		indices[0] = m_numCosAngles*(lowAlt  ) + lowAng  ;
		indices[1] = m_numCosAngles*(lowAlt  ) + lowAng+1;
		indices[2] = m_numCosAngles*(lowAlt+1) + lowAng  ;
		indices[3] = m_numCosAngles*(lowAlt+1) + lowAng+1;
		weights[0] = (  lowAltWeight)*(  lowAngWeight);
		weights[1] = (  lowAltWeight)*(1-lowAngWeight);
		weights[2] = (1-lowAltWeight)*(  lowAngWeight);
		weights[3] = (1-lowAltWeight)*(1-lowAngWeight);
	}

	return ok;
}

bool SKTRAN_PhotonLog_AveKernel::AddToKernel ( const SKTRAN_MCPhoton_Base* photon, size_t order, size_t threadid)
				//ok = ok && m_aveKernel->AddToKernel(mcphoton->m_scatterVector, mcphoton->photonRadiance().GetRecentContribVec(), order-1, mcphoton->m_isGroundScatter, threadid);
{
	bool ok = true;

	size_t indices[4];
	double weights[4];

	ok = ok && FindGridWeights( photon->m_scatterVector, photon->m_isGroundScatter, indices, weights );
	if(ok){
		m_grid[threadid][indices[0]] += photon->photonRadiance().GetRecentContribVec().I()*weights[0];
		m_grid[threadid][indices[1]] += photon->photonRadiance().GetRecentContribVec().I()*weights[1];
		m_grid[threadid][indices[2]] += photon->photonRadiance().GetRecentContribVec().I()*weights[2];
		m_grid[threadid][indices[3]] += photon->photonRadiance().GetRecentContribVec().I()*weights[3];
	}

	return ok;
}


bool SKTRAN_PhotonLog_AveKernel::PrintKernel ( std::string filenameNoExtension )
{
	bool ok = true;

	std::vector<double> toPrint;
	toPrint.resize(m_grid[0].size());

	// Get contributions from each thread
	for(std::vector<double>::iterator it=toPrint.begin(); it<toPrint.end(); ++it){
		*it = 0.0;
	}

	std::vector<double>::iterator dest;

	for(std::vector<std::vector<double> >::iterator it=m_grid.begin(); it<m_grid.end(); ++it){
		dest = toPrint.begin();
		for(std::vector<double>::iterator src=it->begin(); src<it->end(); ++src, ++dest){
			*dest += *src;
		}
	}

	// Print params to file
	FILE * f;
	f = fopen( (filenameNoExtension+".txt").c_str(), "w" );
	dest = toPrint.begin();
	if(ok && NULL!=f){
		for(size_t altidx=0; altidx<(m_numAlts+1); ++altidx){
			for(size_t angidx=0; angidx<m_numCosAngles; ++angidx){
				std::fprintf(f, "%1.8e\t", *dest);
				++dest;
			}
			if(m_numAlts!=altidx) std::fprintf(f,"\n");
		}
		std::fclose(f);
	} else{
		ok = false;
		nxLog::Record(NXLOG_ERROR, "SKTRAN_AveragingKernel::PrintKernel, Cannot write destination file.");
	}

	// Write params to file
	f = fopen( (filenameNoExtension+"_params.txt").c_str(), "w" );
	dest = toPrint.begin();
	if(ok && NULL!=f){
		std::fprintf(f,"m_earthRadius:   %1.15e\n", m_earthRadius);
		std::fprintf(f,"m_minAlt:        %1.15e\n", m_minAlt);
		std::fprintf(f,"m_deltaAlt:      %1.15e\n", m_deltaAlt);
		std::fprintf(f,"m_minCosAngle:   %1.15e\n", m_minCosAngle);
		std::fprintf(f,"m_deltaCosAngle: %1.15e\n", m_deltaCosAngle);
		std::fprintf(f,"m_numAlts:       %u\n"    , (unsigned int)m_numAlts);
		std::fprintf(f,"m_numCosAngles:  %u"      , (unsigned int)m_numCosAngles);
	} else{
		ok = false;
		nxLog::Record(NXLOG_ERROR, "SKTRAN_AveragingKernel::PrintKernel, Cannot write params destination file.");
	}
	std::fclose(f);

	return ok;
}








bool SKTRAN_PhotonLog_RadianceOnLos::ConfigureObserverGeometry   ( const SKTRAN_RayOptical_Base& observerRay, double earthRadius )
{
	bool ok = true;
	m_observer = observerRay.GetObserver();
	
	m_quadPtDistances.resize( observerRay.Storage()->NumQuadraturePoints() );
	for(int tidx=0; tidx<m_weights.size(); ++tidx){
		m_weights[tidx].resize( m_quadPtDistances.size() );
		m_vals[tidx]   .resize( m_quadPtDistances.size() );
	}

	WipeKernel();

	for(int qidx=0; qidx<observerRay.Storage()->NumQuadraturePoints(); ++qidx ){
		m_quadPtDistances[qidx] = observerRay.Storage()->DistanceOfPointFromOrigin( qidx );
	}



	return ok;
}

bool SKTRAN_PhotonLog_RadianceOnLos::ConfigureKernel ( double minAlt, size_t numAlts, double deltaAlt, double minCosAngle, size_t numCosAngles, double deltaCosAngle, size_t numThreads )
{
	bool ok = true;

	m_next_lowerIndex .resize( numThreads );
	m_next_lowerWeight.resize( numThreads );
	m_next_upperWeight.resize( numThreads );

	m_weights.resize(numThreads);
	m_vals.resize(numThreads);

	for(int tidx=0; tidx<numThreads; ++tidx){
		m_weights[tidx].resize( m_quadPtDistances.size() );
		m_vals[tidx]   .resize( m_quadPtDistances.size() );
	}
	
	WipeKernel();

	return ok;
}

void SKTRAN_PhotonLog_RadianceOnLos::WipeKernel ( )
{

	m_next_lowerIndex .resize( m_next_lowerIndex.size(),    0 );
	m_next_lowerWeight.resize( m_next_lowerWeight.size(), 0.0 );
	m_next_upperWeight.resize( m_next_upperWeight.size(), 0.0 );

	for(int tidx=0; tidx<m_weights.size(); ++tidx){
		for(int qidx=0; qidx<m_weights[tidx].size(); ++qidx ){
			m_weights[tidx][qidx]  = 0.0;
			m_vals[tidx][qidx].SetTo(0.0);
		}
	}
}

bool SKTRAN_PhotonLog_RadianceOnLos::AddToKernel ( const SKTRAN_MCPhoton_Base* photon, size_t order, size_t threadid)
				//ok = ok && m_aveKernel->AddToKernel(mcphoton->m_scatterVector, mcphoton->photonRadiance().GetRecentContribVec(), order-1, mcphoton->m_isGroundScatter, threadid);
{
	bool ok = true;

	if(1==order){
		double dist = (photon->m_scatterVector-m_observer).Magnitude();
		size_t ptidx = std::distance( m_quadPtDistances.cbegin(), std::upper_bound(m_quadPtDistances.cbegin(), m_quadPtDistances.cend(), dist ) )-1;
		double upperWeight, lowerWeight;

		if( ptidx < (m_quadPtDistances.size()-1) ){
			upperWeight = (dist-m_quadPtDistances[ptidx]) / (m_quadPtDistances[ptidx+1]-m_quadPtDistances[ptidx]);
			lowerWeight = 1.0 - upperWeight;
		} else{
			ptidx = 0;
			upperWeight = 0.0;
			lowerWeight = 0.0;
		}
		m_next_lowerIndex [threadid] = ptidx;
		m_next_lowerWeight[threadid] = lowerWeight;
		m_next_upperWeight[threadid] = upperWeight;

	} else if(2==order){

		size_t ptidx = m_next_lowerIndex[threadid];
		m_weights[threadid][ptidx  ] += m_next_lowerWeight[threadid];
		m_weights[threadid][ptidx+1] += m_next_upperWeight[threadid];
		m_vals   [threadid][ptidx  ] += photon->photonRadiance().GetRecentContribVec() * m_next_lowerWeight[threadid];
		m_vals   [threadid][ptidx+1] += photon->photonRadiance().GetRecentContribVec() * m_next_upperWeight[threadid];
	} 

	return ok;
}

bool SKTRAN_PhotonLog_RadianceOnLos::PrintKernel ( std::string filenameNoExtension )
{
	bool ok = true;

	// Print radiance contributions along LOS
	{
		std::vector< SKTRAN_Stokes_NC > result_vals;
		std::vector< double > result_weights;
		result_vals.resize( m_vals[0].size() );
		for( int ridx=0; ridx<result_vals.size(); ++ridx) result_vals[ridx].SetTo(0.0);
		result_weights.resize( m_weights[0].size(), 0.0 );

		for(int tidx=0; tidx<m_vals.size(); ++tidx){
			for(int qidx=0; qidx<result_vals.size(); ++qidx){
				result_vals[qidx] += m_vals[tidx][qidx];
				result_weights[qidx] += m_weights[tidx][qidx];
			}
		}
		//for(int qidx=0; qidx<result_vals.size(); ++qidx){
		//	if( 0.0 < result_weights[qidx] ) result_vals[qidx] *= (1.0/result_weights[qidx]);
		//}

		FILE * f;
		f = fopen( (filenameNoExtension+".txt").c_str(), "w" );
		ok = ok && nullptr!=f;
		if(ok){
			for( int qidx=0; qidx<result_vals.size(); ++qidx){
				fprintf(f, "%1.16e, %1.16e, %1.16e, %1.16e, %1.16e, %1.16e", m_quadPtDistances[qidx], result_weights[qidx], result_vals[qidx].I(), result_vals[qidx].Q(), result_vals[qidx].U(), result_vals[qidx].V() );
				if(qidx<(result_vals.size()-1)) fprintf(f, "\n");
			}
		}
		fclose(f);
	}

	return ok;
}







bool SKTRAN_PhotonLog_PhotonsOnLos::ConfigureObserverGeometry   ( const SKTRAN_RayOptical_Base& observerRay, double earthRadius )
{
	bool ok = true;
	m_observer = observerRay.GetObserver();
	
	m_quadPtDistances.resize( observerRay.Storage()->NumQuadraturePoints() );
	for(int tidx=0; tidx<m_weights.size(); ++tidx){
		m_weights[tidx].resize( m_quadPtDistances.size() );
		m_vals[tidx]   .resize( m_quadPtDistances.size() );
	}

	WipeKernel();

	for(int qidx=0; qidx<observerRay.Storage()->NumQuadraturePoints(); ++qidx ){
		m_quadPtDistances[qidx] = observerRay.Storage()->DistanceOfPointFromOrigin( qidx );
	}



	return ok;
}

bool SKTRAN_PhotonLog_PhotonsOnLos::ConfigureKernel ( double minAlt, size_t numAlts, double deltaAlt, double minCosAngle, size_t numCosAngles, double deltaCosAngle, size_t numThreads )
{
	bool ok = true;

	m_next_lowerIndex .resize( numThreads );
	m_next_lowerWeight.resize( numThreads );
	m_next_upperWeight.resize( numThreads );

	m_weights.resize(numThreads);
	m_vals.resize(numThreads);

	for(int tidx=0; tidx<numThreads; ++tidx){
		m_weights[tidx].resize( m_quadPtDistances.size() );
		m_vals[tidx]   .resize( m_quadPtDistances.size() );
	}
	
	WipeKernel();

	return ok;
}

void SKTRAN_PhotonLog_PhotonsOnLos::WipeKernel ( )
{

	m_next_lowerIndex .resize( m_next_lowerIndex.size(),    0 );
	m_next_lowerWeight.resize( m_next_lowerWeight.size(), 0.0 );
	m_next_upperWeight.resize( m_next_upperWeight.size(), 0.0 );

	for(int tidx=0; tidx<m_weights.size(); ++tidx){
		for(int qidx=0; qidx<m_weights[tidx].size(); ++qidx ){
			m_weights[tidx][qidx].clear();
			m_vals[tidx][qidx].clear();
		}
	}
}

bool SKTRAN_PhotonLog_PhotonsOnLos::AddToKernel ( const SKTRAN_MCPhoton_Base* photon, size_t order, size_t threadid)
				//ok = ok && m_aveKernel->AddToKernel(mcphoton->m_scatterVector, mcphoton->photonRadiance().GetRecentContribVec(), order-1, mcphoton->m_isGroundScatter, threadid);
{
	bool ok = true;

	if(1==order){
		double dist = (photon->m_scatterVector-m_observer).Magnitude();
		size_t ptidx = std::distance( m_quadPtDistances.cbegin(), std::upper_bound(m_quadPtDistances.cbegin(), m_quadPtDistances.cend(), dist ) )-1;
		double upperWeight, lowerWeight;

		if( ptidx < (m_quadPtDistances.size()-1) ){
			upperWeight = (dist-m_quadPtDistances[ptidx]) / (m_quadPtDistances[ptidx+1]-m_quadPtDistances[ptidx]);
			lowerWeight = 1.0 - upperWeight;
		} else{
			ptidx = 0;
			upperWeight = 0.0;
			lowerWeight = 0.0;
		}
		m_next_lowerIndex [threadid] = ptidx;
		m_next_lowerWeight[threadid] = lowerWeight;
		m_next_upperWeight[threadid] = upperWeight;

	} else if(2==order){

		size_t ptidx = m_next_lowerIndex[threadid];
		m_weights[threadid][ptidx  ].push_back( m_next_lowerWeight[threadid] );
		m_weights[threadid][ptidx+1].push_back( m_next_upperWeight[threadid] );
		m_vals   [threadid][ptidx  ].push_back( photon->photonRadiance().GetRecentContribVec() * m_next_lowerWeight[threadid] );
		m_vals   [threadid][ptidx+1].push_back( photon->photonRadiance().GetRecentContribVec() * m_next_upperWeight[threadid] );
	} 

	return ok;
}

bool SKTRAN_PhotonLog_PhotonsOnLos::PrintKernel ( std::string filenameNoExtension )
{
	bool ok = true;

	// Print radiance contributions along LOS
	{
		FILE* qf;
		qf = fopen( (filenameNoExtension+"_distances.txt").c_str(), "w" );
        ok = ok && nullptr!=qf;
        if(ok){
		    for( int qidx=0; 0<m_vals.size() && qidx<m_vals[0].size(); ++qidx ){
			    FILE * f;
			    stringstream ss; ss << qidx;
			    f = fopen( (filenameNoExtension+"_q" + ss.str() + ".txt").c_str(), "w" );
			    for( int tidx=0; ok && nullptr!=f && tidx<m_vals.size(); ++tidx ){
				    for( int pidx=0; pidx<m_vals[tidx][qidx].size(); ++pidx ){
					    fprintf(f, "%1.16e, %1.16e, %1.16e, %1.16e, %1.16e", m_vals[tidx][qidx][pidx].I(), m_vals[tidx][qidx][pidx].Q(), m_vals[tidx][qidx][pidx].U(), m_vals[tidx][qidx][pidx].V(), m_weights[tidx][qidx][pidx] );
					    if(tidx<m_vals.size() || pidx<(m_vals[tidx][qidx].size()-1)) fprintf(f,"\n");
				    }
			    }
			    fclose(f);
			    fprintf(qf, "%1.16e", m_quadPtDistances[qidx]);
			    if(qidx<(m_vals[0].size()-1)) fprintf(qf,"\n");
		    }
		    fclose(qf);
        } else{
            nxLog::Record(NXLOG_WARNING, "SKTRAN_PhotonLog_PhotonsOnLos::PrintKernel, Couldn't open file %s.", (filenameNoExtension+"_distances.txt").c_str() );
        }
	}

	return ok;
}








bool SKTRAN_PhotonLog_ScatterPtOnLos::ConfigureObserverGeometry   ( const SKTRAN_RayOptical_Base& observerRay, double earthRadius )
{
	bool ok = true;

	WipeKernel();

	return ok;
}

bool SKTRAN_PhotonLog_ScatterPtOnLos::ConfigureKernel ( double minAlt, size_t numAlts, double deltaAlt, double minCosAngle, size_t numCosAngles, double deltaCosAngle, size_t numThreads )
{
	bool ok = true;

	m_units      .resize ( numThreads );
	m_incoming   .resize ( numThreads );
	m_scattered  .resize ( numThreads );
	
	WipeKernel();

	return ok;
}

void SKTRAN_PhotonLog_ScatterPtOnLos::WipeKernel ( )
{
	for(int tidx=0; tidx<m_units.size(); ++tidx){
		m_units      [tidx] .clear();
		m_incoming   [tidx] .clear();
		m_scattered  [tidx] .clear();
	}
}

bool SKTRAN_PhotonLog_ScatterPtOnLos::AddToKernel ( const SKTRAN_MCPhoton_Base* photon, size_t order, size_t threadid)
				//ok = ok && m_aveKernel->AddToKernel(mcphoton->m_scatterVector, mcphoton->photonRadiance().GetRecentContribVec(), order-1, mcphoton->m_isGroundScatter, threadid);
{
	bool ok = true;

	if(2==order){
		m_units     [threadid].emplace_back( photon->photonOptical()->LookVector() );
		m_incoming  [threadid].emplace_back( photon->photonRadiance().GetDebugVector() );
		m_scattered [threadid].emplace_back( photon->photonRadiance().GetRecentContribVec() );
	} 

	return ok;
}

bool SKTRAN_PhotonLog_ScatterPtOnLos::PrintKernel ( std::string filenameNoExtension )
{
	bool ok = true;

	{
		FILE * f;
		f = fopen( (filenameNoExtension+"_sp.txt").c_str(), "w" );
		for( int tidx=0; ok && nullptr!=f && tidx<m_units.size(); ++tidx ){
			for( int uidx=0; uidx<m_units[tidx].size(); ++uidx ){
				fprintf(f, "%1.6e, %1.6e, %1.6e, %1.6e, %1.6e, %1.6e, %1.6e, %1.6e, %1.6e, %1.6e, %1.6e%s", 
					m_units[tidx][uidx].X(), m_units[tidx][uidx].Y(), m_units[tidx][uidx].Z(),
					m_scattered[tidx][uidx].I(), m_scattered[tidx][uidx].Q(), m_scattered[tidx][uidx].U(), m_scattered[tidx][uidx].V(), 
					m_incoming [tidx][uidx].I(), m_incoming [tidx][uidx].Q(), m_incoming [tidx][uidx].U(), m_incoming [tidx][uidx].V(),
					(tidx<m_units.size() || uidx<(m_units[tidx].size()-1)) ? "\n" : "" );
				//if(tidx<m_units.size() || uidx<(m_units[tidx].size()-1)) fprintf(f,"\n");
			}
		}
		fclose(f);
	}

	return ok;
}








SKTRAN_PhotonLog_StDev::SKTRAN_PhotonLog_StDev ( )
{
	m_historyInterval = 0;
	m_numIntervals    = 0;
}
SKTRAN_PhotonLog_StDev::~SKTRAN_PhotonLog_StDev ( )
{
}
bool SKTRAN_PhotonLog_StDev::ConfigureObserverGeometry ( const SKTRAN_RayOptical_Base& observerRay, double earthRadius ) 
{
	return true;
}

bool SKTRAN_PhotonLog_StDev::ConfigureKernel ( double minAlt, size_t numAlts, double deltaAlt, double minCosAngle, size_t numCosAngles, double deltaCosAngle, size_t numThreads )
{
	m_sumRad_1fo.resize ( numThreads, 0.0 );
	m_sumRad_1so.resize ( numThreads, 0.0 );
	m_sumRad_1ho.resize ( numThreads, 0.0 );
	m_sumRad_2fo.resize ( numThreads, 0.0 );
	m_sumRad_2so.resize ( numThreads, 0.0 );
	m_sumRad_2ho.resize ( numThreads, 0.0 );
	m_sumRad_2ho_intermediate.resize( m_sumRad_1fo.size(), 0.0 );
	m_numPhotons.resize ( numThreads,   0 );
	m_histidx   .resize ( numThreads,   0 );

	m_history_fo.resize    ( numThreads );
	m_history_so.resize    ( numThreads );
	m_history_ho.resize    ( numThreads );
	ConfigureHistoryIntervals ( m_historyInterval, m_numIntervals );

	return true;
}

void SKTRAN_PhotonLog_StDev::WipeKernel ( ) 
{
	m_sumRad_1fo.resize ( m_sumRad_1fo.size(), 0.0 );
	m_sumRad_1so.resize ( m_sumRad_1fo.size(), 0.0 );
	m_sumRad_1ho.resize ( m_sumRad_1fo.size(), 0.0 );
	m_sumRad_2fo.resize ( m_sumRad_1fo.size(), 0.0 );
	m_sumRad_2so.resize ( m_sumRad_1fo.size(), 0.0 );
	m_sumRad_2ho.resize ( m_sumRad_1fo.size(), 0.0 );
	m_sumRad_2ho_intermediate.resize( m_sumRad_1fo.size(), 0.0 );
	m_numPhotons.resize ( m_sumRad_1fo.size(),   0 );
	m_histidx   .resize ( m_sumRad_1fo.size(),   0 );
	ConfigureHistoryIntervals ( m_historyInterval, m_numIntervals );
}

bool SKTRAN_PhotonLog_StDev::AddToKernel ( const SKTRAN_MCPhoton_Base* photon, size_t order, size_t threadid)
				//ok = ok && m_aveKernel->AddToKernel(mcphoton->m_scatterVector, mcphoton->photonRadiance().GetRecentContribVec(), order-1, mcphoton->m_isGroundScatter, threadid);
{
	const SKTRAN_Stokes_NC& thisScattContrib = photon->photonRadiance().GetRecentContribVec();

	switch(order){
	case 1:
		CheckForHistoryPush( threadid );
		m_sumRad_1fo[threadid] += thisScattContrib.I();
		m_sumRad_2fo[threadid] += thisScattContrib.I()*thisScattContrib.I();
		m_sumRad_2ho[threadid] += m_sumRad_2ho_intermediate[threadid]*m_sumRad_2ho_intermediate[threadid];
		m_sumRad_2ho_intermediate[threadid] = 0.0;
		++(m_numPhotons[threadid]);
		break;
	case 2:
		m_sumRad_1so[threadid] += thisScattContrib.I();
		m_sumRad_2so[threadid] += thisScattContrib.I()*thisScattContrib.I();
		break;
	default:
		m_sumRad_1ho[threadid] += thisScattContrib.I();
		m_sumRad_2ho_intermediate[threadid] += thisScattContrib.I();
		break;
	}
	return true;
}

void SKTRAN_PhotonLog_StDev::CheckForHistoryPush( size_t threadid )
{
	if(0 == (m_numPhotons[threadid] % m_historyInterval) && 0!=m_numPhotons[threadid] && m_histidx[threadid]<m_numIntervals ){
		m_history_fo[threadid][m_histidx[threadid]] = m_sumRad_1fo[threadid] / m_numPhotons[threadid];
		m_history_so[threadid][m_histidx[threadid]] =  m_sumRad_1so[threadid] / m_numPhotons[threadid];
		m_history_ho[threadid][m_histidx[threadid]] =  m_sumRad_1ho[threadid] / m_numPhotons[threadid];
		
		++m_histidx[threadid];
	}

	return;
}

bool SKTRAN_PhotonLog_StDev::PrintKernel ( std::string filenameNoExtension ) 
{
	bool ok = true;

	// Print variance traces to file
	{
	// First order
	FILE * f;
	f = fopen( (filenameNoExtension+"_fo.txt").c_str(), "w" );
	ok = ok && nullptr!=f;
	if(ok){
		for( size_t tidx=0; tidx<m_history_fo.size(); ++tidx){
			size_t hidx=0;
			for( ; hidx<m_history_fo[tidx].size(); ++hidx){
				fprintf(f, "%1.6e\t", m_history_fo[tidx][hidx]);
			}
			for( ; hidx<m_numIntervals; ++hidx){
				fprintf(f, "%1.6e\t", 0.0 );
			}
			fprintf(f, "\n");
		}
	}
	fclose(f);
	}
	{
	// Second order
	FILE * f;
	f = fopen( (filenameNoExtension+"_so.txt").c_str(), "w" );
	ok = ok && nullptr!=f;
	if(ok){
		for( size_t tidx=0; tidx<m_history_so.size(); ++tidx){
			size_t hidx=0;
			for( ; hidx<m_history_so[tidx].size(); ++hidx){
				fprintf(f, "%1.6e\t", m_history_so[tidx][hidx]);
			}
			for( ; hidx<m_numIntervals; ++hidx){
				fprintf(f, "%1.6e\t", 0.0 );
			}
			fprintf(f, "\n");
		}
	}
	fclose(f);
	}
	{
	// Higher order
	FILE * f;
	f = fopen( (filenameNoExtension+"_ho.txt").c_str(), "w" );
	ok = ok && nullptr!=f;
	if(ok){
		for( size_t tidx=0; tidx<m_history_ho.size(); ++tidx){
			size_t hidx=0;
			for( ; hidx<m_history_ho[tidx].size(); ++hidx){
				fprintf(f, "%1.6e\t", m_history_ho[tidx][hidx]);
			}
			for( ; hidx<m_numIntervals; ++hidx){
				fprintf(f, "%1.6e\t", 0.0 );
			}
			fprintf(f, "\n");
		}
	}
	fclose(f);
	}
	if(!ok) printf("Couldn't print StDev trace to files %s...", filenameNoExtension.c_str());

	return ok;
}

bool SKTRAN_PhotonLog_StDev::ConfigureHistoryIntervals ( size_t interval, size_t numIntervals )
{
	bool ok = true;

	m_historyInterval = interval;
	m_numIntervals = numIntervals;
	
	std::for_each( m_history_fo.begin(), m_history_fo.end(), [&]( decltype(m_history_fo.front())& h ){ h.resize(m_numIntervals,0.0); } );
	std::for_each( m_history_so.begin(), m_history_so.end(), [&]( decltype(m_history_so.front())& h ){ h.resize(m_numIntervals,0.0); } );
	std::for_each( m_history_ho.begin(), m_history_ho.end(), [&]( decltype(m_history_ho.front())& h ){ h.resize(m_numIntervals,0.0); } );

	return ok;
}

SKTRAN_MCThreadRadianceLogger::SKTRAN_MCThreadRadianceLogger( )
{
	m_minNumSamplesHigherOrder = 100;
	SetMaxOrder(0);
	SetMinFractionHigherOrder(0.1);
	m_parent = nullptr;
	m_primaryWavelengthIndex = 0;
	m_runningSums.resize(1);
	m_chunkCounter = 1;
	m_chunkSize = 1;

	ResetLog();
}

SKTRAN_MCThreadRadianceLogger::SKTRAN_MCThreadRadianceLogger( SKTRAN_MCThreadRadianceLogger& other )
{
	m_minNumSamplesHigherOrder = other.m_minNumSamplesHigherOrder;
	SetMaxOrder ( other.m_hardMaxOrder );
	SetMinFractionHigherOrder( other.m_minFracHO );
	m_parent    = &other;
	m_primaryWavelengthIndex = other.m_primaryWavelengthIndex;
	m_runningSums = other.m_runningSums;
	m_scatterManager = other.m_scatterManager;
	m_chunkCounter = other.m_chunkCounter;
	m_chunkSize = other.m_chunkSize;
	ResetLog( );
}

SKTRAN_MCThreadRadianceLogger::~SKTRAN_MCThreadRadianceLogger ( )
{
	if(nullptr!=m_parent){
		// Merge back into parent
		m_parent->Merge( *this );
	}
}

bool SKTRAN_MCThreadRadianceLogger::SecondaryMeasurementExists() const
{
	return m_scatterManager->SecondaryMeasurement();
}

bool SKTRAN_MCThreadRadianceLogger::ConfigureWavelengths(const std::vector<double>& wavelengths, double currentWavelength, size_t primaryWavelengthIndex)
{
	bool ok = true;
	size_t numwl = wavelengths.size();
	if (numwl > 0)
	{
		m_primaryWavelengthIndex = primaryWavelengthIndex;
		m_runningSums.resize(numwl);
		for (size_t i = 0; i < numwl; i++) m_runningSums[i].wavelength = wavelengths[i];
		ok = ok && primaryWavelengthIndex < numwl;
	}
	else
	{
		m_primaryWavelengthIndex = 0;
		m_runningSums.resize(1);
		m_runningSums[0].wavelength = currentWavelength;
	}

	ResetLog();
	return ok;
}

bool SKTRAN_MCThreadRadianceLogger::ConfigureOptimalScatterSequence(const SKTRAN_OptimalScatterSequenceManager_Base* scatterManager)
{
	bool ok = true;
	ok = ok && scatterManager != nullptr; 
	if (ok)
	{
		m_scatterManager = scatterManager;
		size_t i = 0;
		for (auto&& rSums : m_runningSums)
		{
			m_scatterManager->ConfigureRunningSums(rSums);
			rSums.primary = m_primaryWavelengthIndex == i++;
		}
	}
	ResetLog();
	return ok;
}

void SKTRAN_MCThreadRadianceLogger::ResetLog ( )
{
	for (auto&& wlinfo : m_runningSums)
	{
		for( int idx=0; idx<wlinfo.m_temprad.size(); ++idx){
		//	wlinfo.m_temprad[idx].SetTo(0.0);
		//	wlinfo.m_oinfo[idx].SetToZero();
			wlinfo.radSum[idx].SetTo(0.0);
			wlinfo.radBuffer[idx].SetTo(0.0);
		}
		//for( int idx=0; idx<wlinfo.m_cov.size(); ++idx){
		//	wlinfo.m_cov[idx] = 0.0;
		//}

		std::fill(wlinfo.rad2SumVar.begin(), wlinfo.rad2SumVar.end(), 0.0);
		std::fill(wlinfo.rad2SumCov.begin(), wlinfo.rad2SumCov.end(), 0.0);
		std::fill(wlinfo.numSamples.begin(), wlinfo.numSamples.end(), 0);
	}
}

void SKTRAN_MCThreadRadianceLogger::orderInfo::SetToZero( )
{
	hist1  .SetTo( 0.0 );
	hist2        = 0.0;
	variance     = 0.0;
	numSamples   =   0;
}

bool SKTRAN_MCThreadRadianceLogger::Submit( size_t ossIdx, size_t order, const SKTRAN_MCPhoton_Base* photon )
{
	return m_scatterManager->SubmitSample(ossIdx, order, photon, m_runningSums);
}

bool SKTRAN_MCThreadRadianceLogger::DeclareRayDone( size_t scatterSequenceIndex, SKTRAN_MCVarianceLogger& varianceLogger, size_t threadid)
{
	bool ok = true;

	bool reset = m_runningSums[m_primaryWavelengthIndex].minSamplesComplete && ++m_chunkCounter > m_chunkSize;
	if (reset) m_chunkCounter = 1;

	for (auto&& rSums : m_runningSums)
	{
		ok = ok && m_scatterManager->SortSamples(scatterSequenceIndex, rSums);
	}

	for (auto&& rSums : m_runningSums)
	{
		if (reset)
			ok = ok && m_scatterManager->ProcessVariance(rSums);
	}

	if (reset)
	{
		double variance = 0.0;
		size_t numSamples = 0;
		ok = ok && m_scatterManager->CalculateTargetVariance(m_runningSums[m_primaryWavelengthIndex], variance, numSamples);		
		varianceLogger.UpdateThreadVariance(variance, numSamples, threadid);
	}

	return ok;
}

bool SKTRAN_MCThreadRadianceLogger::DiscardRay()
{
	bool ok = true;

	for (auto&& rSums : m_runningSums)
	{
		for (auto&& rad : rSums.radBuffer)
		{
			rad.SetTo(0.0);
		}
	}

	return ok;
}

void SKTRAN_MCThreadRadianceLogger::SetMaxOrder( size_t hardmax )
{
	m_hardMaxOrder = hardmax;
	m_oinfoMaxIndex = min(MC_NUMDISTINCTORDERS, m_hardMaxOrder) - 1;
}

void SKTRAN_MCThreadRadianceLogger::orderInfo::UpdateVar ( )
{
	variance = 0<numSamples? (hist2 - (pow(hist1.I(),2.0)/((double)numSamples))) * pow(((double)numSamples), -2.0) : 0.0;
}

bool SKTRAN_MCThreadRadianceLogger::SetMinFractionHigherOrder ( double minfrac )
{
	return SetMinFractionHigherOrder(std::vector<double>(1, minfrac));
}

bool SKTRAN_MCThreadRadianceLogger::SetMinFractionHigherOrder(const std::vector<double>& minfrac)
{
	bool ok = true;

	for (auto&& frac : minfrac) ok = ok && (0.0 <= frac) && (frac <= 1.0);
	for (auto&& wlinfo : m_runningSums)
	{
		for (size_t idx = 0; idx < wlinfo.numSamples.size(); ++idx)
		{
			ok = ok && wlinfo.numSamples[idx] == 0; // ok = ok && 0 == wlinfo.m_oinfo[idx].numSamples;
		}
	}

	if (ok) {
		m_minFracHO = minfrac;
	}
	else {
		// nxLog::Record(NXLOG_WARNING, "SKTRAN_MCThreadRadianceLogger::SetMinFractionHigherOrder, Could not set minfrac %f.", minfrac);
	}

	return ok;
}

//size_t SKTRAN_MCThreadRadianceLogger::OptimalMaxOrderScatter( double r )
//{
//	size_t optorder;
//
//	double target;
//	const std::vector<orderInfo>& oinfo = m_runningSums[m_primaryWavelengthIndex].m_oinfo;
//	const std::vector<double>& cov = m_runningSums[m_primaryWavelengthIndex].m_cov;
//
//	if( oinfo[ m_oinfoMaxIndex ].numSamples<m_minNumSamplesHigherOrder || r<m_minFracHO[0] )  // If we haven't sampled enough of the ho rays yet OR we want to make sure we sample the ho space by at least some min amount
//	{
//		optorder = m_hardMaxOrder; // Scatter to higher order
//	} else{
//		r = (r-m_minFracHO[0])/(1.0-m_minFracHO[0]); // Rescale uni(0,1) to account for using it twice
//		double numhscattersave = 20.0; // Assuming 20 orders of scatter for h
//		
//
//		{ 
//		std::array<double,MC_NUMDISTINCTORDERS-1> var;
//		std::array<double,MC_NUMDISTINCTORDERS-1> pro;
//		double runningVar;
//		double runningCov;
//		size_t covidx = 7; 
//		for(size_t idx=1; idx<MC_NUMDISTINCTORDERS; ++idx){
//			runningVar = oinfo[idx].variance;
//			runningCov = 0.0;
//			for(size_t cidx=idx+1; cidx<MC_NUMDISTINCTORDERS && 0<oinfo[cidx].numSamples; ++cidx, ++covidx){
//				//runningCov = runningCov + m_cov[covidx];
//				runningCov += 2.0 * ( cov[covidx] - ( oinfo[idx].hist1.I()*(1.0/oinfo[idx].numSamples) * oinfo[cidx].hist1.I()) ) * pow((double)oinfo[cidx].numSamples,-1.0) *  pow((double)oinfo[idx].numSamples,-1.0);
//			}
//			var[idx-1] = sqrt(max(runningVar+runningCov, 0.0)); // Assign to variance vector
//		}
//		//var.back() = sqrt(m_oinfo[MC_NUMDISTINCTORDERS-1].variance + (runningCov*m_oinfo[MC_NUMDISTINCTORDERS-1].variance/m_oinfo[MC_NUMDISTINCTORDERS-2].variance) );	// Estimate cov amongst higher orders as prop to var
//		for(size_t idx=0;idx<MC_NUMDISTINCTORDERS-2; ++idx){
//			pro[idx] = var[idx] - var[idx+1];
//		}
//		pro.back() = var[MC_NUMDISTINCTORDERS-2];
//		double sum = 0.0;
//		for(size_t idx=0; idx<pro.size(); ++idx){
//			sum += pow(pro[idx],1.0);//* ((double)(idx+2));
//			pro[idx] = sum;
//		}
//		target = r* *(pro.end()-1);
//		auto it = std::upper_bound( pro.begin(), pro.end(), target );
//		if( it< pro.end()-1 ){
//			optorder = std::distance(pro.begin(),it)+2;
//		}else{
//			optorder = m_hardMaxOrder;
//		}
//		//printf("%i\n",optorder);
//		}
//	}
//	
//	for (size_t idx = 0; idx < min(optorder, MC_NUMDISTINCTORDERS); ++idx)
//	{
//		for (auto&& wlinfo : m_runningSums) ++(wlinfo.m_oinfo[idx].numSamples);
//	}
//
//	return optorder;
//}

size_t SKTRAN_MCThreadRadianceLogger::FindIdx( size_t order ) const {
	return order<MC_NUMDISTINCTORDERS ? order-1 : MC_NUMDISTINCTORDERS-1;
}

void SKTRAN_MCThreadRadianceLogger::Merge( const SKTRAN_MCThreadRadianceLogger& other )
{
	#pragma omp critical (__SKTRAN_MCThreadRadianceLogger_Merge_ompCritLabel)
	{
		auto otherwlinfo = other.m_runningSums.cbegin();
		for (auto wlinfo = m_runningSums.begin(); wlinfo != m_runningSums.end(); wlinfo++, otherwlinfo++)
		{
			for (size_t idx = 0; idx < m_scatterManager->NumVarianceTerms(); ++idx)
			{
				wlinfo->numSamples[idx] += otherwlinfo->numSamples[idx];
				wlinfo->radSum[idx] += otherwlinfo->radSum[idx];
				wlinfo->rad2SumVar[idx] += otherwlinfo->rad2SumVar[idx];
			}
			for (size_t idx = 0; idx < m_scatterManager->NumCovarianceTerms(); ++idx) {
				wlinfo->rad2SumCov[idx] += otherwlinfo->rad2SumCov[idx];
			}

			m_scatterManager->ProcessVariance(*wlinfo);
		}
	}
}

void SKTRAN_MCThreadRadianceLogger::orderInfo::AddToMe ( const orderInfo& other )
{
	numSamples += other.numSamples;
	hist1    += other.hist1;
	hist2    += other.hist2;
	UpdateVar( );
}

bool SKTRAN_MCThreadRadianceLogger::OptimalScatterSequenceIndex(double randNum, size_t & ossIdx, size_t & maxOrder)
{
	bool ok = true;
	ok = ok && m_scatterManager->OptimalScatterSequenceIndex(m_runningSums[m_primaryWavelengthIndex], randNum, ossIdx, m_runningSums[m_primaryWavelengthIndex].minSamplesComplete);
	ok = ok && m_scatterManager->Order(ossIdx, maxOrder);
	return ok;
}

SKTRAN_Stokes_NC SKTRAN_MCThreadRadianceLogger::TotalMeasurement ( ) const
{
	return TotalMeasurement(m_primaryWavelengthIndex);
}

SKTRAN_Stokes_NC SKTRAN_MCThreadRadianceLogger::TotalMeasurement(size_t wlidx) const
{	
	//SKTRAN_Stokes_NC result;
	//skRTStokesVector::SetToZero(result);

	////for (size_t idx = 0; idx < MC_NUMDISTINCTORDERS; ++idx)
	////	if (0 < m_runningSums[wlidx].m_oinfo[idx].numSamples)
	////		result += m_runningSums[wlidx].m_oinfo[idx].hist1 * (1.0 / m_runningSums[wlidx].m_oinfo[idx].numSamples);

	//for (size_t idx = 0; idx < m_scatterManager->NumVarianceTerms(); idx++)
	//{
	//	if (0 < m_runningSums[wlidx].numSamples[idx])
	//	{
	//		result += m_runningSums[wlidx].radSum[idx] * (1.0 / m_runningSums[wlidx].numSamples[idx]);
	//	}
	//}

	SKTRAN_Stokes_NC result;
	bool ok = m_scatterManager->CalculateMeasurement(m_runningSums[wlidx], result);

	return result;
}

double SKTRAN_MCThreadRadianceLogger::TargetMeasurement() const
{
	return TargetMeasurement(m_primaryWavelengthIndex);
}

double SKTRAN_MCThreadRadianceLogger::TargetMeasurement(size_t wlidx) const
{
	double measurement;
	bool ok = m_scatterManager->CalculateTargetMeasurement(m_runningSums[wlidx], measurement);
	return measurement;
}

SKTRAN_Stokes_NC SKTRAN_MCThreadRadianceLogger::MeasurementAtScatterIndex( size_t scidx ) const
{
	return MeasurementAtScatterIndex(m_primaryWavelengthIndex, scidx);
}

SKTRAN_Stokes_NC SKTRAN_MCThreadRadianceLogger::MeasurementAtScatterIndex(size_t scidx, size_t wlidx) const
{
	SKTRAN_Stokes_NC result;
	skRTStokesVector::SetToZero(result);

	size_t idx = FindIdx(scidx);
	if (idx >= MC_NUMDISTINCTORDERS)    nxLog::Record(NXLOG_WARNING, "SKTRAN_MCThreadRadianceLogger::Measurement, Returning aggregate higher order data.");
	if (0 < m_runningSums[wlidx].m_oinfo[idx].numSamples) result += m_runningSums[wlidx].m_oinfo[idx].hist1 * (1.0 / m_runningSums[wlidx].m_oinfo[idx].numSamples);
	return result;
}

bool SKTRAN_MCThreadRadianceLogger::ExportStatistics(size_t losIdx, double wl)
{
	char buffer[50];
	int n;

	if (m_runningSums.size() > 1) 		
		n = sprintf(buffer, "-los%03d", losIdx);
	else 
		n = sprintf(buffer, "-los%03d-%.2fnm", losIdx, wl);

	return m_scatterManager->ExportStatistics(m_runningSums, buffer);
}

double SKTRAN_MCThreadRadianceLogger::Variance ( ) const
{
	return Variance(m_primaryWavelengthIndex);
}

double SKTRAN_MCThreadRadianceLogger::Variance(size_t wlidx) const
{
	double result = 0.0;

	bool ok = m_scatterManager->CalculateTotalVariance(m_runningSums[wlidx], result);

	return result;
}

double SKTRAN_MCThreadRadianceLogger::SecondaryMeasurement() const
{
	return SecondaryMeasurement(m_primaryWavelengthIndex);
}

double SKTRAN_MCThreadRadianceLogger::SecondaryMeasurement(size_t wlidx) const
{
	double result = 0.0;

	bool ok = m_scatterManager->CalculateSecondaryMeasurement(m_runningSums[wlidx], result);

	return result;
}

double SKTRAN_MCThreadRadianceLogger::SecondaryVariance() const
{
	return SecondaryVariance(m_primaryWavelengthIndex);
}

double SKTRAN_MCThreadRadianceLogger::SecondaryVariance(size_t wlidx) const
{
	double result = 0.0;

	bool ok = m_scatterManager->CalculateSecondaryVariance(m_runningSums[wlidx], result);

	return result;
}

SKTRAN_MCAirMassFactorLogger::SKTRAN_MCAirMassFactorLogger()
{
	m_parent = nullptr;
}

SKTRAN_MCAirMassFactorLogger::SKTRAN_MCAirMassFactorLogger( SKTRAN_MCAirMassFactorLogger& other )
{
	m_parent = &other;
	Initialize(other);
}

SKTRAN_MCAirMassFactorLogger::~SKTRAN_MCAirMassFactorLogger()
{
	if (nullptr != m_parent) {
		// Merge back into parent
		m_parent->Merge(*this);
	}
}

void SKTRAN_MCAirMassFactorLogger::Merge(const SKTRAN_MCAirMassFactorLogger& other)
{
	#pragma omp critical (__SKTRAN_MCAirMassFactorLogger_Merge_ompCritLabel)
	{
		if (m_numamfcells > 0) {
			m_amfSumsNumSamples += other.m_amfSumsNumSamples;
			m_amfSumsI += other.m_amfSumsI;
			m_amfSumsII += other.m_amfSumsII;
			for (size_t amfidx = 0; amfidx < m_numamfcells; amfidx++) {
				m_amfSumsW[amfidx] += other.m_amfSumsW[amfidx];
				m_amfSumsWI[amfidx] += other.m_amfSumsWI[amfidx];
				m_amfSumsWW[amfidx] += other.m_amfSumsWW[amfidx];
			}
		}
	}
}

void SKTRAN_MCAirMassFactorLogger::Submit( int order, const SKTRAN_MCPhoton_Base* photon )
{
	m_temprad += photon->photonRadiance().GetRecentContribVec();
	for (size_t amfidx = 0; amfidx < m_numamfcells; amfidx++) {
		m_tempwsc[amfidx] += (photon->m_solarSlantColumns[amfidx] + photon->m_scatterSlantColumns[amfidx]) * photon->photonRadiance().GetRecentContribSca();
	}
}

void SKTRAN_MCAirMassFactorLogger::DeclareRayDone()
{
	double radiance = m_temprad.I();
	m_amfSumsNumSamples += 1;
	m_amfSumsI += radiance;
	m_amfSumsII += radiance * radiance;
	for (size_t amfidx = 0; amfidx < m_numamfcells; amfidx++)
	{
		m_amfSumsW[amfidx] += m_tempwsc[amfidx];
		m_amfSumsWI[amfidx] += m_tempwsc[amfidx] * radiance;
		m_amfSumsWW[amfidx] += m_tempwsc[amfidx] * m_tempwsc[amfidx];
	}
	m_temprad.SetTo(0.0);
	std::fill(m_tempwsc.begin(), m_tempwsc.end(), 0.0);
}

void SKTRAN_MCAirMassFactorLogger::DiscardRay()
{
	m_temprad.SetTo(0.0);
	std::fill(m_tempwsc.begin(), m_tempwsc.end(), 0.0);
}

void SKTRAN_MCAirMassFactorLogger::Initialize(size_t numcells, const std::vector<double>& verticalColumns)
{
	m_numamfcells = numcells;
	m_temprad.SetTo(0.0);
	m_tempwsc.resize(m_numamfcells, 0.0);
	m_amfCellColumns.resize(m_numamfcells, 0.0);
	m_amfSumsNumSamples = 0;
	m_amfSumsI = 0.0;
	m_amfSumsII = 0.0;
	m_amfSumsW.resize(m_numamfcells, 0.0);
	m_amfSumsWI.resize(m_numamfcells, 0.0);
	m_amfSumsWW.resize(m_numamfcells, 0.0);

	std::copy(verticalColumns.begin(), verticalColumns.end(), m_amfCellColumns.begin());

	std::fill(m_tempwsc.begin(), m_tempwsc.end(), 0.0);
	std::fill(m_amfSumsW.begin(), m_amfSumsW.end(), 0.0);
	std::fill(m_amfSumsWI.begin(), m_amfSumsWI.end(), 0.0);
	std::fill(m_amfSumsWW.begin(), m_amfSumsWW.end(), 0.0);
}

void SKTRAN_MCAirMassFactorLogger::Initialize(SKTRAN_MCAirMassFactorLogger& other) {
	m_numamfcells = other.m_numamfcells;
	m_temprad.SetTo(0.0);
	m_tempwsc.resize(m_numamfcells, 0.0);
	m_amfSumsNumSamples = 0;
	m_amfSumsI = 0.0;
	m_amfSumsII = 0.0;
	m_amfSumsW.resize(m_numamfcells, 0.0);
	m_amfSumsWI.resize(m_numamfcells, 0.0);
	m_amfSumsWW.resize(m_numamfcells, 0.0);

	m_amfCellColumns = other.m_amfCellColumns;

	std::fill(m_amfSumsW.begin(), m_amfSumsW.end(), 0.0);
	std::fill(m_amfSumsWI.begin(), m_amfSumsWI.end(), 0.0);
	std::fill(m_amfSumsWW.begin(), m_amfSumsWW.end(), 0.0);
}

double SKTRAN_MCAirMassFactorLogger::AirMassFactor(size_t amfidx) const
{
	if (amfidx < m_numamfcells) {
		double              N = (double)m_amfSumsNumSamples;

		// calculate means
		double meanI = m_amfSumsI / N;
		double meanW = m_amfSumsW[amfidx] / N;

		// calculate the mean (a 2nd order correction is commented out)
		return meanW / (meanI * m_amfCellColumns[amfidx]);// -(covWI[amfidx] * std::pow(meanI, -2.0)) + (varI * meanW[amfidx] * std::pow(meanI, -3.0));
	}
	else
	{
		nxLog::Record(NXLOG_ERROR, "SKTRAN_MCAirMassFactorLogger::AirMassFactor, Index %d is out of bounds.", amfidx);
		return 0.0;
	}
}

double SKTRAN_MCAirMassFactorLogger::AirMassFactorVariance(size_t amfidx) const
{
	if (amfidx < m_numamfcells)
	{
		double N = (double)m_amfSumsNumSamples;

		// calculate means, variances and covariances
		double meanI = m_amfSumsI / N;
		double varI = (N*m_amfSumsII - m_amfSumsI * m_amfSumsI) * std::pow(N, -3);
		double meanW = m_amfSumsW[amfidx] / N;
		double varW = (N*m_amfSumsWW[amfidx] - std::pow(m_amfSumsW[amfidx], 2.0)) * std::pow(N, -3);
		double covWI = (N*m_amfSumsWI[amfidx] - m_amfSumsW[amfidx] * m_amfSumsI) * std::pow(N, -3);
	
		// 2nd order approximation (there is no exact expression for the variance of the ratio of two normal distributions)
		return (varW * std::pow(meanI, -2) - 2.0*covWI * meanW * std::pow(meanI, -3) + varI * meanW * meanW * std::pow(meanI, -4)) * std::pow(m_amfCellColumns[amfidx], -2);
	}
	else
	{
		nxLog::Record(NXLOG_ERROR, "SKTRAN_MCAirMassFactorLogger::AirMassFactorVariance, Index %d is out of bounds.", amfidx);
		return 0.0;
	}
}


SKTRAN_MCConvergenceDetector_ThreadIsolated::SKTRAN_MCConvergenceDetector_ThreadIsolated() : 
	m_radiance(0), m_variance(0), m_numSamples(0) {;}

bool SKTRAN_MCConvergenceDetector_ThreadIsolated::SetNumThreads( size_t numthreads ){
	return true;
}

void SKTRAN_MCConvergenceDetector_ThreadIsolated::UpdateThreadData( const SKTRAN_MCThreadRadianceLogger& radlog ){
	m_radiance   = radlog.TotalMeasurement().I(); 
	m_variance   = radlog.Variance();
	m_numSamples = radlog.n(1);
}

bool SKTRAN_MCConvergenceDetector_ThreadIsolated::IsConvergenceReached( ){
	return !( (m_radiance*m_radiance*GetUserDesiredStdev()*GetUserDesiredStdev())<m_variance || m_numSamples<GetMinNumRaysPerThread() );
}


//
//SKTRAN_MCConvergenceDetector_ThreadCrosstalk::SKTRAN_MCConvergenceDetector_ThreadCrosstalk()
//{
//	m_trueUpdateThreadData = nullptr;
//}
//
//SKTRAN_MCConvergenceDetector_ThreadCrosstalk::SKTRAN_MCConvergenceDetector_ThreadCrosstalk( SKTRAN_MCConvergenceDetector_ThreadCrosstalk& parent )
//{
//
//}

