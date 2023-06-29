#include "include/sktran_hr_internals.h"
#include <numeric>
#include <functional>
#include <boost/multi_array.hpp>
#include <iostream>
using boost::multi_array;
using boost::extents;

RadStore_Polarized::RadStore_Polarized (){
	ConfigurePolarizedScattering ( 1, false );
}


void RadStore_Polarized::ConfigurePolarizedScattering ( size_t numPolarizedDiffuseOrders, bool mode_c, double adaptiveFraction ){
	// PV<numPolarizedDiffusions+2>

	// Store the polarization settings because we have to reset them for the adaptive orders calculation
	m_numDiffusePolarizedOrder = numPolarizedDiffuseOrders;
	m_modec = mode_c;
	m_adaptivefraction = adaptiveFraction;

	m_maxPolarizedOrder_scattering = 2<numPolarizedDiffuseOrders ? (numPolarizedDiffuseOrders-2) : 0;
	if( mode_c ){
		m_maxPolarizedOrder_integration = 9000;
	} else{
		m_maxPolarizedOrder_integration = 0<numPolarizedDiffuseOrders ? (numPolarizedDiffuseOrders) : 0;
	}
}


bool RadStore_Polarized::CleanDiffuseIndexes ( )
{
	SKTRAN_Stokes_NC zero;
	zero.SetTo( 0.0 );
	m_totalincomingvector .assign( m_totalincomingvector .size(), zero );
	m_incomingvector      .assign( m_incomingvector      .size(), zero );
	m_outgoingLinPols     .assign( m_outgoingLinPols     .size(), linPolComponents(0.0, 0.0) );

	// Reset the scattering order options
	ConfigurePolarizedScattering(m_numDiffusePolarizedOrder, m_modec, m_adaptivefraction);

	return RadStore_PV1::CleanDiffuseIndexes();
}


bool RadStore_Polarized::AllocateStorage( size_t incomingcounter, size_t outgoingcounter )
{
	bool ok = true;
	m_incomingvector      .resize( incomingcounter );
	m_totalincomingvector .resize( incomingcounter );	
	m_outgoingLinPols     .resize( outgoingcounter );
	m_outgoingLinPols_lastOrder     .resize( outgoingcounter );
	ok = ok && RadStore_PV1::AllocateStorage( incomingcounter, outgoingcounter );
	return ok;
}


bool RadStore_Polarized::StoreFirstOrderIncoming ( const SKTRAN_SourceTermIntegrator_Base* srcintegrator, const std::vector<SKTRAN_Source_Term*>* sources, const SKTRAN_RayOptical_Base* ray, size_t index )
{
	bool ok = true;

	SKTRAN_Stokes_NC radiance;
	ok = ok && srcintegrator->IntegrateSourceTerm( ray, radiance, *sources );
	m_totalincomingvector [ index ] = radiance; // Store only the incoming polarized radiance 
	m_incomingvector[index] = radiance;

	return ok;
}


bool RadStore_Polarized::GetPolarizedComponent_basisImplicit ( size_t radIndex, const nxVector& look_transformed, SKTRAN_Stokes_NC& incomingvec, const HELIODETIC_POINT& queryPoint, const HELIODETIC_UNITVECTOR& inrayQCoords ) const
{
	bool ok = true;
	incomingvec = m_totalincomingvector[radIndex]; // totalinvec is stored in the default basis for the incoming ray, HELIODETIC_BASIS(obs,look)
	return ok;
}


bool RadStore_Polarized::DumpIncomingRadiances ( const SKTRAN_HR_Diffuse_Table_CPU* table, size_t order, double wlen )
{
	double	x,y,z;
	double	radiance;
	char	dsetname[100];

	std::vector<size_t> indices = table->DiagnosticProfIdxs();

	size_t nalts = table->DiagnosticProfAlts();
	size_t npts = table->DiffusePointAt( indices[0] ).NumIncomingRays();
	size_t ndata = 10; //[x, y, z, rad, true_x, true_y, true_z, stokes_x, stokes_y, stokes_z]

	multi_array<double, 3> maindata( extents[nalts][npts][ndata] );

	hsize_t dimsf[3] = { nalts, npts, ndata };

	hid_t H5Fid, H5Sid, H5Did; // H5 id-s for file, dataset and space
	hid_t propList;

	H5Sid = H5Screate_simple( 3, dimsf, NULL );
	propList = H5Pcreate( H5P_DATASET_CREATE );
	H5Pset_layout( propList, H5D_CHUNKED );
	H5Pset_chunk( propList, 3, dimsf );

	H5Fid = H5Fopen( "DiagnosticData.h5", H5F_ACC_RDWR, H5P_DEFAULT );
	if( H5Fid < 0 )
	{
		nxLog::Record( NXLOG_WARNING, "Could not open h5 diagnostic file. Thats not good" );
		return false;
	}

	for( size_t idx = 0; idx < indices.size(); idx++ )
	{
		for( size_t pointidx = indices[idx]; pointidx < indices[idx] + nalts; pointidx++ )
		{
			for( size_t rayidx = 0; rayidx < table->DiffusePointAt( pointidx ).NumIncomingRays(); rayidx++ )
			{
				x = table->DiffusePointAt( pointidx ).IncomingUnitSphere()->UnitVectorAt( rayidx ).X();
				y = table->DiffusePointAt( pointidx ).IncomingUnitSphere()->UnitVectorAt( rayidx ).Y();
				z = table->DiffusePointAt( pointidx ).IncomingUnitSphere()->UnitVectorAt( rayidx ).Z();

				HELIODETIC_UNITVECTOR units[3];
				HELIODETIC_VECTOR trueunit;
				table->DiffusePointAt( pointidx ).Location().LocalUnitVectors( units, 3 );
				trueunit = HELIODETIC_VECTOR( units[0], x ) + HELIODETIC_VECTOR( units[1], y ) + HELIODETIC_VECTOR( units[2], z );
				radiance = DumpIncomingRadiances_GetIncomingRadiance( table->DiffusePointAt( pointidx ).IncomingRadianceIdx( rayidx ) );
				
				double Q = m_totalincomingvector[table->DiffusePointAt( pointidx ).IncomingRadianceIdx( rayidx )].Q();
				double U = m_totalincomingvector[table->DiffusePointAt( pointidx ).IncomingRadianceIdx( rayidx )].U();
				double V = m_totalincomingvector[table->DiffusePointAt( pointidx ).IncomingRadianceIdx( rayidx )].V();
				
				maindata[pointidx - indices[idx]][rayidx][0] = x;
				maindata[pointidx - indices[idx]][rayidx][1] = y;
				maindata[pointidx - indices[idx]][rayidx][2] = z;
				maindata[pointidx - indices[idx]][rayidx][3] = radiance;
				maindata[pointidx - indices[idx]][rayidx][4] = trueunit.X();
				maindata[pointidx - indices[idx]][rayidx][5] = trueunit.Y();
				maindata[pointidx - indices[idx]][rayidx][6] = trueunit.Z();
				maindata[pointidx - indices[idx]][rayidx][7] = Q;
				maindata[pointidx - indices[idx]][rayidx][8] = U;
				maindata[pointidx - indices[idx]][rayidx][9] = V;
			}			
		}
		// write dataset for this profile
		sprintf( dsetname, "in_wlen_%0.2f_ord_%d_prof_%d", wlen, int(order), int( indices[idx]/nalts ) );
		H5Did = H5Dcreate2( H5Fid, dsetname, H5T_NATIVE_DOUBLE, H5Sid, H5P_DEFAULT, propList, H5P_DEFAULT );
		H5Dwrite( H5Did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, maindata.data() );
		H5Dclose( H5Did );
	}
	H5Pclose( propList );
	H5Sclose( H5Sid );

	// and dump data for ground points
	size_t ngp = table->NumGroundPoints();
	size_t ngr = table->GroundPointAt( 0 ).NumIncomingRays();
	hsize_t dims2[2] = { ngp, ngr };

	H5Sid = H5Screate_simple( 2, dims2, NULL );
	propList = H5Pcreate( H5P_DATASET_CREATE );
	H5Pset_layout( propList, H5D_CHUNKED );
	H5Pset_chunk( propList, 2, dims2 );

	multi_array<double, 2> grounddata( extents[ngp][ngr] );

	for( size_t pointidx = 0; pointidx < table->NumGroundPoints(); pointidx++ )
	{
		for( size_t rayidx = 0; rayidx < table->GroundPointAt(pointidx).NumIncomingRays(); rayidx++ )
		{
			grounddata[pointidx][rayidx] = m_incomingscalar[table->GroundPointAt(pointidx).IncomingRadianceIdx(rayidx)];
		}
	}
	sprintf( dsetname, "in_wlen_%0.2f_ord_%d_ground", wlen, int(order) );
	H5Did = H5Dcreate2( H5Fid, dsetname, H5T_NATIVE_DOUBLE, H5Sid, H5P_DEFAULT, propList, H5P_DEFAULT );
	H5Dwrite( H5Did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, grounddata.data() );
	H5Dclose( H5Did );

	H5Pclose( propList );
	H5Sclose( H5Sid );
	H5Fclose( H5Fid );
	maindata.resize( extents[0][0][0] );
	grounddata.resize( extents[0][0] );
	return true;
}


bool RadStore_Polarized::ComputeNextOrderPoint_IntegrateVectors( const SKTRAN_HR_Diffuse_Point& point )
{
	bool ok = true;
	
	SKTRAN_Stokes_NC incoming, temp;
	size_t index;
	double scalarWeight; 
	double q, u;
	for( size_t inidx = 0; inidx < point.NumIncomingRays(); inidx++ )
	{
		incoming.SetTo(0.0);
		auto& diffIndices = point.DiffuseIndices()[inidx];
		for( size_t idx = 0; idx < point.DiffuseIndices()[inidx].diffindex.size(); idx++ )
		{
			index = diffIndices.diffindex[idx].index;
			scalarWeight = m_outgoingscalar[index] * diffIndices.diffindex[idx].weight;
			q = scalarWeight*m_outgoingLinPols[index].Q;
			u = scalarWeight*m_outgoingLinPols[index].U;
			temp.SetTo( scalarWeight, q, u, 0.0 );
			incoming += temp;
		}
		m_totalincomingvector [ point.IncomingRadianceIdx(inidx) ] += incoming;
		m_incomingscalar      [ point.IncomingRadianceIdx(inidx) ]  = SKTRAN_HR_WEIGHT_TYPE(incoming.I());
		m_incomingvector      [ point.IncomingRadianceIdx(inidx) ]  = incoming;

		// really small numbers hurt performance through denormalized floating point arithmetic
		if( m_incomingscalar[point.IncomingRadianceIdx(inidx)] < 1E-14 )
		{
			//m_incomingradiancevector[inidx + incomingidx] = 0.0;
		}
	}
	return ok;
}


bool RadStore_Polarized::ComputeNextOrderPoint_IntegrateScalars( const SKTRAN_HR_Diffuse_Point& point )
{
	return RadStore_PV1::ComputeNextOrderPoint( point );
}


bool RadStore_Polarized::ScatterPoint( const SKTRAN_HR_Diffuse_Point& point, const Avals_Base& Avals )
{
	bool ok = true;

	if( m_order_outgoing <= m_maxPolarizedOrder_scattering ){
		if( !point.IsGroundPoint() ){
			ok = ok && ScatterPoint_Participating( point, Avals );
		} else{
			ok = ok && ScatterPoint_Boundary( point, Avals );
		}
	} else{
        if( !point.IsGroundPoint() ){
    		ok = ok && ScatterPoint_Scalar( point, Avals );
        } else{
            ok = ok && ScatterPoint_Boundary( point, Avals ); // Just always do a polarized scatter at the ground for now until I can figure out what Nick is actually doing at the ground points. 
        }
	}

	return ok;
}



bool RadStore_Polarized::ScatterPoint_Participating( const SKTRAN_HR_Diffuse_Point& point, const Avals_Base& Avals )
{
    bool ok = true;

    const size_t      numincoming = point.NumIncomingRays();
    const size_t      numoutgoing = point.NumOutGoingRays();
    SKTRAN_Stokes_NC  sum;
    SKTRAN_Stokes_NC  incomingrad;

    std::unique_ptr< EtaCalculator_Base > etaCalculator; 
    Avals.CreateEtaCalculator( etaCalculator ); 
    etaCalculator->SetPoint( point );
    std::unique_ptr< Aval_ScatMatIteratorManager_Base > avalIterManager( std::move(Avals.MatrixAvalIteratorManager( point )) ); // TODO: Actually use this 

    HELIODETIC_UNITVECTOR outray; 

    // Unit vectors for the "local" reference frame
    const size_t numQPtVecs = 3;
    HELIODETIC_UNITVECTOR qPtVecs[numQPtVecs];
    point.Location().LocalUnitVectors(qPtVecs, numQPtVecs); 

    for( size_t outidx = 0; outidx < numoutgoing; outidx++ )
    {
        sum.SetTo(0.0);

        nxVector outray_unrotated;
        point.OutgoingRayLocalCoords( outidx, outray_unrotated );
        outray = (HELIODETIC_VECTOR(qPtVecs[0],outray_unrotated.X()) + HELIODETIC_VECTOR(qPtVecs[1],outray_unrotated.Y()) + HELIODETIC_VECTOR(qPtVecs[2],outray_unrotated.Z())).UnitVector(); // Get local incoming LOOK direction

        etaCalculator->UpdateOutgoingIndex( point, outidx );

        for( size_t inidx = 0; inidx < numincoming; inidx++ )
        {
            // Get incoming stokes vector and phase matrix for this scattering angle
            incomingrad = m_incomingvector[point.IncomingRadianceIdx(inidx)];

            // Compute rotation matrices from spherical geometry
            const HELIODETIC_UNITVECTOR inray(point.IncomingRayGlobalCoords(inidx));
            double l1l2 = -(inray & outray ); // prop_1=-ell_1, dot prop_2=ell_2 

            etaCalculator->CalculateEtas( point, inidx ); // Do geometry for this scatter
            etaCalculator->RefToScatt( incomingrad );
            Avals.ApplyValue( point, inidx, outidx, &incomingrad ); // Apply scattering matrices 
            etaCalculator->ScattToRef( incomingrad );

            sum += incomingrad; // Add to total sum
        }

        ScatterPoint_StoreData( point.OutgoingRadianceIdx(outidx), sum );
    }

    return ok;
}


bool RadStore_Polarized::ScatterPoint_Boundary( const SKTRAN_HR_Diffuse_Point& point, const Avals_Base& Avals )
{
	bool					ok = true;
	const size_t			numincoming = point.NumIncomingRays();
	SKTRAN_Stokes_NC		incomingrad;
	HELIODETIC_UNITVECTOR	outray; 
	size_t					fluxidx ;

	NXASSERT((point.IsGroundPoint()));							

    ok = ok && (numincoming==point.NumGroundDownwardFlux()); // so the number of downward flux should equal the number of incoming
    if(ok){
	    for( size_t inidx = 0; inidx < numincoming; inidx++ )
	    {
			    fluxidx         = point.GroundDownwardFluxIdx( inidx );		// This is where we will write the downward flux 
			    incomingrad     = m_incomingvector[point.IncomingRadianceIdx(inidx)];
			    Avals.ApplyValue( point, inidx, 0, &incomingrad );
			    //sum += incomingrad;            // Just taking scalar component so we don't have to rotate into outray's basis
			    //sum.Assign_Q ( 0.0 );
			    //sum.Assign_U ( 0.0 );

			    ScatterPoint_StoreData( fluxidx, incomingrad );
	    }
    } else{
        nxLog::Record(NXLOG_WARNING,"RadStore_Polarized::ScatterPoint_Boundary, number of incoming rays does not equal numdownward flux. Thats not good");
    }
	return ok;
}


void RadStore_Polarized::ScatterPoint_StoreData( size_t outidx, const SKTRAN_Stokes_NC& value )
{
	double norm = 1e-20<value.I() ? 1.0/value.I() : 0.0;
	m_outgoingscalar      [ outidx ] = SKTRAN_HR_WEIGHT_TYPE(value.I());
	m_outgoingLinPols_lastOrder [ outidx ] = linPolComponents( m_outgoingLinPols[outidx].Q, m_outgoingLinPols[outidx].U );
	m_outgoingLinPols     [ outidx ] = linPolComponents( norm*value.Q(), norm*value.U() );
	m_totaloutgoingscalar [ outidx ] += m_outgoingscalar[outidx]; // This is only used to calculate if the polarization ratios have converged
}


bool RadStore_Polarized::ComputeNextOrderPoint( const SKTRAN_HR_Diffuse_Point& point )
{
	bool ok = true;

	if( m_order_incoming <= m_maxPolarizedOrder_integration ){
		ok = ok && ComputeNextOrderPoint_IntegrateVectors( point );
	} else{
		if (!point.IsGroundPoint()) {
			ok = ok && ComputeNextOrderPoint_IntegrateScalars( point );
		}
		else {
			ok = ok && ComputeNextOrderPoint_IntegrateVectors( point );
		}
	}

	return ok;
}


bool RadStore_Polarized::DeclareAllScattered ( ) 
{
	bool ok = RadStore_Scalar::DeclareAllScattered();

	if (m_maxPolarizedOrder_scattering != 0 && m_maxPolarizedOrder_integration != 0)
	{
		size_t numoutgoing = GetNumOutgoing();

		std::vector<double> qDifference;
		std::vector<double> uDifference;

		// Approximate error is going to be (scalar_outgoing) / (scalar_totaloutgoing) * (abs difference in Q/U from last iteration) 

		qDifference.resize(numoutgoing);
		uDifference.resize(numoutgoing);

		for (int outidx = 0; outidx < numoutgoing; outidx++)
		{
			double msFraction = m_outgoingscalar[outidx] / m_totaloutgoingscalar[outidx];

			qDifference[outidx] = msFraction * std::abs(m_outgoingLinPols[outidx].Q - m_outgoingLinPols_lastOrder[outidx].Q);
			uDifference[outidx] = msFraction * std::abs(m_outgoingLinPols[outidx].U - m_outgoingLinPols_lastOrder[outidx].U);
		}

		double maxQ = *std::max_element(std::begin(qDifference), std::end(qDifference));
		double maxU = *std::max_element(std::begin(uDifference), std::end(uDifference));
		if (std::max(maxQ, maxU) < m_adaptivefraction)
		{
			m_maxPolarizedOrder_integration = 0;
			m_maxPolarizedOrder_scattering = 0;
		}
	}


	return ok;
}