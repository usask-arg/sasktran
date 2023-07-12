#include "include/sktran_hr_internals.h"
#include <numeric>
#include <functional>
#include <boost/multi_array.hpp>
using boost::multi_array;
using boost::extents;


RadStore_Scalar::RadStore_Scalar( )
{
	m_order_incoming = 0;
	m_order_outgoing = 0;
}


bool RadStore_Scalar::AllocateStorage( size_t incomingcounter, size_t outgoingcounter )
{
	bool ok = true;
	m_incomingscalar      .resize( incomingcounter );
	m_outgoingscalar      .resize( outgoingcounter );
	m_totaloutgoingscalar .resize( outgoingcounter );
	m_totalincomingscalar .resize( incomingcounter );
	ok = ok && RadStore_Base::AllocateStorage( incomingcounter, outgoingcounter );
	return ok;
}


bool RadStore_Scalar::CleanDiffuseIndexes( )
{
	bool ok = true;
	
	m_incomingscalar      .assign( m_incomingscalar.size(),      0.0 );
	m_outgoingscalar      .assign( m_outgoingscalar.size(),      0.0 );
	m_totaloutgoingscalar .assign( m_totaloutgoingscalar.size(), 0.0 );
	m_totalincomingscalar .assign( m_totalincomingscalar.size(), 0.0 );

	ok = ok && RadStore_Base::CleanDiffuseIndexes( );
	return true;
}


bool  RadStore_Scalar::DeclareFirstOrderInitialized ( ) 
{
	// Could put some error checking in here. 
	m_order_incoming = 1;
	m_order_outgoing = 0;
	return true;
}


bool  RadStore_Scalar::DeclareAllScattered ( ) 
{
	bool ok = true;
	ok = ok && ((m_order_outgoing+1)==m_order_incoming);
	if(ok) m_order_outgoing = m_order_incoming;
	return ok;
}


bool  RadStore_Scalar::DeclareAllIntegrated ( ) 
{
	bool ok = true;
	ok = ok && (m_order_outgoing==m_order_incoming);
	++m_order_incoming;
	return true;
}


size_t RadStore_Scalar::GetDiffuseSourceOrder ( ) const
{
	return m_order_outgoing;
}


bool RadStore_Scalar::StoreFirstOrderIncoming ( const SKTRAN_SourceTermIntegrator_Base* srcintegrator, const std::vector<SKTRAN_Source_Term*>* sources, const SKTRAN_RayOptical_Base* ray, size_t index )
{
	bool ok = true;

	double radiance;
	ok = ok && srcintegrator->IntegrateSourceTerm( ray, radiance, *sources );
	m_incomingscalar [ index ] = SKTRAN_HR_DBL_TO_WEIGHT(radiance);
	m_totalincomingscalar[index] = SKTRAN_HR_DBL_TO_WEIGHT(radiance);
	NXASSERT( (NXFINITE(radiance)) );
	return ok;
}


bool RadStore_Scalar::ComputeNextOrderPoint( const SKTRAN_HR_Diffuse_Point& point )
{
	bool ok = true;

	SKTRAN_HR_WEIGHT_TYPE incoming;
	for( size_t inidx = 0; inidx < point.NumIncomingRays(); inidx++ )
	{
		const size_t incomingidx = point.IncomingRadianceIdx( inidx );

		if( m_incomingscalar[incomingidx] < 1E-10 )
		{
			// if our incoming radiance is 0 in this direction there is no reason to keep calculating
			// the next order
			m_incomingscalar[incomingidx] = 0.0;
			continue;
		}
		incoming = 0.0;
		const SKTRAN_HR_Diffuse_Index_Array& indexarray = point.DiffuseIndices()[inidx];
		const size_t numdiffuseindex = indexarray.diffindex.size();
		for( size_t idx = 0; idx < numdiffuseindex; idx++ )
		{
			SKTRAN_HR_WEIGHT_TYPE value =  indexarray.diffindex[idx].weight * m_outgoingscalar[indexarray.diffindex[idx].index];

			NXASSERT( (NXFINITE(value) ));

			incoming += value;
			
		}
		m_incomingscalar[incomingidx] = incoming;
		m_totalincomingscalar[incomingidx] += incoming;
		NXASSERT( (NXFINITE(incoming) ));
	}

	return ok;
}


bool RadStore_Scalar::DiffuseSource ( const SKTRAN_SourceTermQueryObject_Base& qNoTransform, const SKTRAN_HR_Diffuse_Point& point, const nxVector& lookunit, double& source ) const
{
	bool   ok = true;
	size_t unitsphereindex  [3];
	double unitsphereweight [3];

	ok = ok && point.TriangulateOnOutgoing( lookunit, unitsphereindex, unitsphereweight, 3);

	source = 0.0;
	for( size_t triidx = 0; triidx < 3; triidx++ )
	{
		source += m_totaloutgoingscalar[ point.OutgoingRadianceIdx(unitsphereindex[triidx]) ] * unitsphereweight[triidx];
	}

	return true;
}


bool RadStore_Scalar::DiffuseSource ( const SKTRAN_SourceTermQueryObject_Base& qNoTransform, const SKTRAN_HR_Diffuse_Point& point, const nxVector& look, SKTRAN_Stokes_NC& source ) const
{
	bool ok = true;
	
	double scalarRadiance;
	source.SetTo( 0.0 );
	ok = ok && DiffuseSource( qNoTransform, point, look, scalarRadiance );
	source.Assign_I( scalarRadiance );

	return ok;
}


bool RadStore_Scalar::GroundSource(const SKTRAN_SourceTermQueryObject_Base& qobj, const SKTRAN_HR_Diffuse_Point& point, const nxVector&, double&             source) const
{
	bool							ok = true;
	double							brdf;
	double							factor;
	double							mu_in;
	double							mu_out;
	double							cosdphi;
	size_t							profidx = 0;
	HELIODETIC_UNITVECTOR			localcoords[3];
	HELIODETIC_UNITVECTOR			outgoingray_helio;
	HELIODETIC_UNITVECTOR			outgoingray_local;
	nxVector						outgoinghoriz;
	bool							ok1;

	size_t							tempindex;
	
	source = 0.0;
	const HELIODETIC_POINT&  groundpoint = qobj.GetPoint();					
	outgoingray_helio = qobj.GetLookAway();

	const HELIODETIC_UNITVECTOR& localzenith_helio = groundpoint.LocalZenith();						// Get the local zenith at the point where the ray hits the ground
	mu_out = -1*(localzenith_helio & outgoingray_helio);											// Get the cosine of zenith angle of outbound ray
	NXASSERT((mu_out >= 0.0));																		// and make sure we have it right

	groundpoint.LocalUnitVectors(localcoords, 3);													// Get North, West and Up at the current groundpoint in heliodetic corodinates
	outgoingray_local = groundpoint.TransformToLocalZenithCoords(outgoingray_helio, localcoords);	// Transform the outgoing ray to the local coordinate system (in the local system, up is always (0,0,1)
	outgoinghoriz.SetCoords(outgoingray_local.X(), outgoingray_local.Y(), 0.0);						// Get just the horizontal component as we will need this later for the azimuth determination
	outgoinghoriz.UnitVector();

	

	for (size_t inidx = 0; inidx < point.NumGroundDownwardFlux(); inidx++)							// For each incoming downward flux
	{																								// Loop over incoming directions
		nxVector	incoming(point.IncomingUnitRayLocalCoords(inidx));								// Get the full incoming ray unit vector (away from ground point) in local coords	
		nxVector	incominghoriz(incoming.X(), incoming.Y(), 0);									// Get the horizontal component of the look vector 
		incominghoriz = incominghoriz.UnitVector();													// and get the horizontal component as a unit vector

		mu_in = incoming.Z();																		// cosine of incoming ray zenith angle
		NXASSERT((mu_in >= 0.0));																	// make sure it is acceptable
		cosdphi = outgoinghoriz & incominghoriz;													// get thecosine of the azimuth angle between rays
		ok1 = m_opticaltable->GetBRDF(groundpoint, mu_in, mu_out, cosdphi, &brdf);						// Do the BRDF calculation
		factor = (ok1) ? brdf*point.InCubatureWeight(inidx)*mu_in : 0;								// How much source function does this downward flux add to the ray we are considering 
		tempindex = point.IncomingRadianceIdx(inidx);

		source += m_totalincomingscalar[tempindex] * factor;
	}
	return ok;
}


bool RadStore_Scalar::GroundSource(const SKTRAN_SourceTermQueryObject_Base& qNoTransform, const SKTRAN_HR_Diffuse_Point& point, const nxVector& look, SKTRAN_Stokes_NC& source) const
{
    bool ok = true;
    double scalarRadiance;
    ok = ok && GroundSource( qNoTransform, point, look, scalarRadiance );
    source.SetTo(0.0);
    source.Assign_I( scalarRadiance );
    return ok;
}


bool RadStore_Scalar::DumpIncomingRadiances( const SKTRAN_HR_Diffuse_Table_CPU* table, size_t order, double wlen )
{
	double              x, y, z;
	double              radiance;
	char				dsetname[100];

	std::vector<size_t> indices = table->DiagnosticProfIdxs();

	size_t nalts = table->DiagnosticProfAlts();
	size_t npts = table->DiffusePointAt( indices[0] ).NumIncomingRays();
	size_t ndata = 4; // [x, y, z, rad]

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
				radiance = DumpIncomingRadiances_GetIncomingRadiance( table->DiffusePointAt( pointidx ).IncomingRadianceIdx( rayidx ) );

				maindata[pointidx - indices[idx]][rayidx][0] = x;
				maindata[pointidx - indices[idx]][rayidx][1] = y;
				maindata[pointidx - indices[idx]][rayidx][2] = z;
				maindata[pointidx - indices[idx]][rayidx][3] = radiance;
			}
		}
		// write dataset for this profile
		sprintf( dsetname, "in_wlen_%0.2f_ord_%d_prof_%d", wlen, int( order ), int( indices[idx] / nalts ) );		
		H5Did = H5Dcreate2( H5Fid, dsetname, H5T_NATIVE_DOUBLE, H5Sid, H5P_DEFAULT, propList, H5P_DEFAULT );
		H5Dwrite( H5Did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, maindata.data() );
		H5Dclose( H5Did );
	}
	H5Pclose( propList );
	H5Sclose( H5Sid );

	// and dump data for ground points
	size_t ngp = table->NumGroundPoints();
	size_t ngr = table->GroundPointAt( 0 ).NumIncomingRays();

	multi_array<double, 2> grounddata( extents[ngp][ngr] );
	hsize_t dims2[2] = { ngp, ngr };

	H5Sid = H5Screate_simple( 2, dims2, NULL );
	propList = H5Pcreate( H5P_DATASET_CREATE );
	H5Pset_layout( propList, H5D_CHUNKED );
	H5Pset_chunk( propList, 2, dims2 );

	for( size_t pointidx = 0; pointidx < table->NumGroundPoints(); pointidx++ )
	{
		for( size_t rayidx = 0; rayidx < table->GroundPointAt( pointidx ).NumIncomingRays(); rayidx++ )
		{
			grounddata[pointidx][rayidx] = m_incomingscalar[table->GroundPointAt( pointidx ).IncomingRadianceIdx( rayidx )];
		}
	}
	sprintf( dsetname, "in_wlen_%0.2f_ord_%d_ground", wlen, int( order ) );
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


bool RadStore_Scalar::DumpOutgoingRadiances( const SKTRAN_HR_Diffuse_Table_CPU* table, size_t order, double wlen )
{
	nxVector            outgoingRay;
	double              radiance;
	char				dsetname[100];

	std::vector<size_t> indices = table->DiagnosticProfIdxs();

	size_t nalts = table->DiagnosticProfAlts();
	size_t npts = table->DiffusePointAt( indices[0] ).NumOutGoingRays();
	size_t ndata = 4;

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
			for( size_t rayidx = 0; rayidx < table->DiffusePointAt( pointidx ).NumOutGoingRays(); rayidx++ )
			{
				table->DiffusePointAt( pointidx ).OutgoingRayLocalCoords( rayidx, outgoingRay );
				radiance = m_outgoingscalar[table->DiffusePointAt( pointidx ).OutgoingRadianceIdx( rayidx )];

				maindata[pointidx - indices[idx]][rayidx][0] = outgoingRay.X();
				maindata[pointidx - indices[idx]][rayidx][1] = outgoingRay.Y();
				maindata[pointidx - indices[idx]][rayidx][2] = outgoingRay.Z();
				maindata[pointidx - indices[idx]][rayidx][3] = radiance;
			}
		}
		// write dataset for this profile
		sprintf( dsetname, "out_wlen_%0.2f_ord_%d_prof_%d", wlen, int( order ), int( indices[idx] / nalts ) );
		H5Did = H5Dcreate2( H5Fid, dsetname, H5T_NATIVE_DOUBLE, H5Sid, H5P_DEFAULT, propList, H5P_DEFAULT );
		H5Dwrite( H5Did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, maindata.data() );
		H5Dclose( H5Did );
	}
	H5Pclose( propList );
	H5Sclose( H5Sid );

	// and dump data for ground points
	size_t ngp = table->NumGroundPoints();
	hsize_t dims1[1] = { ngp };

	H5Sid = H5Screate_simple( 1, dims1, NULL );
	propList = H5Pcreate( H5P_DATASET_CREATE );
	H5Pset_layout( propList, H5D_CHUNKED );
	H5Pset_chunk( propList, 1, dims1 );

	multi_array<double, 1> grounddata( extents[ngp] );

	for( size_t pointidx = 0; pointidx < table->NumGroundPoints(); pointidx++ )
	{
		grounddata[pointidx] = m_outgoingscalar[table->GroundPointAt( pointidx ).OutgoingRadianceIdx( 0 )];
	}
	sprintf( dsetname, "out_wlen_%0.2f_ord_%d_ground", wlen, int( order ) );
	H5Did = H5Dcreate2( H5Fid, dsetname, H5T_NATIVE_DOUBLE, H5Sid, H5P_DEFAULT, propList, H5P_DEFAULT );
	H5Dwrite( H5Did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, grounddata.data() );
	H5Dclose( H5Did );

	H5Pclose( propList );
	H5Sclose( H5Sid );
	H5Fclose( H5Fid );
	maindata.resize( extents[0][0][0] );
	grounddata.resize( extents[0] );
	return true;
}



/*-----------------------------------------------------------------------------
 *					RadStore_Scalar::ScatterPoint		 2016- 12- 1*/
/** **/
/*---------------------------------------------------------------------------*/

bool RadStore_Scalar::ScatterPoint( const SKTRAN_HR_Diffuse_Point& point, const Avals_Base& Avals )
{
	return ScatterPoint_Scalar( point, Avals );
}


/*-----------------------------------------------------------------------------
 *					RadStore_Scalar::ScatterPoint_Scalar		 2016- 12- 1*/
/** **/
/*---------------------------------------------------------------------------*/

bool RadStore_Scalar::ScatterPoint_Scalar( const SKTRAN_HR_Diffuse_Point& point, const Avals_Base& Avals )
{
	bool ok = true;
	const size_t          numincoming = point.NumIncomingRays();

	SKTRAN_HR_WEIGHT_TYPE incoming;

    std::unique_ptr< Aval_ScalarIteratorManager_Base > scalarAvalManager ( Avals.ScalarAvalIteratorManager( point ) );
    
	if (point.IsGroundPoint())																// IF we have a ground point
	{																						// then the outgoing signal saved is really the incoming downward flux of each incoming ray =  ( I.cos(theta).domega))
		NXASSERT( (numincoming==point.NumGroundDownwardFlux()) );							// so the number of downward flux should equal the number of incoming
		ok = (numincoming==point.NumGroundDownwardFlux());
		if (ok)
		{
			for( size_t inidx = 0; inidx < numincoming; inidx++ )
			{
				size_t                fluxidx       = point.GroundDownwardFluxIdx( inidx );		// This is where we will write the downward flux 
				size_t				  incomingindex = point.IncomingRadianceIdx  ( inidx );		// This is the incoming radiance
				SKTRAN_HR_WEIGHT_TYPE value;

				incoming = m_incomingscalar[incomingindex];										// Get the incoming radiance
                value = scalarAvalManager->ApplyAndAdvance( incoming );                         // Get the downward flux, multiply by pre-stored values in AVals (See Avals_ScalarStore::CalcScatteringMatrixPoint_Boundary)
				m_outgoingscalar[ fluxidx]        = value;										// Store the downward flux
				m_totaloutgoingscalar[ fluxidx ] += value;										// Store the total downward flux (not sure what this does NDL 2016-12-01)
				NXASSERT( (NXFINITE(value)));
			}
		}
		else
		{
			nxLog::Record(NXLOG_WARNING,"RadStore_Scalar::ScatterPoint_Scalar, number of incoming rays does not equal numdownward flux. Thats not good");
		}
	}
	else
	{
		const size_t  numoutgoing = point.NumOutGoingRays();

		for( size_t outidx = 0; outidx < numoutgoing; outidx++ )
		{
			m_outgoingscalar[point.OutgoingRadianceIdx( outidx )] = 0.0;
		}
		for( size_t inidx = 0; inidx < numincoming; inidx++ )
		{
			if( m_incomingscalar[point.IncomingRadianceIdx( inidx )] < 1E-10 )
			{
                scalarAvalManager->Advance( numoutgoing );
				continue;
			}
			incoming = m_incomingscalar[point.IncomingRadianceIdx( inidx )];
            
			for( size_t outidx = 0; outidx < numoutgoing; outidx++ )
			{
				m_outgoingscalar[point.OutgoingRadianceIdx(outidx)] += scalarAvalManager->ApplyAndAdvance( incoming );
			}
		}
		for( size_t outidx = 0; outidx < numoutgoing; outidx++ )
		{
            NXASSERT( NXFINITE(m_outgoingscalar[point.OutgoingRadianceIdx(outidx)]) );
			m_totaloutgoingscalar[ point.OutgoingRadianceIdx( outidx ) ] += m_outgoingscalar[ point.OutgoingRadianceIdx(outidx) ];
		}
	}

	return ok;
}


