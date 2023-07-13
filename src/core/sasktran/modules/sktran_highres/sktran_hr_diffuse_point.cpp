#include "include/sktran_hr_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Point::SKTRAN_HR_Diffuse_Point		2013-06-20*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_HR_Diffuse_Point::SKTRAN_HR_Diffuse_Point()
{
	m_incomingsphere    = nullptr;
	m_outgoingsphereobj = nullptr;
	m_isgroundpoint     = false;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Point::~SKTRAN_HR_Diffuse_Point		2013-06-20*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_HR_Diffuse_Point::~SKTRAN_HR_Diffuse_Point()
{
	CleanDiffuseIndex();
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Point::ReleaseResources		2013-06-20*/
/** Releases the internal resources, should not have to be called manually 
 *
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Point::ReleaseResources()
{
	if( nullptr != m_incomingsphere    ) m_incomingsphere->Release();
	if( nullptr != m_outgoingsphereobj ) m_outgoingsphereobj->Release();

	m_incomingsphere    = nullptr;
	m_outgoingsphereobj = nullptr;
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Point::DumpFirstOrderRadiance		2013-06-20*/
/** Dumps the first order radiances to a text file with a given file name,
 *  used for debugging purposes
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Point::DumpFirstOrderRadiance( std::string filename )
{
	std::ofstream fileout;

	double x,y,z;
	double radiance;

	fileout.open( filename.c_str() );

	for( size_t rayidx = 0; rayidx < m_incomingsphere->NumUnitVectors(); rayidx++ )
	{
		x = m_rotated_incomingunitvectors[rayidx].X();
		y = m_rotated_incomingunitvectors[rayidx].Y();
		z = m_rotated_incomingunitvectors[rayidx].Z();
		radiance = m_firstorderradiances[rayidx];

		fileout << x << " " << y << " " << z << " " << radiance << std::endl;
	}

	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Point::CreateRotated_GlobalUnitVectors		2013-06-20*/
/** The created unit sphere grid assumes that the south pole of the sphere
 *  looks directly down to the ground,however this is only true if the sphere
 *  is placed on the north pole, therefore we rotate the unit vectors 
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Point::CreateRotated_GlobalUnitVectors()
{
	bool										ok = true;
	double										x,y,z;
	HELIODETIC_UNITVECTOR						xprime;
	HELIODETIC_UNITVECTOR						yprime;
	HELIODETIC_UNITVECTOR						zprime;
	nxVector									v;
	HELIODETIC_UNITVECTOR						unitvectors[3];
	size_t										rayidx;

	m_location.LocalUnitVectors( &unitvectors[0], 3);		// Get the profile's local x',y',z' unit vectors from the parent profile
	xprime = unitvectors[0];
	yprime = unitvectors[1];
	zprime = unitvectors[2];

	m_rotated_incomingunitvectors.resize( m_incomingsphere->NumUnitVectors() );
	m_firstorderradiances.resize( m_rotated_incomingunitvectors.size() );
	if (ok)
	{
		for( rayidx = 0; rayidx < m_incomingsphere->NumUnitVectors(); rayidx++ )
		{
			v = m_incomingsphere->UnitVectorAt( rayidx );
			x = v.X();		// sinzen*cosazi;										// NOte that azimuth is zero degrees at the due south position
			y = v.Y();      // sinzen*sinazi;										// and is 90 degrees or pi/2 at the due East direction (in the system with sun at pole)
			z = v.Z();		// coszen;
			m_rotated_incomingunitvectors[rayidx].SetCoords( x*xprime.X() + y*yprime.X() + z*zprime.X(),				// look = x*xprime + y*yprime + z*zprime
							x*xprime.Y() + y*yprime.Y() + z*zprime.Y(),
							x*xprime.Z() + y*yprime.Z() + z*zprime.Z() );			// Get the look vector in heliodetic coordinates


		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Point::OutgoingRayLocalCoords		2013-06-20*/
/** Here there is an implicit assumption that the outgoing grid is uniform, 
 *  therefore it doesn't matter how we rotate it.  Thus we treat the nxVector
 *  inside the unitsphere container as heliodetic unitvectors for all interpolation
 *  purposes.
 **/
/*---------------------------------------------------------------------------*/

void SKTRAN_HR_Diffuse_Point::OutgoingRayLocalCoords( size_t idx, nxVector& outray ) const
{
	if(nullptr!=m_outgoingsphereobj){
        m_outgoingsphereobj->OutgoingRayLocalCoords( idx, outray );
    }
    else {
        //outray.SetCoords(-m_location.LocalZenith().X(), -m_location.LocalZenith().Y(), -m_location.LocalZenith().Z());
        outray.SetCoords( 0.0, 0.0, -1.0 ); // Reference direction is straight down in local coordinates
    }
}

size_t SKTRAN_HR_Diffuse_Point::NumOutGoingRays() const
{
	size_t result;
	if( m_isgroundpoint )
	{
		throw("Ground points do not have outgoing rays");
		//result =  1; // ground point, doesnt have outgoing sphere
		//result = NumIncomingRays();	// The ground point caches the  downward flux contribution (cos(theta).domega ) of each incoming ray
	} 
	else
	{
		result =  m_outgoingsphereobj->NumOutgoingRays();
	}
	return result;
}


size_t SKTRAN_HR_Diffuse_Point::NumUniqueScatterIncoming ( ) const
{
	return NumIncomingRays();
}

size_t SKTRAN_HR_Diffuse_Point::NumUniqueScatterOutgoing ( ) const
{
	return NumOutGoingRays();
}

size_t SKTRAN_HR_Diffuse_Point::UniqueScatterIncoming ( size_t inidx ) const
{
	return inidx;
}

size_t SKTRAN_HR_Diffuse_Point::UniqueScatterOutgoing ( size_t outidx ) const
{
	return outidx;
}

bool SKTRAN_HR_Diffuse_Point::SetPointIndices( size_t pointidx, size_t incomingcounter, size_t outgoingcounter, size_t scattvalcounter, size_t numdiffuseindices)
{
	m_pointidx    = pointidx;
	m_inidx		= incomingcounter;
	m_outidx		= outgoingcounter;
	m_scatvalidx	= scattvalcounter;
	m_diffuseindexes.resize( numdiffuseindices );
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Point::ConfigureSpheres		 2016- 12- 1*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Point::ConfigureSpheres( const SKTRAN_UnitSphere_V2* incomingsphere, const SKTRAN_HR_OutgoingSphereObject_Base*  outgoingsphereobj, const HELIODETIC_POINT& location, bool isground)
{
	bool	 ok = true;

	if (incomingsphere   != nullptr ) incomingsphere->AddRef();
	if (m_incomingsphere != nullptr ) m_incomingsphere->Release();

	m_incomingsphere= incomingsphere;
	if( !isground )
	{
		if (outgoingsphereobj   != nullptr) outgoingsphereobj->AddRef();
		if (m_outgoingsphereobj != nullptr) m_outgoingsphereobj->Release();
		m_outgoingsphereobj = outgoingsphereobj;
	}
	else
	{
		if (m_outgoingsphereobj != nullptr) m_outgoingsphereobj->Release();
		m_outgoingsphereobj =nullptr;
	}
	m_location = location;
	m_isgroundpoint = isground;
	CreateRotated_GlobalUnitVectors();
	return ok;
}



bool SKTRAN_HR_Diffuse_Point::CleanDiffuseIndex()
{
	bool ok = true;
	if( m_diffuseindexes.size() == 0 )
		return ok;
	for( size_t inidx = 0; inidx < m_incomingsphere->NumUnitVectors(); inidx++ )
	{
		m_diffuseindexes[inidx].diffindex.clear();
		#if defined(_MSC_VER)
			m_diffuseindexes[inidx].diffindex.shrink_to_fit();
		#else
			m_diffuseindexes[inidx].diffindex.swap( m_diffuseindexes[inidx].diffindex );
		#endif
	}
	// m_diffuseindexes.clear();


	return ok;
}

bool CompareDiffuseIndex( const SKTRAN_HR_Diffuse_Index& index1, const SKTRAN_HR_Diffuse_Index& index2 )
{
	return index1.index < index2.index;
}

bool EqualDiffuseIndex( const SKTRAN_HR_Diffuse_Index& index1, const SKTRAN_HR_Diffuse_Index& index2 )
{
	return index1.index == index2.index;
}

bool SKTRAN_HR_Diffuse_Point::CorrectForIndexDuplicates( size_t rayidx )
{
	bool ok = true;
	
	std::vector<SKTRAN_HR_Diffuse_Index>& diffind = m_diffuseindexes[rayidx].diffindex;
	std::sort( diffind.begin(), diffind.end(), CompareDiffuseIndex );
	ok = ok && CorrectForIndexDuplicates_Sorted( rayidx );

	return ok;
}

bool SKTRAN_HR_Diffuse_Point::CorrectForIndexDuplicates_Sorted( size_t sortedray_rayidx )
{
	bool ok = true;
	size_t		idxcount = 0;
	size_t		idx		 = 0;
	std::vector<SKTRAN_HR_Diffuse_Index>& diffind = m_diffuseindexes[sortedray_rayidx].diffindex;
	std::vector<SKTRAN_HR_Diffuse_Index>::iterator it;
	if( diffind.size() == 0 )
		return true;
	
	while( idx < diffind.size()-1 )
	{
		idxcount = 1;
		while( diffind[idx].index == diffind[idx+idxcount].index )
		{
			diffind[idx].weight += diffind[idx+idxcount].weight;
			NXASSERT(( NXFINITE(diffind[idx].weight)));
			++idxcount;
			if( idx + idxcount >= diffind.size() )
			{
				break;
			}
		}
		idx += idxcount;
	}

	it = std::unique( diffind.begin(), diffind.end(), EqualDiffuseIndex );

	diffind.resize( std::distance( diffind.begin(), it ) );
	diffind.shrink_to_fit();


	return ok;
}


bool SKTRAN_HR_Diffuse_Point::TriangulateOnOutgoing( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t maxvertices ) const
{
	return m_outgoingsphereobj->TriangulateOnOutgoing( unit, unit_indexptr, unit_weightptr, maxvertices );
}

double SKTRAN_HR_Diffuse_Point::OutgoingCubatureWeight( size_t outidx ) const
{
	return m_outgoingsphereobj->OutgoingCubatureWeight( outidx );
}
