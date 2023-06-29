#include <skclimatology21.h>
#include <array>
typedef std::vector< std::pair< CLIMATOLOGY_HANDLE, nx2dArray<double> > >::iterator Iterator;


skClimatology_UserDefinedPlane::skClimatology_UserDefinedPlane()
{
	// Set some values for error checking later
	m_normal.SetCoords(0.0,0.0,0.0);
	m_reference.SetCoords(0.0,0.0,0.0);
}

skClimatology_UserDefinedPlane::~skClimatology_UserDefinedPlane()
{

}

bool skClimatology_UserDefinedPlane::UpdateCache( const GEODETIC_INSTANT& placeandtime ) 
{
	// don't have any values to cache
	return true;
}

bool skClimatology_UserDefinedPlane::GetParameter(	const CLIMATOLOGY_HANDLE& species,
													const GEODETIC_INSTANT& placeandtime,
													double* value,
													bool updatecache )
{
	Iterator it = IteratorToProfile( species );
	if( it == std::end( m_profiles ) )
	{
		nxLog::Record(NXLOG_WARNING, "skClimatology_UserDefinedPlane::GetParamater, trying to get paramater for an unsupported species" );
		return false;
	}
	if( !IsInValidState() )
	{
		// TODO: We dont need to check for valid state every call here
		nxLog::Record(NXLOG_WARNING, "skClimatology_UserDefinedPlane::GetParamater, table is not in a valid state" );
		return false;
	}
	bool dolog = m_dologinterp[it - std::begin(m_profiles)];
	double angle = ProjectedAngle( placeandtime );
	*value = InterpolateProfile( angle, placeandtime.heightm, it->second, dolog );
	return true;
}

Iterator skClimatology_UserDefinedPlane::IteratorToProfile( const CLIMATOLOGY_HANDLE& species )
{
	return std::find_if( std::begin( m_profiles ), std::end( m_profiles ),
		[&species] ( const std::pair< CLIMATOLOGY_HANDLE, nx2dArray<double> >& in ) { return in.first == species; } );
}

bool skClimatology_UserDefinedPlane::IsSupportedSpecies( const CLIMATOLOGY_HANDLE& species )
{
	Iterator it = IteratorToProfile( species );

	if( it == std::end( m_profiles ) )
	{
		return false;
	}
	else
	{
		return true;
	}
}

bool skClimatology_UserDefinedPlane::SetHeightGrid( const std::vector<double>& heights )
{
	bool ok = std::is_sorted( std::begin( heights ), std::end( heights ) );
	if( !ok )
	{
		nxLog::Record(NXLOG_WARNING, "skClimatology_UserDefinedPlane::SetHeightGrid, input height grid was not in ascending order" );
		return ok;
	}
	m_heights = heights;
	return ok;
}

bool skClimatology_UserDefinedPlane::SetAngleGrid( const std::vector<double>& angles )
{
	bool ok = std::is_sorted( std::begin( angles ), std::end( angles ) );
	if( !ok )
	{
		nxLog::Record(NXLOG_WARNING, "skClimatology_UserDefinedPlane::SetAngleGrid, input angle grid was not in ascending order" );
		return ok;
	}
	m_angles = angles;
	return ok;
}

bool skClimatology_UserDefinedPlane::SetPlane( const nxVector& normal, const nxVector& referenceinplane )
{
	bool ok = normal.IsValid() && referenceinplane.IsValid();
	if( !ok )
	{
		nxLog::Record(NXLOG_WARNING, "skClimatology_UserDefined::SetPlane, one of normal or referenceinplane was zero" );
		return ok;
	}
	ok = ok && fabs(referenceinplane & normal) < 1E-10;
	if( !ok )
	{
		nxLog::Record(NXLOG_WARNING, "skClimatology_userDefinedPlane::SetPlane, normal and referenceinplane must be orthogonal" );
		return ok;
	}
	m_normal = normal.UnitVector();
	m_reference = referenceinplane.UnitVector();

	return ok;
}


bool skClimatology_UserDefinedPlane::AddSpecies( const nx2dArray<double>& profile, const CLIMATOLOGY_HANDLE& species, bool dolog )
{
	bool ok = true;
	// AddSpecies can also be called to update a profile, so first check if the species already exists in the container
	// and erase it if it does

	for (auto it = std::begin(m_profiles); it != std::end(m_profiles); it++)
	{
		if (it->first == species)
		{
			m_dologinterp.erase(std::begin(m_dologinterp) + std::distance(std::begin(m_profiles), it));
			m_profiles.erase(it);

			// We can only have 1 duplicate so we can stop here
			break;
		}
	}

	m_profiles.emplace_back( std::pair< CLIMATOLOGY_HANDLE, nx2dArray<double> > ( species, profile ) );
	if( dolog )
	{
		nx2dArray<double>& prof = m_profiles[m_profiles.size()-1].second;
		// Need to take log of profile
		for( size_t idx = 0; idx < prof.XSize(); idx++ )
		{
			for( size_t idy = 0; idy < prof.YSize(); idy++ )
			{
				if( prof.At(idx, idy) < 0 )
				{
					nxLog::Record(NXLOG_WARNING, "Adding a negative element to a climatology using log interpolation");
					ok = false;
				}
				prof.At(idx, idy) = log(prof.At(idx, idy));
			}
		}
	}
	m_dologinterp.push_back(dolog);
	return ok;
}

bool skClimatology_UserDefinedPlane::IsInValidState()
{
	bool ok = true;
	size_t numangles, numheights;
	numangles = m_angles.size();
	numheights = m_heights.size();
	ok = ok && (numangles > 0);
	ok = ok && (numheights > 0);

	for( const auto& ele : m_profiles )
	{
		ok = ok && ele.second.YSize() == numangles;
		ok = ok && ele.second.XSize() == numheights;
	}
	ok = ok && m_normal.IsValid();
	ok = ok && m_reference.IsValid();

	return ok;
}

double skClimatology_UserDefinedPlane::ProjectedAngle( const GEODETIC_INSTANT& placeandtime )
{
	// convert the input lat/lon to geographic unitvector
	double lat = placeandtime.latitude;
	double lon = placeandtime.longitude;
	nxVector posunit( nxmath::cosd(lat) * nxmath::cosd(lon),
					  nxmath::cosd(lat) * nxmath::sind(lon),
					  nxmath::sind(lat) );

	// subtract the normal component to get the projection
	nxVector projected = posunit - (posunit & m_normal)*m_normal; // edit 2022-01-27 luf542 `posunit - (posunit & m_normal)*posunit` -> `posunit - (posunit & m_normal)*m_normal` (I don't think this actually changes the angle in the end)

	if( !projected.IsValid() )
	{
		// location was exactly on the normal vector
		return 0.0;
	}

	projected = projected.UnitVector();
	double y = (m_reference.Cross( m_normal ) ) & projected;
	double x = (m_reference & projected );
	
	return nxmath::atan2d(y,x);
}

bool skClimatology_UserDefinedPlane::LinearInterpIndexAndWeight( double val,
																 const std::vector<double>& table,
																 std::array<size_t,2>& index,
																 std::array<double,2>& weight,
																 bool zeropastboundary )
{
	bool ok = true;

	auto it = std::upper_bound( std::begin(table), std::end(table), val );	

	if( it == std::begin(table) )
	{
		// truncate to the beginning of the table
		index[0] = 0;
		index[1] = 0;
		weight[0] = zeropastboundary ? 0 : 1;
		weight[1] = 0;
	}
	else if( it == std::end(table) )
	{
		// truncate to the end of the table
		index[0] = (it-std::begin(table))-1;
		index[1] = 0;
		weight[0] = zeropastboundary ? 0 : 1;
		weight[1] = 0;
	}
	else
	{
		// inside the table
		size_t highindex  = (it - std::begin(table));
		size_t lowindex = highindex - 1;
		index[0] = lowindex;
		index[1] = highindex;

		weight[0] = (table[highindex] - val) / (table[highindex] - table[lowindex]);
		weight[1] = 1.0 - weight[0];
	}

	return ok;
}

double skClimatology_UserDefinedPlane::InterpolateProfile( double angle, double height, const nx2dArray<double>& profile, bool dolog )
{
	std::array<size_t,2>	angleindex, heightindex;
	std::array<double,2>	angleweight, heightweight;

	LinearInterpIndexAndWeight( angle, m_angles, angleindex, angleweight, false );
	LinearInterpIndexAndWeight( height, m_heights, heightindex, heightweight, true );

	double result = 0.0;
	double hresult;
	for( size_t angleidx = 0; angleidx < angleindex.size(); angleidx++ )
	{
		hresult = 0.0;
		for( size_t heightidx = 0; heightidx < heightindex.size(); heightidx++ )
		{
			hresult += profile.At( heightindex[heightidx], angleindex[angleidx] ) * heightweight[heightidx] * angleweight[angleidx];
		}
		if( dolog )
		{
			result += exp(hresult);
		}
		else
		{
			result += hresult;
		}
	}
	return result;
}