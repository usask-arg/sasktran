#pragma once


/*-----------------------------------------------------------------------------
 *					UserDefined_LatLon_Table		2014-4-21*/
/** A 3D user defined class for specifying just one species. 
**/
/*---------------------------------------------------------------------------*/

class UserDefined_LatLon_Table
{
	private:
		nx3dArray<double>						m_profile;
		std::vector<double>						m_alts;
		std::vector<double>						m_lats;
		std::vector<double>						m_lons;

	private:
		bool						LinearInterpWeights							( double val, const std::vector<double>& profile,  double* weights,  size_t* indices,  size_t* numindex ) const;
		void						ReleaseResources							();

	public:
									UserDefined_LatLon_Table		();
								   ~UserDefined_LatLon_Table		();
		bool						LoadProfileFromData							( const std::vector<double>& alts,
																			      const std::vector<double>& lats,
																			      const std::vector<double>& lons,
																				  const nx3dArray<double>& profile 
																			    );
		bool						InterpTable									( double alt,  double lat, double lon, double* value, double badvalue ) const;
		bool						IsDefined									() const { return m_profile.size() > 0;}
};

/*-----------------------------------------------------------------------------
 *					class UserTableSplineEntry		2014-1-30*/
/** **/
/*---------------------------------------------------------------------------*/

class UserDefined3D_LatLonHeightEntry
{
	private:
		UserDefined_LatLon_Table				m_table;
		double												m_badvalue;

	public:
															UserDefined3D_LatLonHeightEntry	()					{ m_badvalue = std::numeric_limits<double>::quiet_NaN();}
														   ~UserDefined3D_LatLonHeightEntry	()					{};
		const UserDefined_LatLon_Table&		Table							() const			{ return m_table;}
		UserDefined_LatLon_Table*				TableVar						()					{ return &m_table;}
		bool												SetBadValue						( double badval)	{ m_badvalue = badval; return true;}
		double												BadValue						() const			{ return m_badvalue;}
};


/*-----------------------------------------------------------------------------
 *					skClimatology_UserDefined3D_LatLonHeight		2014-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

class skClimatology_UserDefined3D_LatLonHeight : public skClimatology
{
	private:
		        std::map< CLIMATOLOGY_HANDLE, UserDefined3D_LatLonHeightEntry>					m_species;
		typedef std::map< CLIMATOLOGY_HANDLE, UserDefined3D_LatLonHeightEntry>::value_type		value_type;
		typedef std::map< CLIMATOLOGY_HANDLE, UserDefined3D_LatLonHeightEntry>::iterator		iterator;


	public:
								skClimatology_UserDefined3D_LatLonHeight ( );
		virtual				   ~skClimatology_UserDefined3D_LatLonHeight ( );
		bool					LoadProfile								( const CLIMATOLOGY_HANDLE&		species,
																		  const std::vector<double>&	alts,
																		  const std::vector<double>&	lons,
																		  const std::vector<double>&	lats,
																		  const nx3dArray<double>&		profile,
																		  double						badval);
	public:
		virtual	bool			UpdateCache								( const GEODETIC_INSTANT&		placeandtime );
		virtual	bool			GetParameter							( const CLIMATOLOGY_HANDLE&		species,  const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache);
		virtual bool			IsSupportedSpecies						( const CLIMATOLOGY_HANDLE&		species );
};


