#pragma once

class skOpticalProperties_HitranChemical;

/*---------------------------------------------------------------------------
 *                  Class hitran_geodetic_point                   2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/

class hitran_geodetic_point
{
	public:
		double m_height;
		double m_latitude;
		double m_longitude;

	public:
					hitran_geodetic_point();
					hitran_geodetic_point( const GEODETIC_INSTANT& point);
		void		FromGeodeticInstant( const GEODETIC_INSTANT& point);
		bool		operator <  (const hitran_geodetic_point& other ) const;
		bool		operator == (const hitran_geodetic_point& other ) const;
};


/*---------------------------------------------------------------------------
 *                 Class Hitran_CrossSection_Cache                  2019-11-06 */
/** This class is used to try and help speed up engines such as HR which
 *	unable to make use of the speed optimized skOpticalProperties_HitranChemical::CalculateCrossSectionsArray.
 *  This class internally calls CalculateCrossSectionsArray whenever any new location
 *	in the atmosphere is encountered.
 optimize **/
/*---------------------------------------------------------------------------*/

class Hitran_CrossSection_Cache
{
	private:
		skOpticalProperties_HitranChemical*						m_parent;
		nx1dArray<double>										m_wavenum;				// Cached Wavenumbers in ascending order.
		const double*											m_wbegin;				// Points to the beginning of m_wavenum;
		const double*											m_wend;					// Points to the end of m_wavenum;

		nx1dArray<double>										m_workerextxs;			// Worker array only used when calculating cross-sections.
		nx1dArray<double>										m_workerscatxs;			// Worker array only used when calculating cross-sections
		nx1dArray<double>*										m_current_absxs;		// Current cached cross-sections 
		nx1dArray<double>										m_blank_entry;			// A blank entry, so we dont have to use null ptrs when there are errors.

	private:
		        std::map< hitran_geodetic_point, nx1dArray<double> >				m_cached_entries;			// Cached cross sections from calls to SetLocation
		typedef std::map< hitran_geodetic_point, nx1dArray<double> >::iterator		iterator;
		typedef std::map< hitran_geodetic_point, nx1dArray<double> >::value_type	value_type;

	private:
		bool													CreateNewEntry					( const GEODETIC_INSTANT& geo_pt, iterator* iter);

	public:
																Hitran_CrossSection_Cache		( skOpticalProperties_HitranChemical* parent);
		bool													SetLocation						( const GEODETIC_INSTANT& geo_pt );
		bool													CalculateCrossSections			( double wavenumber, double *absxs, double* extxs, double* scattxs );
		bool													SetCachedWavenumbers			( const nx1dArray<double>& wavenumbers );									
		bool													HasWavenumberAlreadyInCache		( double wavenum) const;
};



