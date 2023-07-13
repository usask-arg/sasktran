
/*-----------------------------------------------------------------------------
 *					skClimatology_LabowOzoneVMR		2005-7-28							*/
/** \ingroup o3skClimatology
 *	The Labow volume mixing ratio climatology of ozone. This is a climatology of the volume mixing ratio
 *	from 0 km to 60 km in 1 km increments for 18, 10 degree latitudes steps from -85 to +85 for each month of the year.
 *	The climatology seems to be unpublished as all references refer to an unpublished piece of work
 *	by McPeters and Labow in 2002 or 2003. The climatology is hard-coded into the source file so the
 *	caching is fast. We have extended the climatology above 60 km using the differences in scale height
 *	between the neutrals (7 km) and ozone (4.5 km) to calculate a new scale height to extrapolate the signal at 60 km.
 *	We have also extended the volume mixing ratio so it can be converted to an ozone density.  This requires a neutral
 *	density climatology. By default the Labo climatology uses the skClimatology_MSIS90 climatology for this purpose. It is possible to change
 *	the neutral density model used by the climatology.
 *
 *	\par Supported Species
 *	This model supports the following species:-
 *	\n
 *	-# CLIMATOLOGY_O3_VMR
 *	-# CLIMATOLOGY_O3_CM3
 *	\n
 *
 *	\par Range of Valid model output
 *	The model should be valid over the following ranges:
 *	\n
 *	-# All times and dates
 *	-# All latitudes and longitudes
 *	-# The extended version ranges from 0 km to 100 km while the
 *	   intrinisc database ranges from 0 to 60 km.
 *	
**/
/*---------------------------------------------------------------------------*/

class skClimatology_LabowOzoneVMR : public skClimatology
{
	private:
		static  double			g_labowvmrppm[12][18][61];
		skClimatology_MSIS90	m_defaultneutral;				//!< The default model used to infere the neutral density.
		nxTimeStamp				m_mjd0;							//!< Timestamp at the middle of the first month (start of 16th day)
		nxTimeStamp				m_mjd1;							//!< Timestamp at the middle of the second day, (start of 16th day)
		nx2dArray<double>		m_profile0;						//!< double( 61 heights, 18 latitudes) for first month
		nx2dArray<double>		m_profile1;						//!< double( 61 heights, 18 latitudes) for second month
		skClimatology*			m_neutralmodel;					//!< A climatology that gives the neutral air density

	private:
		bool					SetBoundingMonths	( double mjd );

	public:
								skClimatology_LabowOzoneVMR();
		virtual				   ~skClimatology_LabowOzoneVMR();
		bool					SetNeutralDensityClimatology( skClimatology* neutrals) { m_neutralmodel = neutrals; return( nxTRUE ); }
//		bool					DeepCopy			( const skClimatology_LabowOzoneVMR& other );
		virtual	bool			UpdateCache			( const GEODETIC_INSTANT& placeandtime) override;
		virtual	bool			GetParameter		( const CLIMATOLOGY_HANDLE& species,  const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache) override;
		virtual bool			IsSupportedSpecies  ( const CLIMATOLOGY_HANDLE& species ) override;
	//	virtual bool			CreateClone			 (skClimatology** clone)	const override;
};
