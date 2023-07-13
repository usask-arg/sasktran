
/*-----------------------------------------------------------------------------
 *					class Climatology_Zero				2008-2-11*/
/** \ingroup skClimmisc
 *	This class implements a climatology that is exactly zero. Although this
 *	seems pointless it is useful for legacy code that requires climatologies
 *	for specific species even if those species are not required in the
 *	specific problem being solved.
 **/
/*---------------------------------------------------------------------------*/
class skClimatology_Zero : public skClimatology
{

	public:
								skClimatology_Zero	( );
		virtual				   ~skClimatology_Zero	( );

	public:
		virtual	bool			UpdateCache			( const GEODETIC_INSTANT& placeandtime ) override;
		virtual	bool			GetParameter		( const CLIMATOLOGY_HANDLE& species, const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache) override;
		virtual bool			IsSupportedSpecies  ( const CLIMATOLOGY_HANDLE& species ) override;
//		virtual bool			CreateClone			( skClimatology** clone) const override;
};

class skClimatology_Undefined : public skClimatology
{

	public:
								skClimatology_Undefined	( );
		virtual				   ~skClimatology_Undefined	( );

	public:
		virtual	bool			UpdateCache			( const GEODETIC_INSTANT& placeandtime ) override;
		virtual	bool			GetParameter		( const CLIMATOLOGY_HANDLE& species, const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache) override;
		virtual bool			IsSupportedSpecies  ( const CLIMATOLOGY_HANDLE& species ) override;
//		virtual bool			CreateClone			( skClimatology** clone) const override;
};

