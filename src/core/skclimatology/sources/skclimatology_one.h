/*-----------------------------------------------------------------------------
 *					class Climatology_One							2009-06-26*/
/** \ingroup skClimmisc
 *	This class implements a climatology that is exactly 1.0. Although this
 *	seems pointless it is useful under certain situations. For example, where the 
 *	extinction has been fully specified in the optical properties.
 **/
/*---------------------------------------------------------------------------*/
class skClimatology_Constant : public skClimatology
{
	private:
		double					m_constantvalue;

	public:
								skClimatology_Constant	( );
								skClimatology_Constant	(double value );
		virtual				   ~skClimatology_Constant	( );
		void					SetConstantValue	( double constantvalue ) {m_constantvalue = constantvalue;}

	public:
		virtual	bool			UpdateCache			( const GEODETIC_INSTANT& placeandtime ) override;
		virtual	bool			GetParameter		( const CLIMATOLOGY_HANDLE& species, const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache) override;
		virtual bool			IsSupportedSpecies  ( const CLIMATOLOGY_HANDLE& species ) override;
	//	virtual bool			CreateClone			( skClimatology** clone) const override;
};


/*-----------------------------------------------------------------------------
 *					skClimatology_OneTemperatureAndPressure		 2015- 1- 9*/
/** **/
/*---------------------------------------------------------------------------*/

class skClimatology_OneTemperatureAndPressure: public skClimatology
{
	private:
		double					m_T;
		double					m_P;

	public:		
								skClimatology_OneTemperatureAndPressure	( );
		virtual				   ~skClimatology_OneTemperatureAndPressure	( );
		void					SetTemperature		( double kelvin)			{m_T = kelvin;}
		void					SetPressure			( double P_pa)				{m_P = P_pa;}

		virtual	bool			UpdateCache			( const GEODETIC_INSTANT& placeandtime ) override;
		virtual	bool			GetParameter		( const CLIMATOLOGY_HANDLE& species, const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache) override;
		virtual bool			IsSupportedSpecies  ( const CLIMATOLOGY_HANDLE& species ) override;
};

