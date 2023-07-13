#pragma once

/*---------------------------------------------------------------------------
 *           Class  skClimatologyLinearCombination           2020-01-07 */
/** **/
/*---------------------------------------------------------------------------*/

class  skClimatologyLinearCombination : public skClimatology
{
	private:
		std::vector<skClimatology*>			m_climatologies;
		std::vector<double>					m_coefficients;
		std::vector<double>					m_f;
		std::vector<double>					m_h;
		nxPiecewiseLinear					m_altitudecoeffs;


	private:
		void								ReleaseResources						();
		bool								CheckHeightProfile						(const nx1dArray<double>& heightm, const nx1dArray<double>& f  );
		skClimatology&						operator =								( const  skClimatologyLinearCombination& other );  // =delete; Visual Studio 2012 does not like this yet			// Dont allow assignment
											 skClimatologyLinearCombination			( const  skClimatologyLinearCombination& other );  // =delete;			// Dont allow copy constructor

protected:
		virtual	 bool						SetComboCoefficients					( const GEODETIC_INSTANT& pt);

	public:
											 skClimatologyLinearCombination	();
		virtual							   ~ skClimatologyLinearCombination	();
		bool								SetFirstClimatology						( skClimatology* optprop);
		bool								SetSecondClimatology					( skClimatology* optprop);
		bool								SetHeightProfileCoeffsOfFirstProperty	( const nx1dArray<double>& heightm, const nx1dArray<double>& f );

	public:
		virtual	bool						UpdateCache								( const GEODETIC_INSTANT& placeandtime) override;
		virtual	bool						GetParameter							( const CLIMATOLOGY_HANDLE& species,  const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache) override;
		virtual bool						IsSupportedSpecies						( const CLIMATOLOGY_HANDLE& species ) override;

};


