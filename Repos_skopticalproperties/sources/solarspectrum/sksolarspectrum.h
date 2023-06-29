/*-----------------------------------------------------------------------------
 *					skSolarSpectrum									2015- 1- 6*/
/** Base classes for Solar Spectrum classes **/
/*---------------------------------------------------------------------------*/

class skSolarSpectrum :public nxUnknown
{
	private:
		double						m_solardistanceinAU;

	private:
		double						AUCorrection();

	public:
									skSolarSpectrum		();
		virtual					   ~skSolarSpectrum		();
		bool						SetSolarDistanceFromMjd		( double mjd );		
		bool						SetSolarDistanceFromAU		( double au );		
		double						Irradiance					( double wavelen_nm_vacuum );
		virtual double				IrradianceAt1AU				( double wavelen_nm_vacuum ) = 0;
		virtual double				MinValidWavelength			() const = 0;
		virtual double				MaxValidWavelength			() const = 0;
		virtual double				NanometerResolutionFWHM		(double wavelen_nm_vacuum) const = 0;
		virtual double				SampleSpacing				(double wavelen_nm_vacuum) const = 0;
};

/*-----------------------------------------------------------------------------
 *					skSolarSpectrum_TabulatedWavelength		 2015- 1- 6*/
/** **/
/*---------------------------------------------------------------------------*/

class skSolarSpectrum_TabulatedWavelength : public skSolarSpectrum
{
	private:
//		nx2dArray<double>			m_irradiancetable;
		nx1dArray< double>			m_wavelen_nm_invacuum;
		nx1dArray< double>			m_irradiance;

	private:
		double						AUCorrection();

	public:
									skSolarSpectrum_TabulatedWavelength		();
		virtual					   ~skSolarSpectrum_TabulatedWavelength		();
		bool						AttachToTable				( const double* table2d, size_t numpoints );
		bool						AttachToTable				( const double* flux, const double* wavelength, size_t numelements );
		virtual double				IrradianceAt1AU				( double wavelen_nm_vacuum ) override;
		virtual double				MinValidWavelength			() const override;
		virtual double				MaxValidWavelength			() const override;
};

/*-----------------------------------------------------------------------------
 *					skSolarSpectrum_SAO2010							2015- 1- 6*/
/** Implements the SAO2010 solar spectrum reported by Chance & Kurucz in 2010.
 *	The spectrum is a high quality. high resolution solar spectrum going from
 *	200 nm to 1 micron at 0.01 nm spacing and a resolution of 0.04 nm FWHM. 
 *
 *	@article{Chance20101289,
 *	title = "An improved high-resolution solar reference spectrum for earth's atmosphere measurements in the ultraviolet, visible, and near infrared ",
 *	journal = "Journal of Quantitative Spectroscopy and Radiative Transfer ",
 *	volume = "111",
 *	number = "9",
 *	pages = "1289 - 1295",
 *	year = "2010",
 *	note = "Special Issue Dedicated to Laurence S. Rothman on the Occasion of his 70th Birthday. ",
 *	issn = "0022-4073",
 *	doi = "http://dx.doi.org/10.1016/j.jqsrt.2010.01.036",
 *	url = "http://www.sciencedirect.com/science/article/pii/S0022407310000610",
 *	author = "K. Chance and R.L. Kurucz",
 *	keywords = "Solar spectrum",
 *	keywords = "Remote sensing",
 *	keywords = "Atmospheric spectroscopy "
 *	}
 *
 *	http://www.cfa.harvard.edu/atmosphere/
 **/
/*---------------------------------------------------------------------------*/

class skSolarSpectrum_SAO2010 : public skSolarSpectrum_TabulatedWavelength,
							   public skWavelengthToPSF_TableConstant
{
	private:
		void						MakeGetFunctions			();
	public:
									skSolarSpectrum_SAO2010		();
		virtual					   ~skSolarSpectrum_SAO2010		();
		virtual double				NanometerResolutionFWHM		(double wavelen_nm_vacuum) const override;
		virtual double				SampleSpacing				(double wavelen_nm_vacuum) const override;

};


/*-----------------------------------------------------------------------------
 *					skSolarSpectrum_FontelaUVIS3Micron		 2015- 1- 6*/
/** Implements the solar spectrum presented by Fontela et al 2001. We use the 200 nm to 3 micron table
 *	at 0.02 nm spacing with 0.1 nm resolution., file -> NUVisIr3irrad1001Lowion00.conv1a.fits.
 *
 *	See Fontela 2011 paper:
 *	@article {JGRD:JGRD17289,
 *	author = {Fontenla, J. M. and Harder, J. and Livingston, W. and Snow, M. and Woods, T.},
 *	title = {High-resolution solar spectral irradiance from extreme ultraviolet to far infrared},
 *	journal = {Journal of Geophysical Research: Atmospheres},
 *	volume = {116},
 *	number = {D20},
 *	issn = {2156-2202},
 *	url = {http://dx.doi.org/10.1029/2011JD016032},
 *	doi = {10.1029/2011JD016032},
 *	pages = {n/a--n/a},
 *	keywords = {EUV variability, UV variability, infrared variability, solar irradiance spectrum, solar spectral irradiance, visible variability},
 *	year = {2011},
 *	}
 *
 *	http://www.digidyna.com/Results2010/spectra/irradiance/index_spectra_irradiance.html
 *	http://www.digidyna.com/Results2010/spectra/irradiance/NUVisIr3irrad1001Lowion00.conv1a.fits
**/
/*---------------------------------------------------------------------------*/

class skSolarSpectrum_FontelaUVIS3Micron : public skSolarSpectrum_TabulatedWavelength,
							               public skWavelengthToPSF_TableConstant
{
	public:
									skSolarSpectrum_FontelaUVIS3Micron		();
		virtual					   ~skSolarSpectrum_FontelaUVIS3Micron		();
		virtual double				NanometerResolutionFWHM					(double wavelen_nm_vacuum) const override;
		virtual double				SampleSpacing							(double wavelen_nm_vacuum) const override;

};

/*-----------------------------------------------------------------------------
 *					skSolarSpectrum_FontelaUVIS100Micron		 2017- 8- 23*/
/** Implements the solar spectrum presented by Fontela et al 2001. We use the 200 nm to 100 micron table
 *	at 0.2 nm spacing with 1 nm resolution., file -> NUVisIr100irrad1001Lowion00.conv1nm.fits.
 *
 *	See Fontela 2011 paper:
 *	@article {JGRD:JGRD17289,
 *	author = {Fontenla, J. M. and Harder, J. and Livingston, W. and Snow, M. and Woods, T.},
 *	title = {High-resolution solar spectral irradiance from extreme ultraviolet to far infrared},
 *	journal = {Journal of Geophysical Research: Atmospheres},
 *	volume = {116},
 *	number = {D20},
 *	issn = {2156-2202},
 *	url = {http://dx.doi.org/10.1029/2011JD016032},
 *	doi = {10.1029/2011JD016032},
 *	pages = {n/a--n/a},
 *	keywords = {EUV variability, UV variability, infrared variability, solar irradiance spectrum, solar spectral irradiance, visible variability},
 *	year = {2011},
 *	}
 *
 *	http://www.digidyna.com/Results2010/spectra/irradiance/index_spectra_irradiance.html
 *	http://www.digidyna.com/Results2010/spectra/irradiance/NUVisIr100irrad1001Lowion00.conv1nm.fits
**/
/*---------------------------------------------------------------------------*/

class skSolarSpectrum_FontelaUVIS100Micron : public skSolarSpectrum_TabulatedWavelength,
											 public skWavelengthToPSF_TableConstant
{
	public:
									skSolarSpectrum_FontelaUVIS100Micron	();
		virtual					   ~skSolarSpectrum_FontelaUVIS100Micron	();
		virtual double				NanometerResolutionFWHM					(double wavelen_nm_vacuum) const override;
		virtual double				SampleSpacing							(double wavelen_nm_vacuum) const override;

};

