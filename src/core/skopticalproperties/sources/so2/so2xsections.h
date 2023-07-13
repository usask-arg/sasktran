// Vandaele et al., JQSRT 59, 171-184 (1998)


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_SO2_Vandaele2009		 2014- 11- 3*/
/**  \ingroup so2skopticalprop
 *
 * http://satellite.mpic.de/spectral_atlas/cross_sections/Sulfur%20compounds/Inorganic%20S-compounds/SO2_VandaeleHermansFally(2009)_358K_227.275-416.658nm.txt
 * 
 * DATAFILE:	SO2_VandaeleHermansFally(2009)_358K_227.275-416.658nm.txt
 * NAME:	sulfur dioxide
 * FORMULA:	SO2
 * AUTHOR(YEAR):	VandaeleHermansFally(2009)
 * T:	358K
 * ?:	227.275-416.658nm
 * 
 * BIBLIOGRAPHY:
 * C. Hermans, A.C. Vandaele, and S. Fally,
 * "Fourier transform measurements of SO2 absorption cross sections:
 * I. Temperature dependence in the 24000-29000 cm-1 (345-420 nm) region,"
 * J. Quant. Spectrosc. Radiat. Transfer 110, 756-765 (2009); DOI: 10.1016/j.jqsrt.2009.01.031
 * 
 * A.C. Vandaele, C. Hermans, and S. Fally,
 * "Fourier transform measurements of SO2 absorption cross sections:
 * II. Temperature dependence in the 29000-44000 cm-1 (227-345 nm) region,"
 * J. Quant. Spectrosc. Radiat. Transfer 110, 2115-2126 (2009); DOI: 10.1016/j.jqsrt.2009.05.006
 * 
 * COMMENTS:	High-resolution absorption measurements using a Fourier transform spectrometer at a resolution of 2 cm-1
 * T = 298, 318, 338, 358 K,
 * p = 0.05 - 199 Torr,
 * spectral interval 23500-44500 cm-1
 * 
 * Data file at intervals of 0.5 cm-1 from the webpage of the Belgian Institute for Space Aeronomy
 * http://spectrolab.aeronomie.be/index.htm
 **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_SO2_Vandaele2009 : public skOpticalProperties_UserDefinedAbsorption,
	                                         public skWavelengthToPSF_TableConstantWavenumber,					// Describes the instrument used to measure this cross-section
											 public skOpticalProperty_AdditionalStateInfo_TemperatureDependent	// Describes what state parameters affetct the cross-section

{
	private:
												skOpticalProperties_SO2_Vandaele2009	( const skOpticalProperties_SO2_Vandaele2009& other );	// dont allow copy constructor
		skOpticalProperties_SO2_Vandaele2009&	operator =								( const skOpticalProperties_SO2_Vandaele2009& other );	// Dont allow assignment operator
		bool									ConfigureEntries						();


	public:
												skOpticalProperties_SO2_Vandaele2009	();
		virtual								   ~skOpticalProperties_SO2_Vandaele2009	() override {}

};



/*---------------------------------------------------------------------------
 *           Class skOpticalProperties_SO2_Freeman1984           2019-05-14 */
/**
 * http://satellite.mpic.de/spectral_atlas/cross_sections/Sulfur%20compounds/Inorganic%20S-compounds/SO2_Freeman(1984)_213K_238.095-239.860nm.txt 
 * DATAFILE:	SO2_Freeman(1984)_213K_238.095-239.860nm.txt
 * NAME:	sulfur dioxide
 * FORMULA:	SO2
 * AUTHOR(YEAR):	Freeman(1984)
 * T:	213K
 * wavelength:	238.095-239.860nm
 * BIBLIOGRAPHY:
 * D.E. Freeman, K. Yoshino, J.R. Esmond, and W.H. Parkinson,
 * "High resolution cross section measurements of SO2 at 213 K in the wavelength region 172 - 240 nm",
 * Planet. Space Sci. 32, 1125-1134 (1984); DOI: 10.1016/0032-0633(84)90139-9
 * COMMENTS:	High-resolution (FWHM = 0.002 nm) absorption measurements
 * File from  CfA (Harvard-Smithsonian Center for Astrophysics) Molecular Databases
 *  http://cfa-www.harvard.edu/amdata/ampdata/amdata.shtml
**/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_SO2_Freeman1984 : public skOpticalProperties_UserDefinedAbsorption,
											public skWavelengthToPSF_TableConstant,								// Describes the instrument used to measure this cross-section
											public skOpticalProperty_AdditionalStateInfo_NotDependent			// Describes what state parameters affetct the cross-section

{
	private:
												skOpticalProperties_SO2_Freeman1984		(const skOpticalProperties_SO2_Freeman1984& other);	// dont allow copy constructor
		skOpticalProperties_SO2_Freeman1984&	operator =								(const skOpticalProperties_SO2_Freeman1984& other);	// Dont allow assignment operator
		bool									ConfigureEntries						();


	public:
												skOpticalProperties_SO2_Freeman1984	();
		virtual								   ~skOpticalProperties_SO2_Freeman1984	() override {}

};



/*---------------------------------------------------------------------------
 *            Class skOpticalProperties_SO2_Rufus2003             2019-05-14 */
/**
 * http://satellite.mpic.de/spectral_atlas/cross_sections/Sulfur%20compounds/Inorganic%20S-compounds/SO2_Rufus(2003)_295K_220.0-325.2nm.txt
 * DATAFILE:	SO2_Rufus(2003)_295K_220.0-325.2nm.txt
 * NAME:	sulfur dioxide
 * FORMULA:	SO2
 * AUTHOR(YEAR):	Rufus(2003)
 * T:	295K
 * wavelength:	220.0-325.2nm
 * BIBLIOGRAPHY:
 * J. Rufus, G. Stark, P.L. Smith, J.C. Pickering, and A.P. Thorne,
 * "High-resolution photoabsorption cross section measurements of SO2, 2: 220 to 325 nm at 295 K",
 * J. Geophys. Res - Planets 108, Art. No. 5011 (2003); DOI: 10.1029/2002JE001931
 * COMMENTS:	High-resolution (0.0005-0.0032 nm) absorption cross sections measured by using a Fourier transform spectrometer
 * Absorption cross sectionsfrom
 * CfA (Harvard-Smithsonian Center for Astrophysics) Molecular Databases
 * http://cfa-www.harvard.edu/amdata/ampdata/amdata.shtml
**/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_SO2_Rufus2003 : public skOpticalProperties_UserDefinedAbsorption,
										  public skWavelengthToPSF_TableConstant,							// Describes the instrument used to measure this cross-section
										  public skOpticalProperty_AdditionalStateInfo_TemperatureDependent	// Describes what state parameters affetct the cross-section

{
	private:
												skOpticalProperties_SO2_Rufus2003	(const skOpticalProperties_SO2_Rufus2003& other);	// dont allow copy constructor
		skOpticalProperties_SO2_Rufus2003&		operator =							(const skOpticalProperties_SO2_Rufus2003& other);	// Dont allow assignment operator
		bool									ConfigureEntries					();


	public:
												skOpticalProperties_SO2_Rufus2003	();
		virtual								   ~skOpticalProperties_SO2_Rufus2003	() override {}

};




/*---------------------------------------------------------------------------
 *           Class skOpticalProperties_SO2_Bogumil2003            2019-05-14 */
/** **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_SO2_Bogumil2003 : public skOpticalProperties_UserDefinedAbsorption,
										    public skWavelengthToPSF_TableConstant,							// Describes the instrument used to measure this cross-section
										    public skOpticalProperty_AdditionalStateInfo_TemperatureDependent	// Describes what state parameters affetct the cross-section

{
	private:
												skOpticalProperties_SO2_Bogumil2003	(const skOpticalProperties_SO2_Bogumil2003& other);	// dont allow copy constructor
		skOpticalProperties_SO2_Bogumil2003&	operator =							(const skOpticalProperties_SO2_Bogumil2003& other);	// Dont allow assignment operator
		bool									ConfigureEntries					();


	public:
												skOpticalProperties_SO2_Bogumil2003	();
		virtual								   ~skOpticalProperties_SO2_Bogumil2003	() override {}

};
