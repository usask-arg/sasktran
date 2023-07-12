// Vandaele et al., JQSRT 59, 171-184 (1998)


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_NO2_Vandaele1998		 2014- 11- 3*/
/**  \ingroup no2skopticalprop
 *	Tabulates the absorption cross section of NO2 molecules at medium
 *	resolution (2 cm-1) from 238–1000 nm at 2 temperatures. The
 *	cross-sections were measured by Vandaele et al 1998.
 *
 *
 *	\par Temperature Range
 *	The cross-sections were measured at 2 temperatures covering normal stratospheric and tropospheric ranges. The paper recommends
 *	linear interpolation in temperature. The spectra were measured at the following two temperatures.
 *		-# 220K
 *		-# 294K
 *
 *	\par References
 *	Vandaele A.C., C. Hermans, P.C. Simon, M. Carleer, R. Colin, S. Fally, M.F. Mérienne, A. Jenouvrier, and B. Coquart,
 *	Measurements of the NO2 absorption cross-section from 42000 cm-1 to 10000 cm-1 (238-1000 nm) at 220 K and 294 K, 
 *	Jqsrt, 59, 171-184 (1998)
 *
 **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_NO2_Vandaele1998 : public skOpticalProperties_UserDefinedAbsorption,
	                                               public skWavelengthToPSF_TableConstantWavenumber,					// Describes the instrument used to measure this cross-section
												   public skOpticalProperty_AdditionalStateInfo_TemperatureDependent	// Describes what state parameters affetct the cross-section

{
	private:
												skOpticalProperties_NO2_Vandaele1998	( const skOpticalProperties_NO2_Vandaele1998& other );	// dont allow copy constructor
		skOpticalProperties_NO2_Vandaele1998&	operator =								( const skOpticalProperties_NO2_Vandaele1998& other );	// Dont allow assignment operator
		bool									ConfigureEntries						();


	public:
												skOpticalProperties_NO2_Vandaele1998	();
		virtual								   ~skOpticalProperties_NO2_Vandaele1998	() override {}

};

