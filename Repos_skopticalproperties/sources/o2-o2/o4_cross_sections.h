
/*---------------------------------------------------------------------------
*            Class skOpticalProperties_O4_Fally2000            2019-02-14 */
/** **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_O4_Fally2000 :	public skOpticalProperties_UserDefinedAbsorption,
											public skWavelengthToPSF_TableConstantWavenumber,				// Describes the instrument used to measure this cross-section
											public skOpticalProperty_AdditionalStateInfo_NotDependent		// Describes what state parameters affetct the cross-section

{
	private:
		nx1dArray<double>	m_wavenm;
		nx1dArray<double>	m_o4xsc;

	private:
							skOpticalProperties_O4_Fally2000					(const skOpticalProperties_O4_Fally2000& other);	// Dont allow copy constructor
							skOpticalProperties_O4_Fally2000&	operator =		(const skOpticalProperties_O4_Fally2000& other);	// Dont allow assignment operator

	public:
							skOpticalProperties_O4_Fally2000					();
		virtual			   ~skOpticalProperties_O4_Fally2000					() override {}
	};



/*---------------------------------------------------------------------------
 *            Class skOpticalProperties_O4_Thalman2013            2019-12-05 */
/** **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_O4_Thalman2013 :	public skOpticalProperties_UserDefinedAbsorption,
											public skWavelengthToPSF_TableConstant,									// Describes the instrument used to measure this cross-section
											public skOpticalProperty_AdditionalStateInfo_TemperatureDependent		// Describes what state parameters affetct the cross-section

{

	private:
							skOpticalProperties_O4_Thalman2013					(const skOpticalProperties_O4_Thalman2013& other);	// Dont allow copy constructor
							skOpticalProperties_O4_Thalman2013&	operator =		(const skOpticalProperties_O4_Thalman2013& other);	// Dont allow assignment operator

	public:
							skOpticalProperties_O4_Thalman2013					();
	virtual				   ~skOpticalProperties_O4_Thalman2013					() override {}
};


/*---------------------------------------------------------------------------
 *            Class skOpticalProperties_O4_HitranEntry_TempDependent  2020-01-02 */
/** **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_O4_HitranEntry_TempDependent:	public skOpticalProperties_UserDefinedAbsorption,
														public skWavelengthToPSF_TableConstantWavenumber,						// Describes the instrument used to measure this cross-section
														public skOpticalProperty_AdditionalStateInfo_TemperatureDependent		// Describes what state parameters affetct the cross-section
{

	public:
								skOpticalProperties_O4_HitranEntry_TempDependent	(int regionid);
		virtual				   ~skOpticalProperties_O4_HitranEntry_TempDependent	() override {}
		bool					ConfigureAsRegion1									();
		bool					ConfigureAsRegion7									();
		bool					ConfigureAsRegion8									();

};


/*---------------------------------------------------------------------------
 *            Class skOpticalProperties_O4_HitranEntry_NotDependent            2020-01-02 */
/** **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_O4_HitranEntry_NotDependent:	public skOpticalProperties_UserDefinedAbsorption,
														public skWavelengthToPSF_TableConstantWavenumber,									// Describes the instrument used to measure this cross-section
														public skOpticalProperty_AdditionalStateInfo_NotDependent		// Describes what state parameters affetct the cross-section
{
	private:
		nx1dArray<double>		m_nm;
		nx1dArray<double>		m_xs;

	private:
		bool					AddWaveNumberEntry( double t, double* nm, int nmstride, double *xs, int xsstride, int npts );

	public:
								skOpticalProperties_O4_HitranEntry_NotDependent	(int regionid);
		virtual				   ~skOpticalProperties_O4_HitranEntry_NotDependent	() override {}
		bool					ConfigureAsRegion2								();
		bool					ConfigureAsRegion3								();
		bool					ConfigureAsRegion4								();
		bool					ConfigureAsRegion5								();
		bool					ConfigureAsRegion6								();

};

/*---------------------------------------------------------------------------
 *            Class skOpticalProperties_O4_Hitran2016            2019-12-05 */
/** **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_O4_Hitran2016 :	public skOpticalProperties_ListEntries
{
	private:
							skOpticalProperties_O4_Hitran2016					(const skOpticalProperties_O4_Hitran2016& other);	// Dont allow copy constructor
							skOpticalProperties_O4_Hitran2016&	operator =		(const skOpticalProperties_O4_Hitran2016& other);	// Dont allow assignment operator

	public:
							skOpticalProperties_O4_Hitran2016					();
	virtual				   ~skOpticalProperties_O4_Hitran2016					() override;
};
