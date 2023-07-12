
/*---------------------------------------------------------------------------
 *                  Class SKTRAN_Short_Test_Base                  2020-01-31 */
/** **/
/*---------------------------------------------------------------------------*/

class SKTRAN_Short_Test_Base
{
	private:
		enum class										GeometryType { limb, ground, nadir };
	private:
		GeometryType									m_geometrytype;
		double											m_errortolerance;

	protected:
		SKTRAN_LineOfSightArray_V21						m_linesofsight;
		nxVector										m_sun;
		std::vector<double>								m_wavelen;
		SKTRAN_AtmosphericOpticalState_V21				m_opticalstate;

	private:
		bool											MakeWavelengths				();
		bool											MakeLimbBasedLinesOfSight	();
		bool											MakeGroundBasedLinesOfSight ();
		bool											MakeNadirBasedLinesOfSight  ();
		bool											MakeOpticalState			();
		bool											Initialize					( GeometryType geometrytype );

	protected:
		virtual bool									MakeSpecs					() = 0;
		virtual bool									Scatter						( nx2dArray<double>& radiance, size_t scattorder ) = 0;
		virtual bool									LoadHardCodedValues			( nx2dArray<double>& sk_hardcode ) = 0;
		virtual const char*								Name						() const = 0;

	public:
		bool											RunStandardTest				();
		bool											RunGroundBasedTests			();
		bool											RunNadirBasedTests			();
														SKTRAN_Short_Test_Base		(double errortolerance = 1.0E-5) {m_geometrytype = GeometryType::limb; m_errortolerance = errortolerance;}
		virtual										   ~SKTRAN_Short_Test_Base		( )  {}
};


/*---------------------------------------------------------------------------
 *                      Class SKTRAN_SO_Test                      2020-01-31 */
/** **/
/*---------------------------------------------------------------------------*/

class SKTRAN_Short_Test_SO : public SKTRAN_Short_Test_Base
{
	private:
		virtual bool									Scatter						( nx2dArray<double>& radiance, size_t scattorder ) override;
		virtual bool									MakeSpecs					() override { return true;}
		virtual bool									LoadHardCodedValues			( nx2dArray<double>& sk_hardcode ) override;
		virtual const char*								Name						() const override  { return "SO"; }

	public:
														SKTRAN_Short_Test_SO		(double errortolerance = 1.0E-5):SKTRAN_Short_Test_Base(errortolerance) {}
		virtual										   ~SKTRAN_Short_Test_SO		() {}
};	

/*---------------------------------------------------------------------------
 *                      Class SKTRAN_HR_Test                      2020-01-31 */
/** **/
/*---------------------------------------------------------------------------*/


class SKTRAN_Short_Test_HR : public SKTRAN_Short_Test_Base
{
	private:
		SKTRAN_HR_Specs_User*							m_hrspecs;

	private:
		virtual bool									Scatter						( nx2dArray<double>& radiance, size_t scattorder ) override;
		virtual bool									MakeSpecs					() override;
		virtual bool									LoadHardCodedValues			( nx2dArray<double>& sk_hardcode ) override;
		virtual const char*								Name						() const override  { return "HR"; }


	public:
														SKTRAN_Short_Test_HR		(double errortolerance = 1.0E-5):SKTRAN_Short_Test_Base(errortolerance)  { m_hrspecs = nullptr;}
		virtual										   ~SKTRAN_Short_Test_HR		() { if (m_hrspecs != nullptr) delete m_hrspecs;}

};	

/*---------------------------------------------------------------------------
 *                      Class SKTRAN_MC_Test                      2020-01-31 */
/** **/
/*---------------------------------------------------------------------------*/

class SKTRAN_Short_Test_MC : public SKTRAN_Short_Test_Base
{
	private:
		SKTRAN_Specifications_MC*						m_mcspecs;

	private:
		virtual bool									Scatter						( nx2dArray<double>& radiance, size_t scattorder ) override;
		virtual bool									MakeSpecs					() override;
		virtual bool									LoadHardCodedValues			( nx2dArray<double>& sk_hardcode ) override;
		virtual const char*								Name						() const override  { return "MC"; }

	public:
														SKTRAN_Short_Test_MC		(double errortolerance = 1.0E-5):SKTRAN_Short_Test_Base(errortolerance)  { m_mcspecs = nullptr;}
		virtual										   ~SKTRAN_Short_Test_MC		() { if (m_mcspecs != nullptr) delete m_mcspecs;}

};	




/*---------------------------------------------------------------------------
 *                    Class SKTRAN_Short_Test                     2020-01-31 */
/** **/
/*---------------------------------------------------------------------------*/

class SKTRAN_Short_Test
{
	private:
		double											m_errortolerance;
		bool											m_dolimb;					// if true then do the limb tests
		bool											m_doSO;						// if true then do the limb tests
		bool											m_doMC;						// if true then do the limb tests
		bool											m_doHR;						// if true then do the limb tests
		SKTRAN_Short_Test_MC							m_mc;
		SKTRAN_Short_Test_HR							m_hr;
		SKTRAN_Short_Test_SO							m_so;

	public:
														SKTRAN_Short_Test			();
													   ~SKTRAN_Short_Test			();
		bool											RunTests					(bool doso, bool dohr, bool domc, double errortolerance);
		bool											RunGroundBasedTests			(bool doso, bool dohr, bool domc);
		
};


/*---------------------------------------------------------------------------
 *                      Class SKTRAN_Short_test_Inelastic_MC       2020-04-03 */
 /** **/
 /*---------------------------------------------------------------------------*/

class SKTRAN_Short_Test_Inelastic_MC : public SKTRAN_Short_Test_Base
{
private:
	SKTRAN_Specifications_MC*						m_mcspecs;
	bool											m_optimized;
	bool											m_simultaneous;
	bool											m_ring;

private:
	virtual bool									Scatter(nx2dArray<double>& radiance, size_t scattorder) override;
	virtual bool									MakeSpecs() override;
	virtual bool									LoadHardCodedValues(nx2dArray<double>& sk_hardcode) override;
	virtual const char*								Name() const override; 

public:
													SKTRAN_Short_Test_Inelastic_MC(double errortolerance = 1.0E-5, bool optimized = false, bool simultaneous = false, bool ring = false) : SKTRAN_Short_Test_Base(errortolerance) { m_mcspecs = nullptr; m_optimized = optimized; m_simultaneous = simultaneous; m_ring = ring; }
	virtual										   ~SKTRAN_Short_Test_Inelastic_MC() { if (m_mcspecs != nullptr) delete m_mcspecs; }

};