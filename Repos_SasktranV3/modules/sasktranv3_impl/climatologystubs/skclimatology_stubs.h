
/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_BaseEngine		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKClimatology_Stub_Base : public ISKClimatology_Stub
{
	private:
		skClimatology*						m_climatology;

	public:
											ISKClimatology_Stub_Base	( skClimatology* climate);
		virtual 						   ~ISKClimatology_Stub_Base	() override;
		virtual nxUnknown*					RawObjectPointer			() override { return m_climatology;}
		virtual bool						UpdateCache					( const GEODETIC_INSTANT& location ) override;
		virtual bool						GetParameter				( const CLIMATOLOGY_HANDLE& species,  const GEODETIC_INSTANT& location, double* value ) override;
		virtual bool						SetPropertyUserDefined		( const CLIMATOLOGY_HANDLE& species,  double* profile, int numpoints) override;
};


/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_UserDefined					2014-3-7*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKClimatology_Stub_UserDefined : public ISKClimatology_Stub_Base
{
	private:
		skClimatology_UserTableSpline*		m_userdefinedclimatology;
		std::vector<double>					m_currentheightarray;
		bool								m_dologinterpolation;
		bool								m_dopiecewiselinear;
		double								m_currentbadvalue;

	private:
		void								MakeSetFunctions();

	public:
											ISKClimatology_Stub_UserDefined	( skClimatology_UserTableSpline* climate);
		virtual 						   ~ISKClimatology_Stub_UserDefined	() override;
		virtual bool						SetPropertyUserDefined		    ( const CLIMATOLOGY_HANDLE& species,  double* profile, int numheights) override;

};

/*-----------------------------------------------------------------------------
 *					class ISKClimatology_Stub_Constant				2015- 1- 9*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKClimatology_Stub_Constant : public ISKClimatology_Stub_Base
{
	private:
		skClimatology_Constant*				m_ptclimate;

	private:
		void								MakeSetFunctions();

	public:
											ISKClimatology_Stub_Constant( skClimatology_Constant* ptclimate);									
		virtual							   ~ISKClimatology_Stub_Constant();
};

/*-----------------------------------------------------------------------------
 *					class ISKClimatology_Stub_OnePressureTemp		 2015- 1- 9*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKClimatology_Stub_OnePressureTemp : public ISKClimatology_Stub_Base
{
	private:
		skClimatology_OneTemperatureAndPressure*		m_ptclimate;

	private:
		void								MakeSetFunctions();

	public:
											ISKClimatology_Stub_OnePressureTemp( skClimatology_OneTemperatureAndPressure* ptclimate);									
		virtual							   ~ISKClimatology_Stub_OnePressureTemp();
};


/*-----------------------------------------------------------------------------
 *					class ISKClimatology_Stub_Constant				2015- 1- 9*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKClimatology_Stub_LinearCombination : public ISKClimatology_Stub_Base
{
	private:
		skClimatologyLinearCombination*		m_ptclimate;

	private:
		void								MakeSetFunctions();

	public:
											ISKClimatology_Stub_LinearCombination( skClimatologyLinearCombination* ptclimate);									
		virtual							   ~ISKClimatology_Stub_LinearCombination();
};


/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_UserDefined3D		2014-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKClimatology_Stub_UserDefined3D : public ISKClimatology_Stub_Base
{
	private:
		skClimatology_UserDefined3D_LatLonHeight*	m_userdefinedclimatology;
		std::vector<double>							m_currentheightarray;
		std::vector<double>							m_currentlatarray;
		std::vector<double>							m_currentlonarray;
		double										m_currentbadvalue;

	private:
		void								MakeSetFunctions();
	public:
											ISKClimatology_Stub_UserDefined3D	( skClimatology_UserDefined3D_LatLonHeight* climate);
		virtual 						   ~ISKClimatology_Stub_UserDefined3D	() override;
		virtual bool						SetPropertyUserDefined				( const CLIMATOLOGY_HANDLE& species,  double* profile, int numheights) override;
};


/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_OSIRISL2_AEROSOLMODERADIUS_V600		2014-3-31*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKClimatology_Stub_OSIRISL2_AEROSOLMODERADIUS_V600: public ISKClimatology_Stub_Base
{
	private:
		skClimatology_OsirisAerosolModeRadiusV600*	m_climatologyv600;

	private:
		void								MakeSetFunctions();

	public:
											ISKClimatology_Stub_OSIRISL2_AEROSOLMODERADIUS_V600( skClimatology_OsirisAerosolModeRadiusV600* climate);
		virtual							   ~ISKClimatology_Stub_OSIRISL2_AEROSOLMODERADIUS_V600() override;
};



/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_UserDefinedTable		 2016- 9- 26*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKClimatology_Stub_UserDefinedTable : public ISKClimatology_Stub_Base
{
	private:
		CLIMATOLOGY_HANDLE					m_currenthandle;
		skClimatology_UserDefinedTable*		m_userdefinedclim;
		std::vector<double>					m_heights;

	private:
		void								MakeSetFunctions();

	public:
											ISKClimatology_Stub_UserDefinedTable	( skClimatology_UserDefinedTable* clim );
		virtual								~ISKClimatology_Stub_UserDefinedTable	() override { };
		virtual bool						SetPropertyUserDefined					( const CLIMATOLOGY_HANDLE& species,  double* profile, int numheights) override;
};


/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_UserDefinedPlane		 2016- 9- 26*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKClimatology_Stub_UserDefinedPlane : public ISKClimatology_Stub_Base
{
	private:
		skClimatology_UserDefinedPlane*		m_userdefinedclim;
		size_t								m_numheights;
		size_t								m_numangles;
		bool								m_dolog;

	private:
		void								MakeSetFunctions();

	public:
											ISKClimatology_Stub_UserDefinedPlane	( skClimatology_UserDefinedPlane* clim );
		virtual								~ISKClimatology_Stub_UserDefinedPlane	() override { };
		virtual bool						SetPropertyUserDefined					( const CLIMATOLOGY_HANDLE& species,  double* profile, int numheights) override;
};


/*-----------------------------------------------------------------------------
 *					class ISKClimatology_Stub_OnePressureTemp		 2015- 1- 9*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKClimatology_Stub_MSIS : public ISKClimatology_Stub_Base
{
	private:
		skClimatology_MSIS90*		m_msis;

	private:
		void								MakeSetFunctions();

	public:
											ISKClimatology_Stub_MSIS( skClimatology_MSIS90* msis);									
		virtual							   ~ISKClimatology_Stub_MSIS();
};

