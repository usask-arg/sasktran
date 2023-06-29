#include <skopticalproperties21.h>
#include "sources/sasktranif_opticalimpl/skoptprop_stubs.h"

/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_Base::ISKOpticalProperty_Stub_Base		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

ISKOpticalProperty_Stub_Base::ISKOpticalProperty_Stub_Base	( skOpticalProperties* optprop)
{
	m_opticalproperty = optprop;
	if (m_opticalproperty != NULL) m_opticalproperty->AddRef();
}


/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_Base::~ISKOpticalProperty_Stub_Base		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

ISKOpticalProperty_Stub_Base::~ISKOpticalProperty_Stub_Base	()
{
	if (m_opticalproperty != NULL) m_opticalproperty->Release();
}


/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_Base::SkOpticalPropertyPointer		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

//nxUnknown*	ISKOpticalProperty_Stub_Base::RawObjectPointer()
//{ 
//	return m_opticalproperty;
//}

/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_Base::SetAtmosphericState		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKOpticalProperty_Stub_Base::SetAtmosphericState( ISKClimatology_Stub* atmosphere )
{
	skClimatology*	ptr;

	ptr = dynamic_cast<skClimatology*>(atmosphere->RawObjectPointer());
	return (m_opticalproperty != nullptr) && m_opticalproperty->SetAtmosphericState( ptr );
}


/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_Base::SetLocation		 2015- 11- 17*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKOpticalProperty_Stub_Base::SetLocation(const GEODETIC_INSTANT& pt )
{
	bool haschanged;
	return (m_opticalproperty != nullptr) && m_opticalproperty->SetLocation( pt, &haschanged );
}


/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_Base::InternalClimatology_UpdateCache		 2014- 10- 27*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKOpticalProperty_Stub_Base::InternalClimatology_UpdateCache( const GEODETIC_INSTANT& pt )
{
	return (m_opticalproperty != nullptr) && m_opticalproperty->InternalClimatology_UpdateCache(pt);
}


/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_Base::GetExtinctionPerCm		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKOpticalProperty_Stub_Base::CalculateCrossSections( double wavenumber, double* absxs, double* extxs, double* scattxs)
{
	bool	ok;

	ok = (m_opticalproperty != nullptr) && m_opticalproperty->CalculateCrossSections( wavenumber, absxs, extxs, scattxs);
	return ok;
}



/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_Base::CalculateCrossSectionsArray		2014-3-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKOpticalProperty_Stub_Base::CalculateCrossSectionsArray( const double * wavenumber, int numwavenumber, double *absxs,  double* extxs, double* scattxs)
{
	bool	ok;
	ok =	(m_opticalproperty != nullptr) && m_opticalproperty->CalculateCrossSectionsArray( wavenumber, numwavenumber, absxs, extxs, scattxs);
	return ok;
}

bool ISKOpticalProperty_Stub_Base::CalculatePhaseMatrix(const double * wavenumber, const double * cosscatterangle, double* phasematrix)
{
	skRTPhaseMatrix phase;
	bool ok = m_opticalproperty->CalculatePhaseMatrix(*wavenumber, *cosscatterangle, &phase);

	for (int i = 0; i < 16; i++)
	{
		phasematrix[i] = *(phase.ArrayBasePtr() + i);
	}

	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_Base::AddUserDefined		2014-4-2*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKOpticalProperty_Stub_Base::AddUserDefined( double temperature,  double* wavelen_nm, int numwave, double* crosssection, int numcross)
{
	nxLog::Record(NXLOG_WARNING,"ISKOpticalProperty_Stub_Base::SetPropertyObject, this optical property does not support AddUserDefined functionality");
	return false;
}

/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_Base::AddUserDefinedPressure		2021-08-04*/
 /** **/
 /*---------------------------------------------------------------------------*/

bool ISKOpticalProperty_Stub_Base::AddUserDefinedPressure(double* pressure, int numpressure, double* temperature, int numtemperature, double* wavelen_nm, int numwavel, double* crosssection, int numcross, double broadnervmr)
{
	nxLog::Record(NXLOG_WARNING, "ISKOpticalProperty_Stub_Base, this optical property does not support AddUserDefinedPressure functionality");
	return false;
}


/*-----------------------------------------------------------------------------
 *					:ISKOpticalProperty_Stub_Hitran		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

ISKOpticalProperty_Stub_Hitran::ISKOpticalProperty_Stub_Hitran	( skOpticalProperties_HitranChemical* optprop)
	                           :ISKOpticalProperty_Stub_Base (optprop)
{
	m_hitranoptprop = optprop;
	MakeSetPropertyFunctions();
}


/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_Hitran::~ISKOpticalProperty_Stub_Hitran		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

ISKOpticalProperty_Stub_Hitran::~ISKOpticalProperty_Stub_Hitran	()
{
}

/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_Hitran::Set		2014-3-6*/
/** SetWavenumberRange (lowwavenum, highwavenum)
 **/
/*---------------------------------------------------------------------------*/

void ISKOpticalProperty_Stub_Hitran::MakeSetPropertyFunctions()
{

	AddSetVectorFunction( "setwavenumberrange",
		[&, this] ( const double* value, int numpoints)
		{
			bool ok;
			ok =      ( numpoints == 2); 
			ok = ok && m_hitranoptprop->SetWavenumberRange( value[0], value[1] );
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "ISKOpticalProperty_Stub_Hitran::SetProperty(setwavenumberrange). There were errors setting the array properties.");
			}
			return ok;
		}
	);

	AddSetVectorFunction( "enablecachedcrosssections",
		[&, this] ( const double* value, int numpoints)
		{
			bool ok;
			double *	wavenumptr = (double *)(intptr_t)value;
			ok = m_hitranoptprop->EnableCachedCrossSections( wavenumptr, numpoints );
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "ISKOpticalProperty_Stub_Hitran::SetProperty(enablecachedcrosssections). There were errors setting the array properties.");
			}
			return ok;
		}
	);


	AddSetScalarFunction( "setlinetolerance",
		[&, this](double value )
		{
			bool	ok;
			ok = m_hitranoptprop->SetLineTolerance( value );
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKOpticalProperty_Stub_Hitran::SetProperty(setlinetolerance), there were errors setting the line tolerance to %15.8e", (double)value);
			}
			return ok;
		}
	);

	AddSetScalarFunction( "setmaxlinestrength",
		[&, this](double value )
		{
			bool ok;
			ok = m_hitranoptprop->SetUserDefinedMaxLineStrength( value );
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKOpticalProperty_Stub_Hitran::SetProperty(setmaxlinestrength), there were errors setting the micro-window maximum line strength line to %15.8e", (double)value);
			}
			return ok;
		}
	);

	AddSetScalarFunction( "setisotopefilter",
		[&, this](double value )
		{
			bool	ok;
			size_t	filterid;

			filterid = (size_t)value;
			ok = m_hitranoptprop->SetIsotopeFilter( filterid );
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "ISKOpticalProperty_Stub_Hitran::SetProperty(setisotopefilter). There were errors setting the property of SetIsotopeFilter to %15.8e",(double)value);
			}
			return ok;
		}
	);
	AddSetStringFunction( "set_lower_state_global_quanta_filter",
		[&, this](const char *value )
		{
			bool				ok;
				
			ok  = m_hitranoptprop->SetLowerStateGlobalQuantaFilter( std::string(value) );

			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKOpticalProperty_Stub_Hitran::SetProperty(set_lower_state_global_quanta_filter) failed.");
			}
			return ok;
		}
	);
	AddSetStringFunction( "set_upper_state_global_quanta_filter",
		[&, this](const char *value )
		{
			bool				ok;
				
			ok  = m_hitranoptprop->SetUpperStateGlobalQuantaFilter( std::string(value) );

			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKOpticalProperty_Stub_Hitran::SetProperty(set_upper_state_global_quanta_filter) failed.");
			}
			return ok;
		}
	);

	AddSetScalarFunction( "setmicrowindowmargin",
		[&, this](double value )
		{
			bool	ok;

			ok = m_hitranoptprop->SetMicroWindowMargin( value );
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "ISKOpticalProperty_Stub_Hitran::SetProperty(setmicrowindowmargin). There were errors setting the property of SetMicroWindowMargin to %e",(double)value);
			}
			return ok;
		}
	);
	AddSetObjectFunction( "set_self-broadening_climatology", 
		[&, this]( nxUnknown* object )
		{
			bool			ok = false;
			skClimatology*	climate;

			climate = dynamic_cast<skClimatology*>(object);
			ok      = m_hitranoptprop->SetSelfBroadeningClimatology( climate );
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKOpticalProperty_Stub_Hitran::SetProperty( set_self-broadening_climatology ) failed");
			}
			return ok;
		}
	);
	AddSetStringFunction( "set_self-broadening_climatology_handle",
		[&, this](const char *value )
		{
			bool				ok;
			CLIMATOLOGY_HANDLE*	handleptr;
				
			handleptr = FindGlobalClimatologyHandle ( value );					// MAke sure this handle already exists. If it does not then ** IN THE FUTURE *** we should perhaps create the entry
			ok  = (handleptr != nullptr);
			ok  = ok && m_hitranoptprop->SetSelfBroadeningClimatologyHandle( *handleptr );

			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKOpticalProperty_Stub_Hitran::SetProperty(set_self-broadening_climatology_handle) failed as <%s> was not a recognised existing climatology handle.", (const char*)value);
			}
			return ok;
		}
	);
	AddSetObjectFunction( "set_atmospheric_state_climatology", 
		[&, this]( nxUnknown* object )
		{
			bool			ok = false;
			skClimatology*	climate;

			climate = dynamic_cast<skClimatology*>(object);
			ok      = m_hitranoptprop->SetAtmosphericState( climate );
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKOpticalProperty_Stub_Hitran::SetProperty( set_atmospheric_state_climatology ) failed");
			}
			return ok;
		}
	);

}


/*-----------------------------------------------------------------------------
 *					:ISKOpticalProperty_Stub_Aerosol		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

ISKOpticalProperty_Stub_Aerosol::ISKOpticalProperty_Stub_Aerosol( skOpticalProperties_AerosolProfile* aerosol_optprop)
	                           :ISKOpticalProperty_Stub_Base (aerosol_optprop)
{
	m_aerosol_optprop = aerosol_optprop;
	MakeSetPropertyFunctions();
}


/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_Aerosol::~ISKOpticalProperty_Stub_Aerosol		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

ISKOpticalProperty_Stub_Aerosol::~ISKOpticalProperty_Stub_Aerosol	()
{
}


/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_Aerosol::Set		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/
void ISKOpticalProperty_Stub_Aerosol::MakeSetPropertyFunctions()
{

	AddSetObjectFunction( "setparticlesizeclimatology", 
		[&, this]( nxUnknown* object )
		{
			bool			ok = false;
			skClimatology*	climate;

			climate = dynamic_cast<skClimatology*>(object);
			ok      = m_aerosol_optprop->SetParticleSizeClimatology(climate);
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKOpticalProperty_Stub_Aerosol::SetPropertyObject, Error in property SetParticleSizeClimatology");
			}
			return ok;
		}
	);

	AddSetVectorFunction( "refractiveindexdata",
		[&, this] ( const double* value, int numpoints)
		{
			bool ok = true;
			m_refrac_data.resize(numpoints);
			for (int i = 0; i < numpoints; i++)
			{
				m_refrac_data[i] = value[i];
			}

			skRTRefractiveIndex_Tabulated* refrac = new skRTRefractiveIndex_Tabulated();
			ok = ok && refrac->AttachToStatic(&m_refrac_data[0], m_refrac_data.size());

			m_aerosol_optprop->SetRefractiveIndex(refrac);
			return ok;
		}
	);

	AddSetScalarFunction( "heightdependent",
		[&, this](double value)
		{
			bool dependent = value > 0.5;
			m_aerosol_optprop->SetHeightDependent(dependent);
			return true;
		}
	);

}

/*-----------------------------------------------------------------------------
 *					:ISKOpticalProperty_Stub_Baum		 2016- 9- 23*/
/** **/
/*---------------------------------------------------------------------------*/

ISKOpticalProperty_Stub_Baum::ISKOpticalProperty_Stub_Baum( skOpticalProperties_BaumIceCrystals2014* baum )
	:ISKOpticalProperty_Stub_Base(baum)
{
	m_baum_optprop = baum;
	m_effectivesize = new skClimatology_Constant(1.0);
	m_effectivesize->AddRef();
	m_baum_optprop->SetEffectiveSizeClimatology( m_effectivesize );
	MakeSetPropertyFunctions();
}

ISKOpticalProperty_Stub_Baum::~ISKOpticalProperty_Stub_Baum()
{
	m_effectivesize->Release();
}

/*---------------------------------------------------------------------------
 *     ISKOpticalProperty_Stub_Baum::MakeSetPropertyFunctions     2019-10-02 */
/** **/
/*---------------------------------------------------------------------------*/

void ISKOpticalProperty_Stub_Baum::MakeSetPropertyFunctions() 
{

	AddSetScalarFunction( "effectivesizemicrons",
		[&, this](double value )
		{
			m_effectivesize->SetConstantValue( value );
			return value > 0;
		}
	);

	
	AddSetScalarFunction( "usedeltaeddington",
		[&, this](double value )
		{
			m_baum_optprop->SetUseDeltaEddingtonApproximation( 0.0 != value );
			return true;
		}
	);
}

/*-----------------------------------------------------------------------------
 *					:ISKOpticalProperty_Stub_UserDefined		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

ISKOpticalProperty_Stub_UserDefined::ISKOpticalProperty_Stub_UserDefined	( skOpticalProperties_UserDefinedAbsorption* useroptprop)
	                           :ISKOpticalProperty_Stub_Base (useroptprop)
{
	 m_useroptprop = useroptprop;
	 MakeSetPropertyFunctions();
}


/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_UserDefined::~ISKOpticalProperty_Stub_UserDefined		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

ISKOpticalProperty_Stub_UserDefined::~ISKOpticalProperty_Stub_UserDefined	()
{
}


/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_UserDefined::AddUserDefined		2014-4-2*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKOpticalProperty_Stub_UserDefined::AddUserDefined( double temperature,  double* wavelen_nm, int numwave, double* crosssection, int numcross)
{
	nx1dArray<double>	w;
	nx1dArray<double>	xs;
	bool				ok;

	ok =        w.Attach( numwave, wavelen_nm );
	ok = ok && xs.Attach( numcross, crosssection );

	ok = ok && m_useroptprop->AddUserEntry( temperature, w, xs);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"ISKOpticalProperty_Stub_UserDefined::AddUserDefined, There was an error adding the user defined cross-section");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_UserDefined::SetPropertyScalar		2014-4-2*/
/** **/
/*---------------------------------------------------------------------------*/

void ISKOpticalProperty_Stub_UserDefined::MakeSetPropertyFunctions( )
{
	AddSetScalarFunction( "setisscatterer",
		[&, this](double value )
		{
			m_useroptprop->SetIsScatterer( value != 0.0 );
			return true;
		}
	);

	AddSetScalarFunction( "settemperature",
		[&, this](double value )
		{
			m_useroptprop->Set_Temperature( value );
			return true;
		}
	);

	AddSetScalarFunction( "wavelengthtruncation",
		[&, this](double value )
		{
			bool trunc;
			if( value > 0.5 )
				trunc = true;
			else
				trunc = false;
			m_useroptprop->SetQuietWavelengthTruncation( trunc );
			return true;
		}
	);
}

ISKOpticalProperty_Stub_UserDefinedPressure::ISKOpticalProperty_Stub_UserDefinedPressure(skOpticalProperties_UserDefinedAbsorptionPressure* useroptprop) :ISKOpticalProperty_Stub_Base(useroptprop)
{
	m_useroptprop = useroptprop;

	MakeSetPropertyFunctions();
}

ISKOpticalProperty_Stub_UserDefinedPressure::~ISKOpticalProperty_Stub_UserDefinedPressure()
{
}

void ISKOpticalProperty_Stub_UserDefinedPressure::MakeSetPropertyFunctions()
{
	AddSetObjectFunction("broadenerclimatology",
		[&, this](nxUnknown* climatology)
		{
			auto castedClimatology = dynamic_cast<skClimatology*>(climatology);

			m_useroptprop->SetBroadener(castedClimatology);

			return true;
		}
	);

	AddSetStringFunction("broadenerhandle",
		[&, this](const char* str)
		{
			auto handle = FindGlobalClimatologyHandle(str, false);

			m_useroptprop->SetBroadenerHandle(*handle);

			return true;
		}
	);
}

bool ISKOpticalProperty_Stub_UserDefinedPressure::AddUserDefinedPressure(double* pressure, int numpressure, double* temperature, int numtemperature, double* wavelen_nm, int numwavel, double* crosssection, int numcross, double broadnervmr)
{
	size_t nx = numwavel;
	size_t ny = numtemperature / numpressure;
	size_t nz = numpressure;

	nx1dArray<double> pressureArray(nz, pressure);
	nx2dArray<double> temperatureArray(ny, nz, temperature);
	nx3dArray<double> xsArray(nx, ny, nz, crosssection);
	nx1dArray<double> wavelArray(nx, wavelen_nm);

	return m_useroptprop->AddEntry(xsArray, temperatureArray, pressureArray, wavelArray, broadnervmr);
}


ISKOpticalProperty_Stub_UserDefinedScatterConstantHeight::ISKOpticalProperty_Stub_UserDefinedScatterConstantHeight(skOpticalProperties_UserDefinedScatterConstantHeight* useroptprop) :ISKOpticalProperty_Stub_Base(useroptprop) {
	m_useroptprop = useroptprop;
	MakeSetPropertyFunctions();
}

ISKOpticalProperty_Stub_UserDefinedScatterConstantHeight::~ISKOpticalProperty_Stub_UserDefinedScatterConstantHeight()
{

}

void ISKOpticalProperty_Stub_UserDefinedScatterConstantHeight::MakeSetPropertyFunctions() {
	AddSetVectorFunction("wavelengths",
		[&, this](const double* values, int nump) {
			m_numwavelengths = nump;

			nx1dArray<double> wavel(nump, const_cast<double*>(values));

			m_useroptprop->SetWavelengths(wavel);

			return true;
		}
	);

	AddSetVectorFunction("xs_scat",
		[&, this](const double* values, int nump) {
			nx1dArray<double> xs(nump, const_cast<double*>(values));

			m_useroptprop->SetScatteringCrossSection(xs);

			return true;
		}
	);

	AddSetVectorFunction("xs_abs",
		[&, this](const double* values, int nump) {
			nx1dArray<double> xs(nump, const_cast<double*>(values));

			m_useroptprop->SetAbsorptionCrossSection(xs);

			return true;
		}
	);

	AddSetVectorFunction("legendremoments",
		[&, this](const double* values, int nump) {
			size_t numx = m_numwavelengths;
			size_t numy = nump / numx;

			nx2dArray<double> lm(numx, numy, const_cast<double*>(values));

			m_useroptprop->SetLegendreMoments(lm);

			return true;
		}
	);

	AddSetVectorFunction("legendremomentsa1",
		[&, this](const double* values, int nump) {
			size_t numx = m_numwavelengths;
			size_t numy = nump / numx;

			nx2dArray<double> lm(numx, numy, const_cast<double*>(values));

			m_useroptprop->SetLegendrea1(lm);

			return true;
		}
	);

	AddSetVectorFunction("legendremomentsa2",
		[&, this](const double* values, int nump) {
			size_t numx = m_numwavelengths;
			size_t numy = nump / numx;

			nx2dArray<double> lm(numx, numy, const_cast<double*>(values));

			m_useroptprop->SetLegendrea2(lm);

			return true;
		}
	);

	AddSetVectorFunction("legendremomentsa3",
		[&, this](const double* values, int nump) {
			size_t numx = m_numwavelengths;
			size_t numy = nump / numx;

			nx2dArray<double> lm(numx, numy, const_cast<double*>(values));

			m_useroptprop->SetLegendrea3(lm);

			return true;
		}
	);

	AddSetVectorFunction("legendremomentsa4",
		[&, this](const double* values, int nump) {
			size_t numx = m_numwavelengths;
			size_t numy = nump / numx;

			nx2dArray<double> lm(numx, numy, const_cast<double*>(values));

			m_useroptprop->SetLegendrea4(lm);

			return true;
		}
	);

	AddSetVectorFunction("legendremomentsb1",
		[&, this](const double* values, int nump) {
			size_t numx = m_numwavelengths;
			size_t numy = nump / numx;

			nx2dArray<double> lm(numx, numy, const_cast<double*>(values));

			m_useroptprop->SetLegendreb1(lm);

			return true;
		}
	);

	AddSetVectorFunction("legendremomentsb2",
		[&, this](const double* values, int nump) {
			size_t numx = m_numwavelengths;
			size_t numy = nump / numx;

			nx2dArray<double> lm(numx, numy, const_cast<double*>(values));

			m_useroptprop->SetLegendreb2(lm);

			return true;
		}
	);

}



/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_ConvolvedFixedFWHM		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

ISKOpticalProperty_Stub_ConvolvedFixedFWHM::ISKOpticalProperty_Stub_ConvolvedFixedFWHM	( skOpticalProperties_ConvolvedDiscreteWavelenCachedStateFixedFWHM* useroptprop)
	                                       :ISKOpticalProperty_Stub_Base (useroptprop)
{

	m_defaultclimatology.SetStatic();
	m_defaultclimatology.AddRef();
	m_defaultclimatology.SetTemperature( 293.15 );
	m_defaultclimatology.SetPressure   ( 101325.0 );

	m_convolvedoptprop = useroptprop;
	m_default_measurementdetails.SetInstrumentPointSpacing( 0.0001  );			// Use a really narrow almost continuous point spacing
	m_default_measurementdetails.SetInstrumentPSF_FWHM    ( 0.00000 );			// Assume the instrument has infinite resolution
	MakeSetPropertyFunctions();
}

/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_ConvolvedFixedFWHM::~ISKOpticalProperty_Stub_ConvolvedFixedFWHM		 2015- 1- 9*/
/** **/
/*---------------------------------------------------------------------------*/

ISKOpticalProperty_Stub_ConvolvedFixedFWHM::~ISKOpticalProperty_Stub_ConvolvedFixedFWHM	()
{
	m_defaultclimatology.Release();
}

/*---------------------------------------------------------------------------
 * ISKOpticalProperty_Stub_ConvolvedFixedFWHM::MakeSetPropertyFunctions2019-10-02 */
/** **/
/*---------------------------------------------------------------------------*/

void ISKOpticalProperty_Stub_ConvolvedFixedFWHM::MakeSetPropertyFunctions( )
{
	AddSetScalarFunction( "setfwhm",
		[&, this](double value )
		{
			m_convolvedoptprop->SetFWHM( value );
			return true;
		}
	);

	AddSetObjectFunction( "setbasecrosssection",
		[&, this]( nxUnknown* obj )
		{
			return SetHighResProperties( obj );
		}
	);

}

/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_ConvolvedFixedFWHM::SetHighResProperties		 2015- 1- 9*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKOpticalProperty_Stub_ConvolvedFixedFWHM::SetHighResProperties( nxUnknown* userobject)

{
	bool										ok;
	skOpticalProperties*						useroptprop; 
	skWavelengthToPSF_Table*					userpsfinfo;
	skOpticalProperty_AdditionalStateInfo*		useratmosinfo;	

	useroptprop   = dynamic_cast<skOpticalProperties*>     ( userobject );					// By default many cross-section objects
	userpsfinfo   = dynamic_cast<skWavelengthToPSF_Table*> ( userobject );					// will support all three of these 
	useratmosinfo = dynamic_cast<skOpticalProperty_AdditionalStateInfo*> (userobject );	// interfaces through multiple inheritance

	ok = (useroptprop != nullptr);
	if (!ok)
	{
		nxLog::Record( NXLOG_ERROR, "ISKOpticalProperty_Stub_ConvolvedFixedFWHM::ExtractDetailPointersFromObject, The object passed in is not derived from skOpticalProperties. Are you using RawObjectPointer to get access to the underlying object");
	}
	else
	{

		if (userpsfinfo   == nullptr ) userpsfinfo   = &m_default_measurementdetails; 
		if (useratmosinfo == nullptr ) useratmosinfo = &m_default_atmosphericstateinfo;  
		ok =       m_convolvedoptprop->SetHighResOpticalProperties( useroptprop, userpsfinfo, useratmosinfo );
		ok = ok && m_convolvedoptprop->SetAtmosphericState( &m_defaultclimatology );
		if (!ok)
		{
			nxLog::Record(NXLOG_ERROR,"ISKOpticalProperty_Stub_ConvolvedFixedFWHM::SetHighResProperties, there were errors setting the high resolution cross-section properties");
		}
	}
	return ok;
}

