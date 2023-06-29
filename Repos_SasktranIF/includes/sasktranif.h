#if !defined(SASKTRANIF_HEADER_INCLUDE)
#define SASKTRANIF_HEADER_INCLUDE
#pragma once
#include "nxbase_core.h"
#include "nxbase_math.h"
#include "nxbase_geodesy.h"
#include <vector>
#include <functional>
#include <complex>

#include "climatology_handles.h"
#include "sasktran_polarization.h"

bool operator < ( const CLIMATOLOGY_HANDLE& a,  const CLIMATOLOGY_HANDLE& b);	//!, The less than operator for CLIMATOLOGY_HANDLE is defined in our code , so we can use it as a key for std::map


/*---------------------------------------------------------------------------
 *                    Class ISKModuleBase_Stub                    2019-10-01 */
/** **/
/*---------------------------------------------------------------------------*/

class ISKModuleBase_Stub : public nxUnknown
{
	private:
		std::map< std::string, std::function<bool(double)> >			 m_scalarsetfunctions;
		std::map< std::string, std::function<bool(const char*)> >	     m_stringsetfunctions;
		std::map< std::string, std::function<bool(const double*, int)> > m_vectorsetfunctions;
		std::map< std::string, std::function<bool(nxUnknown*)> >		 m_objectsetfunctions;
		std::map< std::string, std::function<bool(double*)> >			 m_scalargetfunctions;
		std::map< std::string, std::function<bool(int)> >				 m_vectorgetfunctions;

	protected:
		std::vector<double>												 m_getpropertybuffer;

	private:
		bool					ParseCommandAndIndex		( const char * input, std::string& cmd, int& index );
		bool					GetPropertyScalar			( const char* propertyname, double* value );

	protected:
		bool					AddSetScalarFunction		( const char * name, std::function<bool(double)> func);
		bool					AddSetStringFunction		( const char * name, std::function<bool(const char*)> func);
		bool					AddSetVectorFunction		( const char*  name, std::function<bool(const double*, int)> func);
		bool					AddSetObjectFunction		( const char*  name, std::function<bool(nxUnknown*)> func);
		bool					AddGetScalarFunction		( const char*  name, std::function<bool(double*)> func);
		bool					AddGetVectorFunction		( const char*  name,  std::function<bool(int)> func);

	public:
								ISKModuleBase_Stub			();
		virtual				   ~ISKModuleBase_Stub			();
		virtual bool			SetPropertyScalar			( const char* propertyname, double value );
		virtual bool			SetPropertyArray			( const char* propertyname, const double* value, int numpoints );
		virtual bool			SetPropertyObject			( const char* propertyname, nxUnknown* object );
		virtual bool			SetPropertyString			( const char* propertyname, const char* str);
		virtual bool			GetProperty					( const char* propertyname, const double** value, int* numpoints );
};

/*-----------------------------------------------------------------------------
 *					ISkClimatology_Stub		2014-2-8*/
/** This is the base class implemented by each skClimatology in each
 *	DLL
 **/
/*---------------------------------------------------------------------------*/

class ISKClimatology_Stub : public ISKModuleBase_Stub
{
	public:
								ISKClimatology_Stub			(){}
		virtual 			   ~ISKClimatology_Stub			(){}
		virtual nxUnknown*		RawObjectPointer				() = 0;
		virtual bool			UpdateCache					( const GEODETIC_INSTANT& location )                                                    = 0;
		virtual bool			GetParameter				( const CLIMATOLOGY_HANDLE& species,  const GEODETIC_INSTANT& location, double* value ) = 0;
		virtual bool			SetPropertyUserDefined		( const CLIMATOLOGY_HANDLE& species,  double* profilevalues, int numpoints) = 0;
};


/*-----------------------------------------------------------------------------
 *					ISkOpticalProperty_Stub		2014-2-8*/
/** This is the base class implemented by each skOpticalProperty in each
 *	DLL
 **/
/*---------------------------------------------------------------------------*/

class ISKOpticalProperty_Stub : public ISKModuleBase_Stub
{
	public:
								ISKOpticalProperty_Stub		(){}
		virtual 			   ~ISKOpticalProperty_Stub		(){}
		virtual nxUnknown*		RawObjectPointer				() = 0;
		virtual bool			SetAtmosphericState			( ISKClimatology_Stub* atmosphere ) = 0;
		virtual bool			SetLocation					( const GEODETIC_INSTANT& pt ) = 0;
		virtual bool			InternalClimatology_UpdateCache	( const GEODETIC_INSTANT& pt)								= 0;
		virtual bool			CalculateCrossSections		( double wavenumber,		                    double* absxs,  double* extxs, double* scattxs) = 0;
		virtual bool			CalculateCrossSectionsArray	( const double * wavenumber, int numwavenumber, double *absxs,  double* extxs, double* scattxs) = 0;
		virtual bool			CalculatePhaseMatrix		( const double * wavenumber, const double * cosscatterangle, double* phasematrix ) = 0;
		virtual bool			AddUserDefined				( double temperature,  double* wavelen_nm, int numwave, double* crosssection, int numcross) = 0;
		virtual bool            AddUserDefinedPressure      ( double* pressure, int numpressure, double* temperature, int numtemperature, double* wavelen_nm, int numwavel, double* crosssection, int numcross, double broadnervmr ) = 0;
};


/*-----------------------------------------------------------------------------
 *					ISKEmission_Stub		2014-2-8*/
/** This is the base class implemented by each skEmission in each
 *	DLL
 **/
/*---------------------------------------------------------------------------*/

class ISKEmission_Stub : public ISKModuleBase_Stub
{
	public:
								ISKEmission_Stub				(){};
		virtual 			   ~ISKEmission_Stub				(){};
		virtual nxUnknown*		RawObjectPointer				() = 0;
//		virtual bool			SetAtmosphericState				( ISKClimatology_Stub* atmosphere)               = 0;		// This should go in through an object specific property
		virtual bool			UpdateLocation					( const GEODETIC_INSTANT& pt, bool isground )      = 0;
		virtual bool			UpdateCache						( const GEODETIC_INSTANT& pt )      = 0;
		virtual bool			IsotropicEmission				( double wavenumber, double* isotropicradiance) = 0;
		virtual bool			IsotropicEmissionArray			( const double * wavenumber, int numwavenumber, double * isotropicradiance, int numisotrop) = 0;
};

/*-----------------------------------------------------------------------------
 *					ISKBrdf_Stub		2014-2-8*/
/** This is the base class implemented by each SKBRDF in each
 *	DLL
 **/
/*---------------------------------------------------------------------------*/

class ISKBrdf_Stub : public ISKModuleBase_Stub
{
	public:
								ISKBrdf_Stub				(){};
		virtual 			   ~ISKBrdf_Stub				(){};
		virtual nxUnknown*		RawObjectPointer				() = 0;
		virtual bool			BRDF							( double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* brdf) = 0;
};

/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum_Stub		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKSolarSpectrum_Stub : public ISKModuleBase_Stub
{
	public:
								ISKSolarSpectrum_Stub				(){};
		virtual 			   ~ISKSolarSpectrum_Stub				(){};
		virtual nxUnknown*		RawObjectPointer					() = 0;
		virtual bool			Irradiance							( double        wavelen_nm_vacuum, double* irradiance ) = 0;
		virtual bool			IrradianceArray						( const double* wavelen_nm_vacuum, double* irradiance, int numpoints) = 0;
		virtual bool			IrradianceAt1AU						( double        wavelen_nm_vacuum, double* irradiance ) = 0;
		virtual bool			IrradianceAt1AUArray				( const double* wavelen_nm_vacuum, double* irradiance, int numpoints ) = 0;
		virtual bool			NanometerResolutionFWHM				( double        wavelen_nm_vacuum, double* resolution_nm_fwhm) = 0;
		virtual bool			NanometerResolutionFWHMArray		( const double* wavelen_nm_vacuum, double* resolution_nm_fwhm, int numpoints) = 0;
		virtual bool			SampleSpacing						( double wavelen_nm_vacuum, double* sample_spacing) = 0;
		virtual bool			SampleSpacingArray					( const double* wavelen_nm_vacuum, double* sample_spacing, int numpoints) = 0;
		virtual bool			SetSolarDistanceFromMjd				( double mjd ) = 0;
		virtual bool			SetSolarDistanceFromAU				( double au ) = 0;
		virtual bool			MinValidWavelength					( double* minwavelength_nm) = 0;
		virtual bool			MaxValidWavelength					( double* minwavelength_nm) = 0;

};


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub									2014-2-8*/
/** The ISKEngine stub class. Each engine in each DLL implements a class
 *	derived from this class. This allows each radiative transfer engine
 *	to be fully encapsulated within the interface. The ISKEngine class
 *	loads the stub class in its constructor. This object uses nxUnknown
 *	for its lifetime management: this lets us avoid using the
 *	C/C++ runtime library from one compiler to deallocate code allocated
 *	with another compiler which would happen with std::sharedptr.
 **/
/*---------------------------------------------------------------------------*/

class ISKEngine_Stub : public ISKModuleBase_Stub
{
	public:
								ISKEngine_Stub				(){};
		virtual 			   ~ISKEngine_Stub				(){};
		virtual nxUnknown*		RawObjectPointer			() {return nullptr;} // engines cannot be accessed
		virtual bool			AddLineOfSight				( double mjd,  const nxVector& observer, const nxVector& lookvector, int* losindex  ) = 0;
		virtual bool			AddSpecies					( const CLIMATOLOGY_HANDLE& species, ISKClimatology_Stub* climatology, ISKOpticalProperty_Stub* opticalproperty) = 0;
		virtual bool			AddEmission					( const EMISSION_HANDLE& species,    ISKEmission_Stub*    emission ) = 0;
		virtual bool			SetAtmosphericState			( ISKClimatology_Stub* climatology ) = 0;
		virtual bool			SetAlbedo					( double albedo ) = 0;
		virtual bool			SetBRDF						( ISKBrdf_Stub* brdf ) = 0;
		virtual bool			SetPolarizationMode			( int polarizationmode) = 0;
		virtual bool			SetWavelengths				( const double* wavelen, int numwavelen ) = 0;
		virtual bool			InitializeModel				() = 0;
		virtual bool			CalculateRadiance			( const double**			   radiance,  int* numwavelens, int* numlinesofsight) = 0;
		virtual bool			CalculateStokesVector		( const ISKStokesVector** radiancep, int* numwavelens, int* numlinesofsight) = 0;
		virtual bool			GetWeightingFunctions	    ( const double**          wf,        int* numwavel,    int* numlinesofsight, int* numwf ) = 0;
};


/*-----------------------------------------------------------------------------
 *					ISKGeodetic_Stub							2014- 5- 6*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKGeodetic_Stub : public ISKModuleBase_Stub
{

   	public:
							ISKGeodetic_Stub			( ){};
		virtual			   ~ISKGeodetic_Stub			( ){};
      	virtual bool		FromGeodetic				( double latitude, double longitude,  double Height ) = 0;				//!< Set the current location using the specified geodetic coordinates
      	virtual bool		FromGeocentric				( const nxVector& geocentric  ) = 0;									//!< Set the current location from the specified geocentric X,Y,Z vector. (all in meters).
		virtual bool		FromTangentPointLocation	( const nxVector& r, const nxVector& lookv ) = 0;						//!< Set the current location from the implied tangent point
		virtual bool		FromTangentAltitude			( double required_height, const nxVector& spacecraftlocation, const nxVector& boresightplane, nxVector* requiredlookvector) = 0;	//!< Set the current location from the tangent point at the specified height. Also return the \e look \e vector necessary to do this.
      	virtual nxVector	GeodeticWest				( ) = 0;			//!< Get the topocentric unit vectors at the current location
      	virtual nxVector	GeodeticSouth				( ) = 0;			//!< Get the topocentric unit vectors at the current location
      	virtual nxVector	GeodeticUp					( ) = 0;			//!< Get the topocentric unit vectors at the current location
		virtual nxVector	Location					( ) = 0;
		virtual double		GeodeticLongitude			( ) = 0;
		virtual	double		GeodeticLatitude			( ) = 0;
		virtual double		Height						( ) = 0;
		virtual bool		GetShellHeightLocation		( double H, const nxVector& observerposition, const nxVector& look, nxVector* entrypoint, nxVector* exitpoint) = 0;
		virtual nxVector	OsculatingSpheroidCenter    ( ) = 0;
		virtual double		OsculatingSpheroidRadius    ( ) = 0;

};

class ISKMie_Stub : public ISKModuleBase_Stub
{
	public:
		ISKMie_Stub() {};
		virtual ~ISKMie_Stub() {};

		// Function to set the inputs
		virtual bool Calculate(double lambda, double radius, double refrac_real, double refrac_imag) = 0;

		// Output functions
		virtual double								Qext() = 0;
		virtual double								Qsca() = 0;
		virtual double								Qabs() = 0;
		virtual double								Cext() = 0;
		virtual double								Csca() = 0;
		virtual double								Cabs() = 0;
		virtual nx1dArray< std::complex<double> >*  S1() = 0;
		virtual nx1dArray< std::complex<double> >*	S2() = 0;
		virtual nx2dArray<double>*					PMom() = 0;
		virtual std::complex<double>				SForward() = 0;
		virtual std::complex<double>				SBackward() = 0;
		virtual std::complex<double>				TForward(int i) = 0;
		virtual std::complex<double>				TBackward(int i) = 0;
		virtual double								Spike() = 0;
};

#endif
