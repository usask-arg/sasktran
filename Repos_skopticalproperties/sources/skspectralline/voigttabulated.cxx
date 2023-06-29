#include <skopticalproperties21.h>
#include <limits>
#include <boost/thread.hpp>

using namespace nxcgs;
/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtTabulated		2013-3-15*/
/** \ingroup spectrallineinternals
 *	 This is a storage buffer used for each spectral line that implements
 *	skSpectralLineShape_VoigtTabulated. Basically it just saves the doppler
 *	and Lorenz **/
/*---------------------------------------------------------------------------*/


static  const double  sqrtln2	 = 0.8325546111576977563531646448952;				// SQRT( LN(2) )
static  const double  twoln2     = 1.3862943611198906188344642429164;				// 2*log(2);
static  const double  oneAMU     = 1.66053886E-24;									// One atomic mass unit in grams 

skSpectralLineShape_VoigtTabulatedLUT skSpectralLineShapeStorageBuffer_VoigtTabulated::m_lookuptable;

/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtTabulated::skSpectralLineShapeStorageBuffer_VoigtTabulated		2013-3-15*/
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLineShapeStorageBuffer_VoigtTabulated::skSpectralLineShapeStorageBuffer_VoigtTabulated()
{
	ResetLineParams();
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtTabulated::~skSpectralLineShapeStorageBuffer_VoigtTabulated		2014-2-13*/
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLineShapeStorageBuffer_VoigtTabulated::~skSpectralLineShapeStorageBuffer_VoigtTabulated()
{
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtTabulated::ResetLineParams		2013-3-15*/
/** **/
/*---------------------------------------------------------------------------*/

void skSpectralLineShapeStorageBuffer_VoigtTabulated::ResetLineParams()
{
	m_nu00 = 0.0;
	m_aD   = 0.0;			// The doppler halh width
	m_aL   = 0.0;			// The Lorenz half width.
	m_Y    = std::numeric_limits<double>::quiet_NaN();
}

/*---------------------------------------------------------------------------
 *'					VoigtProfile::SetLineParams                      2002-9-17
 *
 *	nu00              = lines central wavelength in wavenumbers
 *	pressure          = air atmospheric pressure in dynes, cgs 
 *	partialpressure   = partial pressure of the species in dynes 
 *	temperature       = temperature in Kelvin
 *	mass              = molecular mass in grams
 *	LorentzWidthMult  = Lorentz Half Width at 1 atmosphere.
 *
 *	Calculates the Half width Half max of the Doppler (m_aD) and the 
 *	Half width half max of the Lorentz width (m_aL)
 *
 *	Doppler HWHM from Page 92, Exploration of the Solar System by Infra-red
 *	remote sensing.  Note their formula gives full width Half Max and we
 *`	want Half width half Max.
 *
 *	Lorentz pressure broadening from "The Hitran Molecular Spectroscopic Database
 *	and HAWKS (Hitran Atmospheric Workstation): 1996 Edition, Appendix equation A12
 *
 *-------------------------------------------------------------------------*/

bool skSpectralLineShapeStorageBuffer_VoigtTabulated::SetLineParams( double nu00,
																  double pressure,
																  double partialpressure,
																  double temperature,
																  double tref,
																  double tempcoeff,
																  double mass,
																  double airhalfwidth,
																  double selfhalfwidth)
{
	const double  pref = 1013250.0;										// Reference pressure in Dynes

	m_nu00 = nu00;
	m_aD   = nu00*sqrt(twoln2*KBOLTZMAN*temperature/(mass))/CLIGHT;
	m_aL   = pow(tref/temperature,tempcoeff)*( (airhalfwidth*(pressure-partialpressure)/pref)+(selfhalfwidth*partialpressure/pref) );
	m_Y    = sqrtln2*m_aL/m_aD;
	return true;
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_VoigtTabulated::skSpectralLineShape_VoigtTabulated		2013-3-14*/
/** The overriding purpose of this class is to convert a tab delimited (.txt) lookup table used by 
	the Voigt profile approximation engine into the binary format expected by the class HitranVoigtKBinary.
	This is original documentation that needs editing to make it accurately reflect the code but it is still useful.

	EXPLANATION OF VOIGT BROADENING APPROXIMATION
	---------------------------------------------

	The Voigt distribution V(deltaNu) for a distribution centred about wavenumber nu0 at a wavenumber nu = nu0 + deltaNu can be approximated as:

		V(deltaNu) = sqrt( ln(2)/pi ) * 1/gammaD * K(x,y)

	where gammaD is the Doppler distribution half-width (a parameter that depends on both the wavenumber in question & the temperature)
	and K(X,Y) is given by

		K = Re[w(z)]

	In turn, w(z) is given by

		w(z) = e^(-z^2) * ( 1 - erf(-i*z) )

	where, for our approximation,

		z = x + iy
		x = deltaNu/gammaD * sqrt(ln(2))
		y =  gammaL/gammaD * sqrt(ln(2))

	where gammaL is the Lorentz distribution half-width (a isotope-specific parameter given by the HITRAN database and modified by atmospheric conditions).

	Thus, once we have determined our current Doppler and Lorentz halfwidths, we can find values for x and y;
	these values are then used to calculate a K(x,y), which allows us to calculate V(deltaNu), the Voigt distribution function at the given point.


	THE TABLE
	---------
	The lookup table contains precalculated values of K for the following range of x,y:

	0 <= x <= 3.9
	0 <= y <= 3.0

	Outside of the bounds of the table (0<x<3.9, 0<y<3), K(x,y) can be computed by analytic functions given in Abramowitz and Stegun (see cite note below).
	These analytic functions are hard-coded into the sister class HitranVoigtKBinary.

	The table found to work accurately with this code used a resolution of 0.001 in both the x and y ranges; a table of lower resolution can be found in:

	Abramowitz & Stegun 
	Handbook of Mathematical Functions, 9th Printing
	Table 7.9: Error Function for Complex arguments, pp. 325-328.

	Note that the table in Abramowitz & Stegun contains both the real and imaginary parts of the solution,
	while the VoigtK lookup table only contains the real part, as the imaginary part is unneeded for our approximation.


	THE TEXT FILE
	-------------

	An extremely low resolution version of the table format within the text value can be found below; this table provides
	K-values for x = 0-3, y = 0-2:

	VoigtK[X,Y]	0	1	2
	0	1	0.427584	0.255396
	1	0.367879	0.304744	0.218493
	2	0.0183156	0.14024	0.147953
	3	0.00012341	0.0653178	0.0927108

	The first entry on the first line is a title to remind the user of what data is stored in the table.  
	Following that are the column labels, tab-delimited for easy computer reading.  On the ensuing lines are the row labels,
	followed by the actual tabular data corresponding to each row and column.  Again, these values are tab-delimited.

	Note that the table is indexed similar to standard matrix notition; X, the first variable, denotes the row, 
	while Y, the second variable, denotes the column.

	The table is assumed to have equally spaced indicies in both dimensions; however, these dimensions need not be the same for x and y.
	This assumption comes into play in the sister class HitranVoigtKBinary, which reads the first, second, and final labels for the rows and columns.
	The resolution is determined by the difference between entries 1 and 2, while the class calculates the total number of entries by 
	using this resolution along with the first and final labels.


	GENERATING THE TEXT FILE USING A MATLAB SCRIPT
	----------------------------------------------
	% declaring and allocating
	X = (0:1:3)';
	Y = (0:1:3)';
	Z = zeros( length(X), length(Y) );
	K = zeros( length(X), length(Y) );

	%a column-by-column calculation to hopefully increase efficiency
	for j = 1:length(Y)
		Z(:,j) = X + 1i*Y(j);
	    
		S = zeros(length(X), 1); %our summation approximation to erf()
		for n = 0:150
			S = S + (-1)^n * (-1i*Z(:,j)).^(2*n+1) /factorial(n)/(2*n+1);
		end
	    
		%grabbing the real part
		K(:,j) = real( exp( Y(j)^2 - X.^2 ) .* ( cos(2*X.*Y(j)) - 1i.*sin(2*X.*Y(j)) ) .* (1 - 2*S/sqrt(pi) ) ); 
	end

	%outputting to a text file
	fid = fopen( 'XYlookupK01.txt', 'w');
	fprintf(fid, '%s', 'VoigtK[X,Y]');

	for j = 1:length(Y)
		fprintf( fid, '\t%g', Y(j) );
	end
	for i = 1:length(X)
		fprintf( fid, '\n%g', X(i) );
		for j = 1:length(Y)
			fprintf( fid, '\t%g', K(i,j) );
		end
	end

	fclose(fid);


	CLASS OPERATIONS
	----------------

	After opening the lookup table text file, HitranVoigtKText acquires the X and Y parameters (X and Y are assumed to start at zero, and so we only store
	the resolution and the number of records in the table).  

	Now that the size of the table is known, HitranVoigtKText is able to allocate sufficient memory.  It then loads the entire table into memory, 
	into a single 2D array.  This array, along with the table parameters, is then written to the output binary file.  Finally, the table memory is freed.

	Thus, the output binary file leads with 2 ints and 2 doubles (numXelements, numYelements, xResolution, yResolution), followed by the table data 
	in an uniterrupted stream of doubles.

	A basic outline of use:
	1) Class is initialized with null values for the input/output file paths;
	2) I/O file pathways are set using the public member function SetIOFilePaths(char voigtKTextPathIn[], char voigtKBinPathIn[]);
	3) The conversion is performed using the public member function TextToBinary(), creating an output binary file at the specified point.

	EXAMPLE CODE:
	{
		HitranVoigtKText		voigtKText;

		voigtKText.SetIOFilePaths( "C:\\Hitran\\XYlookupK.txt", "C:\\HitranBin\\XYlookupK.bin" );
		voigtKText.TextToBinary();
	}
	(the above assumes that <XYlookupK.txt> is stored at C:\Hitran, and that C:\HitranBin is a valid location to write to.)

**/
/*---------------------------------------------------------------------------*/

skSpectralLineShape_VoigtTabulatedLUT::skSpectralLineShape_VoigtTabulatedLUT()
{
	m_xresolution = 0.001;
	m_yresolution = 0.001;
	m_numx  = (size_t)(3.9/m_xresolution + 0.5)+1;
	m_numy  = (size_t)(3.0/m_yresolution + 0.5)+1;
	m_table.resize(m_numy);
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_VoigtTabulated::Yindex		2013-3-15*/
/** **/
/*---------------------------------------------------------------------------*/

size_t skSpectralLineShape_VoigtTabulatedLUT::Yindex( double Y ) const
{
	NXASSERT(Y >=0.0);
	return (size_t)(Y/m_yresolution + 0.5);				// Do nearest neighbour interpolation
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_VoigtTabulated::Xindex		2013-3-15*/
/** **/
/*---------------------------------------------------------------------------*/

size_t skSpectralLineShape_VoigtTabulatedLUT::Xindex( double X ) const
{
	NXASSERT(X >=0.0);
	return (size_t)(X/m_xresolution + 0.5);				// Do nearest neighbour interpolation
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_VoigtTabulated::CalculateEntries		2013-3-14*/
/** 
/*	This calculates the voigt profile for a whole slew of X,Y values between 0 and 3 (or 3.9).
 *	I have spot checked the values against equivalent matlab  code and it was giving the
 *	same answer so I'm prettu confident its good.
 *	The equivalent matlab code is,

	X = (0:0.001:3.9)';
	Y = (0:0.001:3)';
	Z = zeros( length(X), length(Y) );
	K = zeros( length(X), length(Y) );

	%a column-by-column calculation to hopefully increase efficiency
	for j = 1:length(Y)
		Z(:,j) = X + 1i*Y(j);

		S = zeros(length(X), 1); %our summation approximation to erf()
		for n = 0:150
			S = S + (-1)^n * (-1i*Z(:,j)).^(2*n+1) /factorial(n)/(2*n+1);
		end

		%grabbing the real part
		K(:,j) = real( exp( Y(j)^2 - X.^2 ) .* ( cos(2*X.*Y(j)) - 1i.*sin(2*X.*Y(j)) ) .* (1 - 2*S/sqrt(pi) ) ); 
	end

**/
/*---------------------------------------------------------------------------*/

bool skSpectralLineShape_VoigtTabulatedLUT::CalculateXEntries( std::vector<double>& K, size_t yindex ) const
{
	double					Y;
	double					X;
	size_t					ix;
	size_t					n;
	double					minus1n;
	double					twosqrtpi;
	std::vector<double>		factorialn;

	Y = yindex*m_yresolution;
	K.resize(m_numx);
	factorialn.resize(151);
	factorialn.at(0) = 1;
	for (n=1; n < 151; n++)
	{
		factorialn.at(n) = factorialn.at(n-1)*n;
	}
	minus1n = 1;
	for (n=0; n < 151; n++)
	{
		factorialn.at(n) = 1.0/( minus1n*factorialn.at(n)*(2*n+1) );
		minus1n         *= -1;
	}
	twosqrtpi = 2.0/sqrt(nxmath::Pi);
	for (ix=0; ix < m_numx; ix++ )
	{
		X = ix*m_xresolution;
		std::complex<double>	ZSTAR(Y,-X);
		std::complex<double>	S    (0,0);
		std::complex<double>	Q;

		for (n=0; n < 151; n++)
		{
			S += pow(ZSTAR, (int)(2*n+1) )*factorialn.at(n);
		}
		std::complex<double>	R( cos(2.0*X*Y), -sin(2.0*X*Y));
	
		Q        = R * (1.0 - twosqrtpi*S); 
		K.at(ix) = Q.real()*exp( Y*Y - X*X );
	}
	return true;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_VoigtTabulatedLUT::LoadDirectoryNameFromRegistry		2013-3-15*/
/** **/
/*---------------------------------------------------------------------------*/

nxString skSpectralLineShape_VoigtTabulatedLUT::LoadDirectoryNameFromRegistry( )
{
	nxRegistryConfiguration		config( "USask-ARG", "skOpticalProperties/VoigtTables/",nxRegistryConfiguration::GLOBAL_INI, true);
	nxString					filename;
	bool						ok;

	ok = config.LocateDirectoryFromKey( "CacheDirectory", &filename, true, true, "Enter location for caching voigt table parameters");
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skSpectralLineShape_VoigtTabulated::LoadDirectoryNameFromRegistry, error loading aersol cache directory name from registry, using TEMP");
		filename = getenv("TEMP");
	}
	return filename;
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_VoigtTabulatedLUT::FullCacheName		2013-3-15*/
/** **/
/*---------------------------------------------------------------------------*/

nxString skSpectralLineShape_VoigtTabulatedLUT::FullCacheName( size_t yindex)
{
	nxString	dirname;
	nxString	basename;
	nxString	fullfilename;

	dirname = LoadDirectoryNameFromRegistry();
	dirname.EnsureLastCharIsDirectoryChar();
	basename.sprintf("voigttables/voigtXY%06Iu.dat",(size_t)yindex);	
	fullfilename = dirname + basename;
	fullfilename.MakeDirectorySeparatorsOSConsistent();
	return fullfilename;
}


/*-----------------------------------------------------------------------------
 *				skSpectralLineShape_VoigtTabulatedLUT::ReadCacheFile		2009-6-4*/
/*---------------------------------------------------------------------------*/

bool skSpectralLineShape_VoigtTabulatedLUT::ReadCacheFile( std::vector<double>&	K, const char* filename )
{
	bool			ok;
	nxFile			f;
	double*			rowptr;
	unsigned int	numx;

	ok = nxDirectory::FileExists(filename);					// See if the cache exists
	if (ok)
	{
		K.resize(m_numx);
		rowptr = &K.front();

		f.Open( filename, "rb");
		ok = f.IsOpen();
		if (ok)
		{
			ok =       (f.Read( &numx,	sizeof(numx),      1 ) == 1);
			ok = ok && (numx == m_numx);
			ok = ok && (f.Read( rowptr,	sizeof(double), m_numx) == m_numx);
			f.Close();
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"skSpectralLineShape_VoigtTabulated::ReadCacheFile, Error reading Voigt table cache cache <%s>. Thats a problem!",(const char*)filename);
			}
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *				skSpectralLineShape_VoigtTabulatedLUT::ReadCacheFile		2009-6-4*/
/*---------------------------------------------------------------------------*/

bool skSpectralLineShape_VoigtTabulatedLUT::WriteCacheFile( const std::vector<double>&	K, const char* userfilename )
{
	bool			ok;
	nxFile			f;
	const double*	rowptr;
	unsigned int	numx;
	nxString		filename( userfilename );
	filename.MakeDirectorySeparatorsOSConsistent();
	nxFileSpec		spec(filename);


	ok = nxDirectory::CreateADirectory( spec.FullDirSpec() );
	ok = ok && (K.size() == m_numx);
	if (ok)
	{
		numx   = (unsigned int)m_numx;
		rowptr = &K.front();
		f.Open( filename, "wb");
		ok = f.IsOpen();
		if (ok)
		{
			ok =       (f.Write( &numx,	 sizeof(numx),      1 ) == 1);
			ok = ok && (f.Write( rowptr, sizeof(double), m_numx) == m_numx);
			f.Close();
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skSpectralLineShape_VoigtTabulated::WriteCacheFile, Error writing Voigt table cache <%s>. Thats a problem!",(const char*)filename);
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_VoigtTabulatedLUT::CheckTableEntry		2013-3-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineShape_VoigtTabulatedLUT::CheckTableEntry(size_t yindex)
{
	bool	ok;
	static std::mutex		cachemutex;

	ok = (m_table.at(yindex).size() == m_numx);
	if (!ok)
	{
		nxString	fullname;

		std::lock_guard<std::mutex>		lock( cachemutex);				// Lock the cache mutex so no other thread is reading/or writing at this time
		fullname = FullCacheName( yindex);
		ok = ReadCacheFile( m_table.at(yindex), fullname );
		if (!ok)															// If !ok then we need to create it. Lets watch out for thread sync.
		{																	// only allow 1 thread at a time through this section
			ok = ReadCacheFile( m_table.at(yindex), fullname );				// Try and read it againas a previous thread might have just made it
			if (!ok)														// But if it failed then we need to create it so do it.
			{
				ok =       CalculateXEntries( m_table.at(yindex), yindex);
				ok = ok && WriteCacheFile( m_table.at(yindex), fullname);
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING,"skSpectralLineShape_VoigtTabulated::CheckTableEntry, There were errors creating the voigt table cache thats not good.");
				}
			}
		}
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_VoigtTabulated::skSpectralLineShape_VoigtTabulated		2014-2-12*/
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLineShape_VoigtTabulated::skSpectralLineShape_VoigtTabulated()
{
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_VoigtTabulated::NormalizedShape		2013-3-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineShape_VoigtTabulated::NormalizedShape( double X, double Y, double* W, skSpectralLineShape_VoigtTabulatedLUT* lut) const
{
	static const double		g_voigtCoeff[10]	= {0.5124242,0.2752551,0.05176536,2.724745,0.4613135,0.1901635,0.09999216,1.7844927,0.002883894,5.5253437};
	size_t					yindex;
	size_t					xindex;
	bool					outofrange;
	bool					ok = true;


	yindex     = lut->Yindex(Y);
	xindex     = lut->Xindex(X);
	outofrange = (yindex >= lut->NumY()  || (xindex >= lut->NumX()) );

	if (outofrange)
	{
		std::complex<double>	i(0,1);
		std::complex<double>	Z(X,Y);
		std::complex<double>	W_Z;
		std::complex<double>	Z2 = Z*Z;

		if ( X > 6.0 || Y > 6.0 )
		{
			W_Z = i * Z * ( g_voigtCoeff[0]/( Z2 - g_voigtCoeff[1] ) + g_voigtCoeff[2]/( Z2 - g_voigtCoeff[3]) );
		}
		else if ( X > 3.9 || Y > 3.0 )
		{ 
			W_Z = i * Z * ( g_voigtCoeff[4]/( Z2 - g_voigtCoeff[5] ) + g_voigtCoeff[6]/( Z2 - g_voigtCoeff[7]) + g_voigtCoeff[8]/( Z2 - g_voigtCoeff[9]) );
		}
		*W = W_Z.real();
	}
	else 
	{
		ok      = lut->CheckTableEntry(yindex);														// Make sure table is loaded
		*W      = ok ? lut->At(yindex,xindex) : std::numeric_limits<double>::quiet_NaN();	// And get the nearest neighbour value if ok, else return NaN.
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					ExtractStorageBuffer		2013-3-15*/
/** **/
/*---------------------------------------------------------------------------*/

static skSpectralLineShapeStorageBuffer_VoigtTabulated*	ExtractStorageBuffer	( skSpectralLineShapeStorageBuffer* storagebuffer)
{
	return skSpectralLineShapeStorageBuffer_Type< skSpectralLineShapeStorageBuffer_VoigtTabulated >::Cache( storagebuffer);
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_VoigtTabulated::CreateStorageBuffer		2013-3-7*/
/** Create an instance of the storage buffer. It looks big and ugly
 *	but its actually quite simple
 **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineShape_VoigtTabulated::CreateStorageBuffer( skSpectralLineShapeStorageBuffer**	storagebuffer )
{
	skSpectralLineShapeStorageBuffer_Type< skSpectralLineShapeStorageBuffer_VoigtTabulated >*	buffer;
	bool																						ok;

	buffer = new skSpectralLineShapeStorageBuffer_Type< skSpectralLineShapeStorageBuffer_VoigtTabulated >;		// Allocate the storage object
	ok = (buffer != NULL);																						// 
	if ( ok)
	{
		buffer->AddRef();
	}
	*storagebuffer = dynamic_cast<skSpectralLineShapeStorageBuffer*>( buffer );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_VoigtTabulated::ConfigureLineParameters		2013-3-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineShape_VoigtTabulated::ConfigureLineParameters( const skSpectralLine*						spectralline,
																  double								temperature,
																  double								pressure,
																  const GEODETIC_INSTANT&				geopt,
															      skClimatology*						/*atmosphericstate*/,
															      skSpectralLineShapeStorageBuffer*		storagebuffer )
{
	skSpectralLineShapeStorageBuffer_VoigtTabulated*	cache;
	bool			ok;
	double			tref;
	double			nu00;
	double			partialpressure;
	double			tempcoeff;
	double			mass;
	double			airhalfwidth;
	double			selfhalfwidth;
	const double	pref    = 1013250.0;									// Reference pressure in Dynes
	double			patm;

	cache = ExtractStorageBuffer( storagebuffer );					// set m_cache to point to storage buffer
	ok =       (cache != NULL);
	if (ok)
	{
		tref            = spectralline->Tref();				// reference temperature for temp coeff
		nu00            = spectralline->Nu();				// Get the wave number of the spectral line
		tempcoeff       = spectralline->Nair();				// Get the temperature coeffiecient
		partialpressure = spectralline->ParentMolecule()->PartialPressure( geopt, pressure, temperature );
		mass            = spectralline->ParentMolecule()->MassAMU();
		airhalfwidth    = spectralline->GammaAir();
		selfhalfwidth   = spectralline->GammaSelf();
		pressure        = pressure*10.0;			// Convert Pascals to Dynes
		partialpressure = partialpressure*10.0;		// Convert Pascals to Dynes.
		mass            = mass*oneAMU;				// Convert the AMU to grams
		patm            = pressure/pref;
		nu00           += spectralline->Deltaair()*patm;

		ok = cache->SetLineParams(  nu00, pressure, partialpressure, temperature, tref, tempcoeff, mass, airhalfwidth, selfhalfwidth );
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skSpectralLineShape_VoigtTabulated::ConfigureLineParameters, There was an error configuring the line parameters for the Voigt Humlik code");
		cache->ResetLineParams();
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_VoigtTabulated::LineShapeFunction		2013-3-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineShape_VoigtTabulated::LineShapeFunction(	double								nu,						// Calculate line shape function at this wavenumber
															double*								uservalue,				// returns the line shape value in this variable
															const skSpectralLine*				spectralline,			// Spectral line for line shape function
															skSpectralLineShapeStorageBuffer*	storagebuffer ) const	// Storage buffer from an earlier call to ConfigureStorageBuffer
{
	skSpectralLineShapeStorageBuffer_VoigtTabulated*	cache;
	double												deltaNu;
	double												X;
	double												Y;
	bool												ok;
	static const double									sqrtln2divpi = 0.46971863934982566688617016420509; // Square root ( log2)/Pi)

	cache      =       ExtractStorageBuffer( storagebuffer );					//get the cache associated with this spectral line
	ok         =       (cache != NULL);											// make sutre it is good
	ok         = ok && cache->IsConfigured();									// and make sure it is properly configured
	if (ok)																		// If it is
	{																			// then
		deltaNu	   = fabs( nu - cache->Nu00() );								// Get the info 
		X		   = deltaNu*sqrtln2/cache->aD();								// required to calculate the 
		Y          = cache->Y();												// Voigt profile
		ok         = NormalizedShape( X, Y, uservalue, cache->LUTVar());		// and do it.
		*uservalue *= spectralline->LineIntensity()*sqrtln2divpi/cache->aD();	// klj588, Add the proper normalization
	}
	if (!ok)
	{
		*uservalue = std::numeric_limits<double>::quiet_NaN();
	}
	return ok;
}

