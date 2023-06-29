#include "../sktran_common.h"
//#include "modules/sktran_engine/include/indexofrefraction.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
using namespace boost::numeric::ublas;

/*-----------------------------------------------------------------------------
 *			skRTRefractiveIndex_Profile::~skRTRefractiveIndex_Profile		2012-8-27*/
/** 
 *	Code from http://www.crystalclearsoftware.com/cgi-bin/boost_wiki/wiki.pl?LU_Matrix_Inversion
 *	Following code inverts the matrix input using LU-decomposition with backsubstitution of unit
 *	vectors.  Reference: "Numerical Recipes in C, 2nd ed.", by Press, Teukolsky, Vetterling, & Flannery
 **/
/*---------------------------------------------------------------------------*/

static bool LocalFunction_MatrixInverse( matrix<double>& input, matrix<double>& inverse )
{
	using namespace boost::numeric::ublas;
	
	typedef permutation_matrix<std::size_t> pmatrix;

	// create a working copy of the input:
	matrix<double>		A(input);

	// create a permutation matrix for the LU-factorization:
	pmatrix pm(A.size1());

	// perform LU-factorization:
	int res = (int)lu_factorize(A,pm);
	if ( res != 0 ) return false;

	// create identity matrix of "inverse":
	inverse.assign( identity_matrix<double>(A.size1()) );

	// backsubstitute to get the inverse:
	lu_substitute( A, pm, inverse );

	return true;
}
/*-----------------------------------------------------------------------------
 *			skRTRefractiveIndex_Profile::skRTRefractiveIndex_Profile		2012-8-27*/
/** **/
/*---------------------------------------------------------------------------*/

skRTRefractiveIndex_Profile::skRTRefractiveIndex_Profile()
{
	m_calculator = new skRTRefractiveIndex_MoistAir;
}


/*-----------------------------------------------------------------------------
 *			skRTRefractiveIndex_Profile::~skRTRefractiveIndex_Profile		2012-8-27*/
/** **/
/*---------------------------------------------------------------------------*/

skRTRefractiveIndex_Profile::~skRTRefractiveIndex_Profile()
{
	if (m_calculator != NULL)	delete m_calculator;
}


/*-----------------------------------------------------------------------------
 *			skRTRefractiveIndex_Profile::CalculateProfile		2012-8-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTRefractiveIndex_Profile::CalculateProfile(SKTRAN_AtmosphericOpticalState_V21 *opticalstate, const SKTRAN_GridDefRayTracingShells_V21 *raytracingspecs, double wavelen_nm, GEODETIC_INSTANT referencepoint)
{
	bool					ok;
	//double					*shellheight;
	skClimatology			*species_temperature_pressure;
	size_t					heightidx;
	size_t					numshells;
	double					wavenum = ( 1 / wavelen_nm ) * 1.0E07;
	size_t					numbad;
	

	// (re-)calculate the index of refraction profile based on the opticalstate and ray tracing grid

	ok =	( referencepoint.heightm >= 0.0 ) && 
			( referencepoint.latitude >= -90.0 && referencepoint.latitude <= 90.0 ) &&
			( referencepoint.longitude >= -180.0 && referencepoint.longitude <= 360.0 ) &&
			( referencepoint.mjd >= 10000.0 );	// mjd is invalid for less than 10000 for some reason

	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING, "skRTRefractiveIndex_Profile::CalculateProfile, Reference point given not valid.");
	}

	if (ok)
	{
		m_referencepoint = referencepoint;

		ok			= ok && opticalstate->SetWavelength( wavelen_nm );
		ok			= ok && opticalstate->SetTimeAndLocation( referencepoint, true );

		numshells	= raytracingspecs->NumShells();

		if ( numshells == 0 )
		{
			nxLog::Record( NXLOG_WARNING, "skRTRefractiveIndex_Profile::CalculateProfile, Please set the shell altitudes for the refractivity grid.");
		}
		if ( ok &&  numshells != 0 )
		{
			m_heights.resize( numshells );
			m_temperature.resize( numshells );
			m_temperature.resize( numshells );
			m_refractiveindex.resize( numshells );

			//shellheight = raytracingspecs->begin();
			heightidx	= 0;
			
			// may need to add in check that h <= 120 km, since pressure/temperature profiles only valid till 120 km
			/*while ( shellheight != raytracingspecs->end() )
			{
				m_heights.at(heightidx) = *shellheight;	// midpoint of cell
				++heightidx;
				++shellheight;
			}*/

			for(size_t heightidx=0; heightidx< raytracingspecs->NumGridPoints(); heightidx++)
				m_heights.at(heightidx) = raytracingspecs->At(heightidx);
			
			opticalstate->GetAtmosphericStateModel(&species_temperature_pressure);

			if (ok)
			{
				m_pressure.resize(m_heights.size());
				m_temperature.resize( m_heights.size());

				species_temperature_pressure->GetHeightProfile( SKCLIMATOLOGY_PRESSURE_PA, m_referencepoint, &m_heights.front(), (int)m_heights.size(), &m_pressure.front(), true, &numbad );
				species_temperature_pressure->GetHeightProfile( SKCLIMATOLOGY_TEMPERATURE_K, m_referencepoint, &m_heights.front(), (int)m_heights.size(), &m_temperature.front(), true, &numbad );

				UpdateRefractiveIndex( opticalstate, wavenum );	
			}
		}
	}

	InitializeCubicSplineInterpolation( );	// once the index of refraction profile is found

	return ok;
}

double skRTRefractiveIndex_Profile::ExponentialLinearInterp(double altitude)
{
	auto upper = std::upper_bound(std::begin(m_heights), std::end(m_heights), altitude);

	if (upper == std::begin(m_heights))
	{
		// below the lowest altitude just return the first element
		return m_refractiveindex[0];
	}
	else if(upper == std::end(m_heights))
	{
		// above the highest altitude, just return the last element
		return m_refractiveindex[m_refractiveindex.size() - 1];
	}
	else
	{
		int upper_index = (int)std::distance(std::begin(m_heights), upper);
		int lower_index = upper_index - 1;

		double dh = m_heights[upper_index] - m_heights[lower_index];

		double w1 = (m_heights[upper_index] - altitude) / dh;
		double w2 = (1 - w1);

		return exp(w1 * log(m_refractiveindex[lower_index]) + w2 * log(m_refractiveindex[upper_index]));
	}
}



/*-----------------------------------------------------------------------------
 *			skRTRefractiveIndex_Profile::UpdateRefractiveIndex		2012-8-27*/
/** 
 *	This function is only called when it is verified that reference point, pressure(total) 
 *	and temperature are valid, and are saved as members to this class, so it should be safe.
 *
 *	If there is H20 present in the atmospheric state being used, then it will add
 *	effects on the refractivity due to water vapour.  Otherwise, it just calculates the index
 *	of refraction profile based on the wavenumber, temperature, and pressure.
 **/
/*---------------------------------------------------------------------------*/

void skRTRefractiveIndex_Profile::UpdateRefractiveIndex(SKTRAN_AtmosphericOpticalState_V21 *opticalstate, double wavenum )
{
	bool					ok;
	skClimatology			*species;	// search for other species, for now, only other one is water vapour
	std::vector<double>		wvp;
	size_t					numbad;
	size_t					numshells;
	double					Rv = 461.495; // specific gas constant for water vapour (J/(kg K))
	std::complex<double>	n;

	numshells = m_heights.size();

	const CLIMATOLOGY_HANDLE waterName = SKCLIMATOLOGY_H2O_CM3;
	ok = opticalstate->GetSpeciesClimatology( waterName, &species );	// get humidity profile, if it's there

	if (ok)
	{
		wvp.resize( m_heights.size() );
		species->GetHeightProfile( SKCLIMATOLOGY_H2O_CM3, m_referencepoint, &m_heights.front(), (int)m_heights.size(), &wvp.front(), true, &numbad );

		for (size_t i=0; i<numshells; i++)
		{
			// wvp is not yet correct, it is still a density:
			wvp[i] *= 1000000; // in per cm^3, change to per m^3
			wvp[i] *= Rv * m_temperature[i];
			wvp[i] *= 2.99151e-26; //kg/molecule
			// wvp is now correct for this height
			m_calculator->Set_WaterVapourPressure( wvp[i] );
			m_calculator->Set_TotalPressure( m_pressure[i] );
			m_calculator->Set_Temperature( m_temperature[i] );
			
			n = m_calculator->RefractiveIndex( wavenum );
			m_refractiveindex[i] = n.real();
		}
	}
	else
	{
		for (size_t i=0; i<numshells; i++)
		{
			m_calculator->Set_TotalPressure( m_pressure[i] );
			m_calculator->Set_Temperature( m_temperature[i] );
			
			n = m_calculator->RefractiveIndex( wavenum );
			m_refractiveindex[i] = n.real();
		}
	}

}


/*-----------------------------------------------------------------------------
 *			skRTRefractiveIndex_Profile::InitializeCubicSplineInterpolation		2012-7-6*/
/** 
 *	Function to find the vector M necessary to find the coefficients ai, bi, 
 *	ci, di for the equation yk = ai*(xk-xi)^3 + bi*(xk-xi)^2 + ci*(xk-xi) + di
 *  on the interval x_i <= x_k < x_{i+1}.  M is calculated for a quadratic 
 *	dropoff at the end points, so that M(1) = M(end) = 0.
 *
 *	This is intended for interpolating n(r) as a function of r.  It may be best
 *	to interpolate ln( n ) vs r, as n(r) is approximately exponential.
 *
 *	'num' is the size of y (i.e. n(r)), and h is the (constant) spacing of x,
 *	i.e. x_{i+1} - x_i ( in this case, r[i+1] - r[i] )
**/
/*---------------------------------------------------------------------------*/

bool skRTRefractiveIndex_Profile::InitializeCubicSplineInterpolation()
{
	using namespace boost::numeric::ublas;	// work in the boost 1.46.1 linear algebra namespace

	double						h;
	size_t						num;
	bool						ok;

	num = m_heights.size();
	h	= m_heights[2] - m_heights[1];	// assuming heights are evenly spaced

	matrix<double>				X(num-2,num-2);
	matrix<double>				Xinv(num-2,num-2);
	matrix<double>				Y(num-2,1);	// actually a vector
	matrix<double>				M(num-2,1);

	for ( size_t i=0; i < num-2; i++ )
	{
		for ( size_t j=0; j < num-2; j++ )
		{
			if ( i == j )
				X(i,j) = 4;
			else if (i == j+1 || j == i+1)
				X(i,j) = 1;
			else
				X(i,j) = 0;
		}

		Y(i,0) = m_refractiveindex[i] - 2*m_refractiveindex[i+1] + m_refractiveindex[i+2];

	}

	X(1,1) = 5;
	X(num-3,num-3) = 5;

	Y = ( 6/(h*h) )*Y;

	ok = LocalFunction_MatrixInverse( X, Xinv );	// take the inverse

	if (ok)
	{
		// Compute M:
		M = prod( Xinv, Y );

		// Now place those elements into member m_M:
		m_M.resize( num );

		for ( size_t k=0; k < num-2; k++ )
			m_M[k+1] = M(k,0);	
	
		m_M[0]		= m_M[1];		// parabolic run-out spline
		m_M[num-1]	= m_M[num-2];	// parabolic run-out spline
	}
	else
	{
		printf("skRTRefractiveIndex_Profile::InitializeCubicSplineInterpolation, Could not compute inverse of matrix.\n" );
	}

	return ok;
}


