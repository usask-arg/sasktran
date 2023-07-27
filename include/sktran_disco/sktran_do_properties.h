#pragma once
#include "sktran_disco/sktran_do.h"
#include "sktran_do_specs.h"

namespace sasktran_disco
{
	// An object which inherits the requested read-only properties and safely
	// copies them to this object. This object should be inherited by an 
	// object which would like access to read-only properties.
	template<class... Ts>
	class ReadOnlyProperties: public Ts...
	{
	protected:
		ReadOnlyProperties() = delete;

		// Only constructible from a super property.
		template<class SuperProperty>
		ReadOnlyProperties(const SuperProperty& super):
			Ts(static_cast<const Ts*>(&super))...
		{
		}
	};

#pragma region "Read only properties"
	// Requirements for a properties:
	// Public: Copy constructor from const pointer
	// Protected: Data members. Default constructor.

    template <int NSTOKES, int CNSTR=-1>
	class BasicProperties
	{
	public:
		BasicProperties(const BasicProperties* other):
			M_NSTR(other->M_NSTR),
			M_NLYR(other->M_NLYR),
			M_MU(other->M_MU),
			M_WT(other->M_WT),
			M_LP_MU(other->M_LP_MU),
			M_LAYER_CONSTRUCTION(other->M_LAYER_CONSTRUCTION),
			M_USE_PSEUDO_SPHERICAL(other->M_USE_PSEUDO_SPHERICAL),
			M_USE_LOS_SPHERICAL(other->M_USE_LOS_SPHERICAL),
            M_USE_GREENS_FUNCTION(other->M_USE_GREENS_FUNCTION),
			M_NUM_SZA(other->M_NUM_SZA),
			M_SS_ONLY(other->M_SS_ONLY),
            M_SZA_REL_SEP(other->M_SZA_REL_SEP)
		{
			// empty
		}
	protected:
		BasicProperties():
			M_NSTR(0),
			M_NLYR(0),
			M_MU(nullptr),
			M_WT(nullptr),
			M_LP_MU(nullptr),
			M_LAYER_CONSTRUCTION(SKTRAN_DO_UserSpec::LayerConstructionMethod::uniform_pressure),
			M_USE_PSEUDO_SPHERICAL(true),
			M_USE_LOS_SPHERICAL(false),
            M_USE_GREENS_FUNCTION(false),
			M_NUM_SZA(2),
			M_SS_ONLY(false),
            M_SZA_REL_SEP(0.05)
		{
			// empty
		}
		// Number of streams
		const uint				M_NSTR;

		// Number of atmospheric layers
		const uint				M_NLYR;

		// Default way to place the layers
		const SKTRAN_DO_UserSpec::LayerConstructionMethod M_LAYER_CONSTRUCTION;

		// Whether to use the pseudo-spherical approximation
		const bool M_USE_PSEUDO_SPHERICAL;

		// Whether or not to include line of sight sphericity corrections
		const bool M_USE_LOS_SPHERICAL;

        // Whether or not to use the greens function method for the particular solution
        const bool M_USE_GREENS_FUNCTION;

		// How many SZAs to use in the spherical LOS calculation
		const size_t M_NUM_SZA;

		// True if only SS terms are to be included
		const bool M_SS_ONLY;

		// Pointer to vector of stream angles.
		const VectorDim1<double>*	M_MU;

		// Pointer to vector of quadrature weights
		const VectorDim1<double>*	M_WT;

		// Legendre polynomials evaluated at stream angles. Access is [m][i][l] 
		// where m is AEOrder, i is stream index, l is polynomial order.
		const VectorDim3<LegendrePhaseContainer<NSTOKES>>*	M_LP_MU;

        const double M_SZA_REL_SEP;
	};

    template <int NSTOKES, int CNSTR=-1>
	class SolarProperties
	{
	public:
		SolarProperties(const SolarProperties* other):
			M_CSZ(other->M_CSZ),
			M_SAZ(other->M_SAZ),
			M_SOLAR_DIFFUSE_INTENSITY(other->M_SOLAR_DIFFUSE_INTENSITY),
			M_SOLAR_DIRECT_INTENSITY(other->M_SOLAR_DIRECT_INTENSITY),
			M_LP_CSZ(other->M_LP_CSZ)
		{
			// empty
		}
	protected:
		SolarProperties():
			M_CSZ(std::nan("1")),
			M_SAZ(std::nan("1")),
			M_SOLAR_DIFFUSE_INTENSITY(std::nan("1")),
			M_SOLAR_DIRECT_INTENSITY(std::nan("1")),
			M_LP_CSZ(nullptr)
		{
			 // empty
		}

		// Direct solar intensity at the TOA.
		const double								M_SOLAR_DIRECT_INTENSITY;
		// Diffuse solar intensity at the TOA.
		const double								M_SOLAR_DIFFUSE_INTENSITY;
		// Cosine of the solar zenith angle.
		const double								M_CSZ;
		// Solar azimuth angle in radians.
		const double								M_SAZ;	
		// Legendre polynomials evaluated at M_CSZ. Accessed by [m][l] where
		// m is the AEOrder and l is the polynomial order.
		LegendrePolynomials<NSTOKES>*				M_LP_CSZ;
	};

	class UserSpecProperties
	{
	public:
		UserSpecProperties(const UserSpecProperties* other):
			m_userspec(other->m_userspec)
		{
			// empty
		}
	protected:
		UserSpecProperties():
			m_userspec(nullptr)
		{
			//empty
		}
		// Pointer to user specifications.
		const SKTRAN_DO_UserSpec*			m_userspec;
	};

	class TestProperties
	{
	public:
		TestProperties(const TestProperties* other):
			m_testing(other->m_testing)
		{
			// empty
		}
	protected:
		TestProperties():
			m_testing(false)
		{
			// empty
		}
		const bool m_testing;
	};
#pragma endregion
}