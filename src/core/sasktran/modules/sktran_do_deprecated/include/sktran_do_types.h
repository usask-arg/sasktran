#pragma once
#include "modules/sktran_do_deprecated/include/sktran_do.h"
#include <boost/container/map.hpp>


// Absolute difference between the solar secant and an eigenvalue before a taylor expansion is used in the solution
#define SKTRAN_DO_GREENS_EPS 1e-4

namespace sktran_do_detail
{
#pragma region "Index aliases"

	// General purpose unsigned integer
	typedef unsigned int uint;

	// Specifies a specific order of the azimuth expansion
	typedef unsigned int AEOrder;

	// Specifies a polynomial order 
	typedef unsigned int LPOrder;

	// Specifies a solution index
	typedef unsigned int SolutionIndex;

	// Specifies a stream index
	typedef unsigned int StreamIndex;

	// Specifies an atmospheric layer index. TOA layer index is 0, ground 
	// layer index is NLYR-1.
	typedef unsigned int LayerIndex;

	// Specifies the boundary between two atmospheric layer. A boundary 
	// index, p, has a top layer with LayerIndex=p-1 and bottom layer with 
	// LayerIndex=p.
	typedef unsigned int BoundaryIndex;
#pragma endregion

#pragma region "Type aliases"
	// A Legendre polynomial
	typedef double LPoly;

	// Alias for 1 dimensional std::vector
	template<class StoredType>
	using VectorDim1 = std::vector<StoredType>;

	// Alias for 2 dimensional std::vector
	template<class StoredType>
	using VectorDim2 = std::vector<VectorDim1<StoredType>>;

	// Alias for 3 dimensional std::vector
	template<class StoredType>
	using VectorDim3 = std::vector<VectorDim2<StoredType>>;
#pragma  endregion

#pragma region "Common enumerations"
	// Specifies a location in an optical layer
	enum struct Location { CEILING, INSIDE, FLOOR };

	// Specifies a propagating direction
	enum class Propagating { UP, DOWN };

#pragma endregion

#pragma region "Linearization"
	// Derivative of the layer quantities with respect to a parameter
	// These form the input parameters to the model to calculate the derivatives for
    // The core LIDORT algorithm needs to know the changes in SSA, Optical Depth, and Legendre Coeff.
    template <int NSTOKES, int CNSTR>
	class LayerInputDerivative
	{
	public:
		LayerInputDerivative(uint nstr, LayerIndex p) {
			d_legendre_coeff.resize(nstr);
			layer_index = p;

			setZero();
		}

		std::vector<LegendreCoefficient<NSTOKES>> d_legendre_coeff;
		double d_optical_depth;
		double d_SSA;
		double d_albedo;

		LayerIndex layer_index;
        // Mapping between engine weighting functions and internal weighting functions
		std::vector<std::pair<uint, double>> group_and_triangle_fraction;
        // Engine weighting function parameters
		std::vector<std::tuple<double, double, double>> alt_and_widths;
		std::vector<double> extinctions;
		CLIMATOLOGY_HANDLE handle;
	private:
		void setZero();
	};

    template <int NSTOKES, int CNSTR>
	class InputDerivatives {
	public:
		InputDerivatives() {
            m_geometry_configured = false;
        };

		const VectorDim1<LayerInputDerivative<NSTOKES>>& layerDerivatives() const { return m_layerderivs; }
		VectorDim1<LayerInputDerivative<NSTOKES>>& layerDerivatives() { return m_layerderivs; }

		const LayerInputDerivative<NSTOKES>& operator [](int i) const { return m_layerderivs[i]; }

		void addDerivative(LayerInputDerivative<NSTOKES>&& deriv) {
			m_layerderivs.emplace_back(deriv);
		}

        inline LayerInputDerivative<NSTOKES>& addDerivative(uint nstr, LayerIndex p) {
            m_layerderivs.emplace_back(nstr, p);

            return m_layerderivs.back();
        }

		inline size_t numDerivative() const {
			return m_layerderivs.size();
		}

		inline size_t numDerivativeLayer(LayerIndex p) const {
			if (m_layerderivs.size() == 0) {
				// We are not calculating derivatives at all
				return 0;
			}
			return m_numderivlayer[p];
		}

		inline size_t layerStartIndex(LayerIndex p) const {
			if (m_layerderivs.size() == 0) {
				// Not calculating derivatives
				return 0;
			}
			return m_layerstartindex[p];
		}

		inline LayerInputDerivative<NSTOKES>& addIfNotExists(uint nstr, LayerIndex p, const CLIMATOLOGY_HANDLE& handle, size_t size_hint = 0) {
			for (uint i = 0; i < m_layerderivs.size(); ++i) {
				if (m_layerderivs[i].layer_index == p && m_layerderivs[i].handle == handle) {
					return m_layerderivs[i];
				}
			}
			m_layerderivs.emplace_back(LayerInputDerivative<NSTOKES>(nstr, p));
			auto& layerDeriv = m_layerderivs.back();
			layerDeriv.handle = handle;

            if(size_hint > 0) {
                layerDeriv.group_and_triangle_fraction.reserve(size_hint);
                layerDeriv.extinctions.reserve(size_hint);
                layerDeriv.alt_and_widths.reserve(size_hint);
            }
			return layerDeriv;
		}

		void sort(uint nlyr) {
			std::sort(std::begin(m_layerderivs), std::end(m_layerderivs), 
				[](const LayerInputDerivative<NSTOKES>& a, const LayerInputDerivative<NSTOKES>& b) -> bool {
				return a.layer_index < b.layer_index;
			});
			// Construct the layer start indicies
			m_layerstartindex.resize(nlyr);
			m_numderivlayer.resize(nlyr, 0);
			for (uint i = 0; i < m_layerderivs.size(); ++i) {
				m_numderivlayer[m_layerderivs[i].layer_index]++;
			}

			m_layerstartindex[0] = 0;
			for (uint i = 0; i < m_numderivlayer.size() - 1; i++) {
				m_layerstartindex[i + 1] = m_layerstartindex[i] + m_numderivlayer[i];
			}
		}

        bool is_geometry_configured() const { return m_geometry_configured; }
        void set_geometry_configured() { m_geometry_configured = true; }

	private:
		VectorDim1<LayerInputDerivative<NSTOKES>> m_layerderivs;
		std::vector<size_t> m_layerstartindex;
		std::vector<size_t> m_numderivlayer;

        bool m_geometry_configured;
	};


	// Derivative of a quantity with respect to the layer parameters
    template <int NSTOKES, int CNSTR=-1>
	class LayerFundamentalDerivative
	{
	public:
		LayerFundamentalDerivative(uint nstr);
		LayerFundamentalDerivative()
		{
			value = 0.0;
			d_by_opticalDepth = 0.0;
			d_by_SSA = 0.0;
		};

		void resize(uint nstr) {
			d_by_legendre_coeff.resize(nstr);
		}

		void reduce(const LayerInputDerivative<NSTOKES>& layer_deriv, double& output) const;

		LegendreCoefficientContainer<NSTOKES> d_by_legendre_coeff;
		double d_by_opticalDepth;
		double d_by_SSA;
        double value;
	};

    // A dual is a combination of a value, and the derivatives of value with respect to ALL quantities
	template <typename T>
	struct Dual {
		T value;
		Eigen::VectorX<T> deriv;

		Dual(size_t numderiv) {
			resize(numderiv);
		}
		Dual() {
			value = 0.0;
		};

		void resize(size_t numderiv, bool setzero = true) {
			deriv.resize(numderiv);
            if(setzero) {
                deriv.setZero();
                value = 0.0;
            }
		}
	};
	

    // Maps a dual defined with respect to layer quantities to the dual defined with respect to weighting function (
    // engine input) quantities
    template <int NSTOKES, typename T>
	inline Dual<T> convert_dual_to_wf(const Dual<T>& dual, const InputDerivatives<NSTOKES>& in_deriv,
		size_t numwf) {

		Dual<T> result(numwf);

		for (uint l = 0; l < in_deriv.numDerivative(); ++l)
		{
			const auto& qty = in_deriv.layerDerivatives()[l];
			for (uint k = 0; k < qty.group_and_triangle_fraction.size(); ++k) {
				result.deriv[qty.group_and_triangle_fraction[k].first] += dual.deriv[l] * qty.group_and_triangle_fraction[k].second;
			}
		}

		result.value = dual.value;
		return result;
	}

    // Operators involving Duals, note that these are usually slow because they involve reallocs so they are
    // used only sparingly
	template <typename T>
	inline Dual<T> operator+(const Dual<T>& lhs, const Dual<T>& rhs) {
		Dual<T> result;

		result.deriv.resize(lhs.deriv.size());

		result.value = lhs.value + rhs.value;
		for (uint i = 0; i < result.deriv.size(); ++i) {
			result.deriv[i] = lhs.deriv[i] + rhs.deriv[i];
		}

		return result;
	}

	template <typename T>
	inline Dual<T> operator*(const Dual<T>& lhs, const Dual<T>& rhs) {
		Dual<T> result;

		result.deriv.resize(lhs.deriv.size());

		result.value = lhs.value * rhs.value;
		for (uint i = 0; i < result.deriv.size(); ++i) {
			result.deriv[i] = lhs.deriv[i] * rhs.value + lhs.value * rhs.deriv[i];
		}

		return result;
	}

	template <typename T>
	inline Dual<T> operator/(const Dual<T>& lhs, const Dual<T>& rhs) {
		Dual<T> result;

		result.deriv.resize(lhs.deriv.size());

		result.value = lhs.value / rhs.value;
		for (uint i = 0; i < result.deriv.size(); ++i) {
			result.deriv[i] = lhs.deriv[i] / rhs.value - rhs.deriv[i] * lhs.value / (rhs.value * rhs.value);
		}

		return result;
	}

	template <typename T>
	inline Dual<T> operator*(T lhs, const Dual<T>& rhs) {
		Dual<T> result;

		result.deriv.resize(rhs.deriv.size());

		
		result.value = lhs * rhs.value;
		for (uint i = 0; i < result.deriv.size(); ++i) {
			result.deriv[i] = lhs * rhs.deriv[i];
		}

		return result;
	}

	template <typename T>
	inline Dual<T> operator*(const Dual<T>& lhs, T rhs) {
		return rhs * lhs;
	}

	template <typename T>
	inline Dual<T> operator/(const Dual<T>& lhs, T rhs) {
		return (1 / rhs) * lhs;
	}

	template <typename T>
	inline Dual<T> operator+(const Dual<T>& lhs, double rhs) {
		Dual<T> result;

		result.deriv.resize(lhs.deriv.size());

		result.value = lhs.value + rhs;
		for (uint i = 0; i < result.deriv.size(); ++i) {
			result.deriv[i] = lhs.deriv[i];
		}

		return result;
	}

	template <typename T>
	inline Dual<T> operator+(double lhs, const Dual<T>& rhs) {
		return rhs + lhs;
	}

	template <typename T>
	inline Dual<T> operator+=(Dual<T>& lhs, double rhs) {
		lhs.value += rhs;
		
		return lhs;
	}

	template <typename T>
	inline Dual<T> operator*=(Dual<T>& lhs, double rhs) {
		lhs.value *= rhs;
		for (uint i = 0; i < lhs.deriv.size(); ++i) {
			lhs.deriv[i] *= rhs;
		}

		return lhs;
	}

	template <typename T>
	inline Dual<T> operator+=(Dual<T>& lhs, const Dual<T>& rhs) {
		lhs.value += rhs.value;
		for (uint i = 0; i < lhs.deriv.size(); ++i) {
			lhs.deriv[i] += rhs.deriv[i];
		}

		return lhs;
	}

	template <typename T>
	inline Dual<T> operator-(const Dual<T>& lhs, const Dual<T>& rhs) {
		return lhs + (-1.0) * rhs;
	}

	template <typename T>
	inline Dual<T> operator-(const Dual<T>& lhs, double rhs) {
		return lhs + (-1.0) * rhs;
	}

	template <typename T>
	inline Dual<T> operator-(double lhs, const Dual<T>& rhs) {
		return lhs + (-1.0) * rhs;
	}

    // A layer dual is a Dual where the derivative is 0 except for inside a single layer.  We define a separate
    // class to represent this since operations involving LayerDuals are significantly faster
	template <typename T>
	struct LayerDual {
		T value;
		uint layer_start;
		LayerIndex layer_index;
		Eigen::VectorX<T> deriv;

		LayerDual(size_t numderiv, LayerIndex l, uint layerstartindex)
		{
			resize(numderiv);
			layer_index = l;
			layer_start = layerstartindex;
			value = 0.0;
		}

		LayerDual() {
			value = 0.0;
		};

		LayerDual<T>& operator=(const LayerDual<T>& other) {
			this->value = other.value;
			this->deriv.noalias() = other.deriv;

			return *this;
		}

		void resize(size_t numderiv) {
			deriv.resize(numderiv);
			deriv.setZero();
		}
	};

    // Operators involving Duals and LayerDuals, once again these are slow and are only used sparingly
	template <typename T>
	inline Dual<T> operator*(const Dual<T>& lhs, const LayerDual<T>& rhs) {
		Dual<T> result;

		result.deriv.resize(lhs.deriv.size());

		result.value = lhs.value * rhs.value;
		if (lhs.deriv.size() > 0) {
			result.deriv = rhs.value * lhs.deriv;
		}
		// Add in layer derivs
		if (rhs.deriv.size() > 0) {
			result.deriv(Eigen::seq(rhs.layer_start, rhs.layer_start + rhs.deriv.size() - 1)) += lhs.value * rhs.deriv;
		}

		return result;
	}

	template <typename T>
	inline Dual<T> operator*(const LayerDual<T>& lhs, const Dual<T>& rhs) {
		return rhs * lhs;
	}

	template <typename T>
	inline LayerDual<T> operator*(const LayerDual<T>& lhs, const LayerDual<T>& rhs) {
		LayerDual<T> result(lhs.deriv.size(), lhs.layer_index, lhs.layer_start);

		result.value = lhs.value * rhs.value;
		if (result.deriv.size() > 0) {
			result.deriv.noalias() = lhs.value * rhs.deriv + lhs.deriv * rhs.value;
		}

		return result;
	}

	template <typename T>
	inline LayerDual<T> operator+(const LayerDual<T>& lhs, const LayerDual<T>& rhs) {
		LayerDual<T> result(lhs.deriv.size(), lhs.layer_index, lhs.layer_start);

		result.value = lhs.value + rhs.value;
		if (result.deriv.size() > 0) {
			result.deriv.noalias() = lhs.deriv + rhs.deriv;
		}

		return result;
	}

	template <typename T>
	inline Dual<T> operator+(const Dual<T>& lhs, const LayerDual<T>& rhs) {
		Dual<T> result;

		result.deriv.resize(lhs.deriv.size());

		result.value = lhs.value + rhs.value;

		if (result.deriv.size() > 0) {
			result.deriv = lhs.deriv;
		}
		if (rhs.deriv.size() > 0) {
			result.deriv(Eigen::seq(rhs.layer_start, rhs.layer_start + rhs.deriv.size() - 1)) += rhs.deriv;
		}

		return result;
	}

	template <typename T>
	inline Dual<T> operator+(const LayerDual<T>& lhs, const Dual<T>& rhs) {
		return rhs + lhs;
	}

	template <typename T>
	inline LayerDual<T> operator*(double lhs, const LayerDual<T>& rhs) {
		LayerDual<T> result(rhs.deriv.size(), rhs.layer_index, rhs.layer_start);

		result.value = lhs * rhs.value;
		if (result.deriv.size() > 0) {
			result.deriv.noalias() = lhs * rhs.deriv;
		}

		return result;
	}

	template <typename T>
	inline LayerDual<T> operator+(double lhs, const LayerDual<T>& rhs) {
		LayerDual<T> result(rhs.deriv.size(), rhs.layer_index, rhs.layer_start);

		result.value = lhs + rhs.value;
		if (result.deriv.size() > 0) {
			result.deriv.noalias() = rhs.deriv;
		}

		return result;
	}

	template <typename T>
 	inline LayerDual<T> operator*=(LayerDual<T>& lhs, double rhs) {
		lhs.value *= rhs;
		if (lhs.deriv.size() > 0) {
			lhs.deriv *= rhs;
		}

		return lhs;
	}

	template <typename T>
	inline LayerDual<T> operator+=(LayerDual<T>& lhs, const LayerDual<T>& rhs) {
		lhs.value += rhs.value;
		if (lhs.deriv.size() > 0) {
			lhs.deriv += rhs.deriv;
		}

		return lhs;
	}

	template <typename T>
	inline LayerDual<T> operator-(const LayerDual<T>& lhs, const LayerDual<T>& rhs) {
		return lhs + (-1.0)* rhs;
	}

	template <typename T>
	inline LayerDual<T> operator+(const LayerDual<T>& lhs, double rhs) {
		LayerDual<T> result(lhs.deriv.size(), lhs.layer_index, lhs.layer_start);

		result.value = lhs.value + rhs;
		result.deriv.noalias() = lhs.deriv;

		return result;
	}

	template <typename T>
	inline LayerDual<T> operator-(const LayerDual<T>& lhs, double rhs) {
		return lhs + -1.0*rhs;
	}

	template <typename T>
	inline LayerDual<T> operator-(double lhs, const LayerDual<T>& rhs) {
		return lhs + -1.0*rhs;
	}

	template <typename T>
	inline LayerDual<T> operator/(const LayerDual<T>& lhs, double rhs) {
		return (1 / rhs) * lhs;
	}

	template <typename T>
	inline LayerDual<T> operator/(const LayerDual<T>& lhs, const LayerDual<T>& rhs) {
		LayerDual<T> result(lhs.deriv.size(), lhs.layer_index, rhs.layer_start);

		result.value = lhs.value / rhs.value;
		if (result.deriv.size() > 0) {
			result.deriv.noalias() = lhs.deriv / rhs.value - rhs.deriv * (lhs.value / (rhs.value*rhs.value));
		}

		return result;
	}

	template <typename T>
	inline Dual<T> operator*=(Dual<T>& lhs, const LayerDual<T>& rhs) {
		if (lhs.deriv.size() > 0) {
			lhs.deriv *= rhs.value;
		}

		if (rhs.deriv.size() > 0) {
			lhs.deriv(Eigen::seq(rhs.layer_start, rhs.layer_start + rhs.deriv.size() - 1)) += lhs.value * rhs.deriv;
		}
		lhs.value *= rhs.value;

		return lhs;
	}

	template <typename T>
	inline VectorLayerDual<T> operator*(const VectorLayerDual<T>& lhs, const LayerDual<T>& rhs) {
		VectorLayerDual<T> result(lhs);

		result.value = lhs.value * rhs.value;

		if (result.deriv.size() > 0) {
			result.deriv = rhs.value * lhs.deriv;
		}
		
		for (uint i = 0; i < result.deriv.rows(); ++i) {
			for (uint j = 0; j < result.deriv.cols(); ++j) {
				result.deriv(i, j) += lhs.value(j)*rhs.deriv(i);
			}
		}

		return result;
	}

	namespace dual {
		// Functions operating on duals
		template <typename T>
		inline sktran_do_detail::Dual<T> exp(sktran_do_detail::Dual<T> x) {
			sktran_do_detail::Dual<T> result;

			result.deriv.resize(x.deriv.size());

			result.value = std::exp(x.value);
			for (uint i = 0; i < x.deriv.size(); ++i) {
				result.deriv[i] = std::exp(x.value) * x.deriv[i];
			}

			return result;
		}
	}

	namespace layerdual {
		// Functions operating on layerduals
		template <typename T>
		inline sktran_do_detail::LayerDual<T> exp(sktran_do_detail::LayerDual<T> x) {
			sktran_do_detail::LayerDual<T> result(x.deriv.size(), x.layer_index, x.layer_start);

			result.value = std::exp(x.value);
			if (result.deriv.size() > 0) {
				result.deriv = std::exp(x.value) * x.deriv;
			}

			return result;
		}
	}



#pragma endregion


#pragma region "Math helpers"
	// Kronecker delta function.
	// Returns : 1.0 if i == j, else return 0.0
	inline double kronDelta(uint i, uint j)
	{
		return (i == j) ? 1.0 : 0.0;
	}
	
	static const double PI = EIGEN_PI;

	// Used to efficiently calculate the triple product of vectors of 
	// legendre polynomials.
	template <int NSTOKES, int CNSTR=-1>
	class LPTripleProduct
	{
	public:
		// Calculates an auxilary result for the given Legendre polynomial 
		// vectors.
		LPTripleProduct(AEOrder m,
			const std::vector<LegendreCoefficient<NSTOKES>>& lephase,
			const std::vector<LegendrePhaseContainer<NSTOKES>>& lp1,
			const std::vector<LegendrePhaseContainer<NSTOKES>>& lp2) :
			m_association_order(m),
			m_aux(std::piecewise_construct, std::forward_as_tuple((uint) lephase.size()), std::forward_as_tuple((uint) lephase.size())),
            m_nstr(lephase.size())
		{
			calculate(lephase, lp1, lp2);
		}

		LPTripleProduct(AEOrder m, uint size) : 
			m_association_order(m),
			m_aux(std::piecewise_construct, std::forward_as_tuple(size), std::forward_as_tuple(size)),
            m_nstr(size)
		{}

        LPTripleProduct(uint size) :
                m_aux(std::piecewise_construct, std::forward_as_tuple(size), std::forward_as_tuple(size)),
                m_nstr(size)
        {}

        void calculate(AEOrder m, const std::vector<LegendreCoefficient<NSTOKES>>& lephase, const std::vector<LegendrePhaseContainer<NSTOKES>>& lp1, const std::vector<LegendrePhaseContainer<NSTOKES>>& lp2) {
            m_association_order = m;
            calculate(lephase, lp1, lp2);
        }


        void calculate(const std::vector<LegendreCoefficient<NSTOKES>>& lephase, const std::vector<LegendrePhaseContainer<NSTOKES>>& lp1, const std::vector<LegendrePhaseContainer<NSTOKES>>& lp2);

		// Calculate the result of the triple product given the number of 
		// negations.
		inline Eigen::Matrix<double, NSTOKES, NSTOKES> negations(uint num) {
			if(num % 2 == 0) {
				return Eigen::Matrix<double, NSTOKES, NSTOKES>(m_aux.first.value);
			} else {
				return Eigen::Matrix<double, NSTOKES, NSTOKES>(m_aux.second.value);
			}
		}

		inline sktran_do_detail::TripleProductDerivativeHolder<NSTOKES> negations_derivative(uint num) {
            sktran_do_detail::TripleProductDerivativeHolder<NSTOKES> result(m_nstr);

            negations_derivative_emplace(num, result);

			return result;
		}

		void negations_derivative_emplace(uint num, sktran_do_detail::TripleProductDerivativeHolder<NSTOKES>& holder);

	private:
		// A pair of the result of the triple product of even and odd order
		// polynomials. Triple product of even polynomials are stored in first,
		// the triple product of odd polynomials are stored in second.
		using AuxilaryResult = std::pair<sktran_do_detail::TripleProductDerivativeHolder<NSTOKES>, sktran_do_detail::TripleProductDerivativeHolder<NSTOKES>>;

		AuxilaryResult m_aux;
		AEOrder m_association_order;
        uint m_nstr;
	};


#pragma  endregion

#pragma region "Useful DO classes"
	// Storage object for simplified line of sight. Efficiencies can be made in
	// the discrete ordinates algorithm if lines of sight are sorted by cosine 
	// of zenith angles.
	struct LineOfSight
	{
		double coszenith;
		double azimuth;
		double cos_scattering_angle;
		std::vector<double>* wf;
		double* radiance;
		uint unsorted_index;
		double observeraltitude;


		// Sorts a vector of lines of sight while maintaining the proper 
		// radiance and weight function return orders.
		inline static void sort(VectorDim1<LineOfSight>& los, VectorDim1<double>& radiance, VectorDim2<double>* wf, int NSTOKES)
		{
			for(uint ridx = 0; ridx < los.size(); ++ridx) {
				los[ridx].radiance = &radiance[ridx*NSTOKES];
				los[ridx].wf = (wf != nullptr) ? &((*wf)[ridx]) : nullptr;
				los[ridx].unsorted_index = ridx;
			}
			std::sort(los.begin(), los.end(), [](const LineOfSight& lhs, const LineOfSight& rhs) {
				return lhs.coszenith > rhs.coszenith;
			});
		}
	};

	// Contains the homogeneous and particular solutions for a single layer and azimuth order.
	// Note that both the homogeneous and particular solutions
	// are decomposed in sum/difference operators to reduce the order of the problem from NSTR to NSTR/2
    template <int NSTOKES, int CNSTR=-1>
	class RTEGeneralSolution
	{

	// If NSTOKES == 1 all of the eigenvalues are real, but for NSTOKES > 1 we can
	// have complex eigenvalues and complex homogeneous solution vectors
	using HomogType = typename std::conditional<NSTOKES != 5, double, std::complex<double>>::type;

    using HomogVector = typename std::conditional<CNSTR == -1, Eigen::VectorX<HomogType>, Eigen::Vector<HomogType, CNSTR/2 * NSTOKES>>::type;
    using HomogMatrix = typename std::conditional<CNSTR == -1, Eigen::VectorX<HomogType>, Eigen::Vector<HomogType, CNSTR/2 * NSTOKES * CNSTR/2 * NSTOKES>>::type;

    using Vector = typename std::conditional<CNSTR == -1, Eigen::VectorX<double>, Eigen::Vector<double, CNSTR/2 * NSTOKES>>::type;

    using HomogVectorLayerDual = typename std::conditional<CNSTR == -1, VectorLayerDual<HomogType>, VectorLayerDual<HomogType, CNSTR/2 * NSTOKES * CNSTR/2 * NSTOKES>>::type;
    using EigVectorLayerDual = typename std::conditional<CNSTR == -1, VectorLayerDual<HomogType>, VectorLayerDual<HomogType, CNSTR/2 * NSTOKES>>::type;

    using VectorLayerDualType = typename std::conditional<CNSTR == -1, VectorLayerDual<double>, VectorLayerDual<double, CNSTR/2 * NSTOKES>>::type;
    using VectorDualType = typename std::conditional<CNSTR == -1, VectorDual<double>, VectorDual<double, CNSTR/2 * NSTOKES>>::type;


    public:
		RTEGeneralSolution():
			M_NSTR(0)
		{
			m_use_green_function = false;
		}
		void resize(size_t nstr, size_t nderivlayer, LayerIndex l, uint layer_start, uint nderivtotal) {
			const_cast<uint&>(M_NSTR) = (uint)nstr;
			m_eigval.resize(nstr/2 * NSTOKES, nderivlayer, l, layer_start);

			m_homog_plus.resize((nstr / 2) * (nstr / 2) * NSTOKES * NSTOKES, nderivlayer, l, layer_start);
			m_homog_minus.resize((nstr / 2) * (nstr / 2) * NSTOKES * NSTOKES, nderivlayer, l, layer_start);
			m_particular_plus.resize(nstr / 2 * NSTOKES, nderivtotal);
			m_particular_minus.resize(nstr / 2 * NSTOKES, nderivtotal);

            m_green_A_minus.resize( nstr / 2 * NSTOKES, nderivlayer, l, layer_start);
            m_green_A_plus.resize( nstr / 2 * NSTOKES, nderivlayer, l, layer_start);

            m_Gplus_top.resize( nstr / 2 * NSTOKES, nderivtotal);
            m_Gplus_bottom.resize( nstr / 2 * NSTOKES, nderivtotal);
            m_Gminus_top.resize( nstr / 2 * NSTOKES, nderivtotal);
            m_Gminus_bottom.resize( nstr / 2 * NSTOKES, nderivtotal);

		}
		inline HomogType eigval(SolutionIndex a) const {
			return m_eigval.value(a);
		}
		inline HomogType homog_plus(StreamIndex j, SolutionIndex a) const {
			return m_homog_plus.value(j + (M_NSTR*NSTOKES / 2) * a);
		}
		inline HomogType homog_minus(StreamIndex j, SolutionIndex a) const {
			return m_homog_minus.value(j + (M_NSTR*NSTOKES / 2) * a);
		}

		inline HomogType d_homog_minus(StreamIndex j, SolutionIndex a, uint deriv) const {
			return m_homog_minus.deriv(deriv, j + (M_NSTR*NSTOKES / 2) * a);
		}

		inline HomogType d_homog_plus(StreamIndex j, SolutionIndex a, uint deriv) const {
			return m_homog_plus.deriv(deriv, j + (M_NSTR*NSTOKES / 2) * a);
		}

		inline double particular_plus(StreamIndex j) const {
			return m_particular_plus.value(j);
		}
		inline double particular_minus(StreamIndex j) const {
			return m_particular_minus.value(j);
		}

		inline HomogVector& eigval() {
			return m_eigval.value;
		}
		inline const HomogVector& eigval() const {
			return m_eigval.value;
		}

		inline EigVectorLayerDual& dual_eigval() {
			return m_eigval;
		}
		inline const EigVectorLayerDual& dual_eigval() const {
			return m_eigval;
		}

		inline HomogMatrix& homog_plus() {
			return m_homog_plus.value;
		}
		inline const HomogMatrix& homog_plus() const {
			return m_homog_plus.value;
		}

		inline HomogVectorLayerDual& dual_homog_plus() {
			return m_homog_plus;
		}
		inline const HomogVectorLayerDual& dual_homog_plus() const {
			return m_homog_plus;
		}

		inline HomogType* homog_plus(SolutionIndex a) {
			return m_homog_plus.value.data() + (M_NSTR*NSTOKES / 2) * a;
		}
		inline const HomogType* homog_plus(SolutionIndex a) const {
			return m_homog_plus.value.data() + (M_NSTR*NSTOKES / 2) * a;
		}

		inline HomogMatrix& homog_minus() {
			return m_homog_minus.value;
		}
		inline const HomogMatrix& homog_minus() const {
			return m_homog_minus.value;
		}

		inline HomogVectorLayerDual& dual_homog_minus() {
			return m_homog_minus;
		}
		inline const HomogVectorLayerDual& dual_homog_minus() const {
			return m_homog_minus;
		}

		inline HomogType* homog_minus(SolutionIndex a) {
			return m_homog_minus.value.data() + (M_NSTR*NSTOKES / 2) * a;
		}
		inline const HomogType* homog_minus(SolutionIndex a) const {
			return m_homog_minus.value.data() + (M_NSTR*NSTOKES / 2) * a;
		}

		inline Vector& particular_plus() {
			return m_particular_plus.value;
		}
		inline const Vector& particular_plus() const {
			return m_particular_plus.value;
		}

		inline VectorDualType& dual_particular_plus() {
			return m_particular_plus;
		}
		inline const VectorDualType& dual_particular_plus() const {
			return m_particular_plus;
		}

		inline Vector& particular_minus() {
			return m_particular_minus.value;
		}
		inline const Vector& particular_minus() const {
			return m_particular_minus.value;
		}

		inline VectorDualType& dual_particular_minus() {
			return m_particular_minus;
		}
		inline const VectorDualType& dual_particular_minus() const {
			return m_particular_minus;
		}

        inline VectorLayerDualType& dual_green_A_minus() {
            return m_green_A_minus;
        }

        inline const VectorLayerDualType& dual_green_A_minus() const {
            return m_green_A_minus;
        }

        inline VectorLayerDualType& dual_green_A_plus() {
            return m_green_A_plus;
        }

        inline const VectorLayerDualType& dual_green_A_plus() const {
            return m_green_A_plus;
        }

        inline VectorDualType& dual_Gplus_top() {
            return m_Gplus_top;
        }

        inline const VectorDualType& dual_Gplus_top() const {
            return m_Gplus_top;
        }

        inline VectorDualType& dual_Gplus_bottom() {
            return m_Gplus_bottom;
        }

        inline const VectorDualType& dual_Gplus_bottom() const {
            return m_Gplus_bottom;
        }

        inline VectorDualType& dual_Gminus_top() {
            return m_Gminus_top;
        }

        inline const VectorDualType& dual_Gminus_top() const {
            return m_Gminus_top;
        }

        inline VectorDualType& dual_Gminus_bottom() {
            return m_Gminus_bottom;
        }

        inline const VectorDualType& dual_Gminus_bottom() const {
            return m_Gminus_bottom;
        }

		inline const uint nstr() const { return M_NSTR; }

        void set_use_green_function( bool use ) {
            m_use_green_function = use;
        }

        bool use_green_function() const {
            return m_use_green_function;
        }

	private:
		const uint M_NSTR;

        EigVectorLayerDual m_eigval;            // k,  See eq (16)

        HomogVectorLayerDual  m_homog_plus;	     // W+, See eq (17)
        HomogVectorLayerDual  m_homog_minus;       // W-, See eq (17)
        VectorDualType m_particular_plus;		 // Z+, See eq (23)
        VectorDualType m_particular_minus;		 // Z-, See eq (23)

        VectorLayerDualType m_green_A_plus;
        VectorLayerDualType m_green_A_minus;

        VectorDualType m_Gplus_top;
        VectorDualType m_Gplus_bottom;
        VectorDualType m_Gminus_top;
        VectorDualType m_Gminus_bottom;

        bool m_use_green_function;
	};

	// Stores some cached quantities that are needed between different steps
	// of the RTE solution
    template <int NSTOKES, int CNSTR=-1>
	class RTESolutionCache
	{
        using Matrix = typename std::conditional<CNSTR == -1, Eigen::MatrixXd, Eigen::Matrix<double, CNSTR/2 * NSTOKES, CNSTR/2 * NSTOKES>>::type;

		public:
			RTESolutionCache() :
				M_NSTR(0)
			{
				// empty
			}
			RTESolutionCache(uint nstr) :
				M_NSTR(nstr),
				m_s_plus(nstr/2 * NSTOKES, nstr/2 * NSTOKES),
				m_s_minus(nstr/2 * NSTOKES, nstr/2 * NSTOKES),
				m_eigmtx(nstr/2 * NSTOKES, nstr/2 * NSTOKES)
			{
				// empty
			}
			void resize(uint nstr) {
				const_cast<uint&>(M_NSTR) = nstr;
				m_s_plus.resize(nstr/2 * NSTOKES, nstr/2 * NSTOKES);
				m_s_minus.resize(nstr/2 * NSTOKES, nstr/2 * NSTOKES);
				m_eigmtx.resize(nstr/2 * NSTOKES, nstr/2 * NSTOKES);
			}

			inline Matrix& s_plus() {
				return m_s_plus;
			}
			inline const Matrix& s_plus() const {
				return m_s_plus;
			}

			inline Matrix& s_minus() {
				return m_s_minus;
			}
			inline const Matrix& s_minus() const {
				return m_s_minus;
			}

			inline Matrix& eigmtx() {
				return m_eigmtx;
			}
			inline const Matrix& eigmtx() const {
				return m_eigmtx;
			}

		private:
			const uint M_NSTR;

			Matrix m_s_plus;
			Matrix m_s_minus;
			Matrix m_eigmtx;
	};

	// Stores the homogeneous coefficients L, M for a single layer solution
	struct RTEBoundarySolution
	{
		VectorDual<double>	L_coeffs;
		VectorDual<double>  M_coeffs;
	};


    template <int NSTOKES, int CNSTR=-1>
	struct LayerSolution
	{
		using HomogType = typename std::conditional<NSTOKES != 5, double, std::complex<double>>::type;

		void configure(size_t nstr, LayerIndex l, const InputDerivatives<NSTOKES>& input_deriv) {
			// Core solution contains no cross derivatives
			value.resize(nstr, input_deriv.numDerivativeLayer(l), l, (uint)input_deriv.layerStartIndex(l),
				(uint)input_deriv.numDerivative());
			
			// Cache uses a different format and is a little messy 
			cache.resize((uint)nstr);

			boundary.L_coeffs.resize(nstr / 2 * NSTOKES, input_deriv.numDerivative());
			boundary.M_coeffs.resize(nstr / 2 * NSTOKES, input_deriv.numDerivative());

			layer_index = l;

			configureDerivatives(nstr, layer_index, input_deriv);
		}

		void configureDerivatives(size_t nstr, LayerIndex layer_index, const InputDerivatives<NSTOKES>& input_deriv)
		{
			if (input_deriv.numDerivative() == 0) {
				// We have no derivatives, just return
				return;
			}

			// Terms that don't contain cross derivatives
			auto numderiv = input_deriv.numDerivativeLayer(layer_index);
			d_cache.resize(numderiv);
			for (uint i = 0; i < numderiv; ++i) {
				d_cache[i].resize((uint)nstr);
			}
		}

		// TODO: Remove all of these accessor functions below, they shouldn't really be used anymore
		Dual<double> dual_particular_plus(const InputDerivatives<NSTOKES>& in_deriv, StreamIndex j) const {
			Dual<double> result(in_deriv.numDerivative());

			result.value = value.particular_plus(j);
			result.deriv = value.dual_particular_plus().deriv(Eigen::all, j);

			return result;
		}

		Dual<double> dual_particular_minus(const InputDerivatives<NSTOKES>& in_deriv, StreamIndex j) const {
			Dual<double> result(in_deriv.numDerivative());

			result.value = value.particular_minus(j);
			result.deriv = value.dual_particular_minus().deriv(Eigen::all, j);

			return result;
		}

		LayerDual<HomogType> dual_homog_minus(const InputDerivatives<NSTOKES>& in_deriv, StreamIndex j, SolutionIndex a) const {
			LayerDual<HomogType> result(in_deriv.numDerivativeLayer(layer_index), layer_index, (uint)in_deriv.layerStartIndex(layer_index));

			result.value = value.homog_minus(j, a);
			result.deriv = value.dual_homog_minus().deriv(Eigen::all, j + (value.nstr() / 2) * a);

			return result;
		}

		LayerDual<HomogType> dual_homog_plus(const InputDerivatives<NSTOKES>& in_deriv, StreamIndex j, SolutionIndex a) const {
			LayerDual<HomogType> result(in_deriv.numDerivativeLayer(layer_index), layer_index, (uint)in_deriv.layerStartIndex(layer_index));

			result.value = value.homog_plus(j, a);
			result.deriv = value.dual_homog_plus().deriv(Eigen::all, j + (value.nstr() / 2) * a);

			return result;
		}

		LayerIndex layer_index;
		RTEGeneralSolution<NSTOKES, CNSTR> value;
		RTESolutionCache<NSTOKES, CNSTR> cache;
		std::vector<RTESolutionCache<NSTOKES, CNSTR>> d_cache;
		RTEBoundarySolution boundary;
	};

#pragma endregion

	struct LOSDiagnostics
	{
		/// Reference point. Atmospheric properties are pulled from the vertical column above the reference point.
		GEODETIC_INSTANT reference_point;
		/// Zenith angle of the line of sight at the reference point [radians].
		double local_look_zentih;
		/// Solar zenith angle at the reference point [radians].
		double local_solar_zenith;
		/// Relative difference in azimuth between the line of sight and the sun [radians].
		double local_relative_azimuth;
	};
	struct DiscreteAtmosphereDiagnostics
	{
		/// Optical depths of discrete-atmosphere layers. TOA ~ index=0. layer_optical_depths[layer]
		std::vector<double> layer_optical_depths;
		/// Single-scatter-albedo of discrete-atmosphere layers. TOA ~ index=0. layer_ssa[layer]
		std::vector<double> layer_ssa;
		/// Phase function expansion terms for discrete-atmosphere layers. TOA ~ index=0. layer_phasef_expansion[layer][m]
		std::vector<std::vector<std::vector<double>>> layer_phasef_expansion;

		std::vector<double> layer_boundary_altitudes;
	};
	struct SingleRTSDiagnostic
	{
		SingleRTSDiagnostic(uint nstr, uint nlyr) {
			intensity_components.resize(nstr, 0.0);
			reflection_components.resize(nstr, 0.0);
			layer_participating_source_term.resize(nlyr, std::vector<double>(nstr, 0.0));
		}
		/// Computed intensity (at TOA in LOS)
		double intensity; //
		/// Components of the azimuth expansion (at TOS in LOS) I = cos(m*az_diff) * I_components. intensity_components[m = azimuthexpansion index]
		std::vector<double> intensity_components;
		/// Components of the ground-reflected azimuth expansion (at ground in LOS). reflection_components[m = azimuthexpansion index]
		std::vector<double> reflection_components; //
		/// Participating source term for discrete-atmosphere layers. This is the source-term integrated through the given layer. layer_participating_source_term[layer][m]
		std::vector<std::vector<double>> layer_participating_source_term; //
		/// Exit/observer optical depth
		double observer_od;
	};

	struct RTSDiagnostics
	{
		void initialize(uint p_nstr, uint p_nlyr, uint p_nlos, uint nptrbs) {
			nlos = p_nlos;
			nstr = p_nstr;
			nlyr = p_nlyr;

			los_diagnostic.resize(nlos, SingleRTSDiagnostic(nstr, nlyr));
			atmo_diagnostics.layer_optical_depths.resize(nlyr);
			atmo_diagnostics.layer_ssa.resize(nlyr);
			atmo_diagnostics.layer_phasef_expansion.resize(nlyr, std::vector<std::vector<double>>(nstr, std::vector<double>(nstr)));
			atmo_diagnostics.layer_boundary_altitudes.resize(nlyr + 1);

			ptrb_altitudes.reserve(nptrbs*nlyr);
			ptrb_heights.reserve(nptrbs*nlyr);
			ptrb_eps.reserve(nptrbs*nlyr);
			ptrb_ssa_qty.reserve(nptrbs*nlyr);
			ptrb_optd_qty.reserve(nptrbs*nlyr);
			ptrb_group.reserve(nptrbs*nlyr);
		}

		void push_back_ptrb() {
			ptrb_altitudes.push_back(std::nan("1"));
			ptrb_heights.push_back(std::nan("1"));
			ptrb_eps.push_back(std::nan("1"));
			ptrb_ssa_qty.push_back(std::nan("1"));
			ptrb_optd_qty.push_back(std::nan("1"));
			ptrb_group.push_back(std::nan("1"));

		}
		uint nlos;
		uint nstr;
		uint nlyr;

		/// Diagnostics for the discretized atmosphere
		DiscreteAtmosphereDiagnostics atmo_diagnostics;
		/// Diagnostic for lines of sight. Index is line of sight index.
		std::vector<SingleRTSDiagnostic> los_diagnostic;

		std::vector<double> ptrb_altitudes;
		std::vector<double> ptrb_heights;
		std::vector<double> ptrb_eps;
		std::vector<double> ptrb_ssa_qty;
		std::vector<double> ptrb_optd_qty;
		std::vector<double> ptrb_group;
	};

    struct LegendreSumMatrixCommon {
        uint linear_index(StreamIndex m, StreamIndex j) const {
            uint first = m;
            uint second = j;
            uint offset;

            if (first < N && second < N) {
                // Quadrant one
                offset = 0;
                if (first > second) {
                    std::swap(first, second);
                }
            }
            else if (first >= N && second >= N) {
                // Quadrant 4
                offset = 0;
                first -= N;
                second -= N;
                if (first > second) {
                    std::swap(first, second);
                }
            }
            else {
                offset = N * N * 2;
                if (first >= N) {
                    first -= N;
                }
                if (second >= N) {
                    second -= N;
                }
                if (first > second) {
                    std::swap(first, second);
                }
            }
            return offset + first * N - ((first + 1)*first) / 2 + second;

        }

        uint N;
        double M_SSA;
    };

    template <int NSTOKES, int CNSTR=-1>
    struct LegendreSumMatrixStorage : LegendreSumMatrixCommon {
        std::vector<sktran_do_detail::TripleProductDerivativeHolder<NSTOKES>> storage;

        // This is actually unique ownership but we can't use a smart ptr because it is inside a VectorDim2
        // where the internal VectorDim1 is default copy? But it is never copied or moved so we just treat it as
        // unique
        LPTripleProduct<NSTOKES>* triple_product;

        LegendreSumMatrixStorage() : triple_product(nullptr) {}

        ~LegendreSumMatrixStorage() {
            if(triple_product != nullptr) {
                delete triple_product;
            }
        }

        void resize(uint nstr) {
            if(storage.size() > 0)
                return;
            const uint nval = nstr * nstr;
            N = nstr / 2;

            storage.resize(nval, nstr);
            if(triple_product == nullptr) {
                triple_product = new LPTripleProduct<NSTOKES>(nstr);
            }
        }

        TripleProductDerivativeHolder<NSTOKES> dual(StreamIndex m, StreamIndex j) const {
            const uint index = linear_index(m, j);

            return storage[index];
        }

        void inplace_dual(StreamIndex m, StreamIndex j, TripleProductDerivativeHolder<NSTOKES>& dual) const {
            const uint index = linear_index(m, j);

            dual = storage[index];
        }
    };

	class InvalidConfiguration: public std::exception
	{
	public:
		explicit InvalidConfiguration(const char* message):
			m_msg(message)
		{
		}
		explicit InvalidConfiguration(const std::string& message):
			m_msg(message)
		{
		}
		virtual ~InvalidConfiguration() throw () {}
		virtual const char* what() const throw() {
			return m_msg.c_str();
		}
	private:
		std::string m_msg;
	};
	
	class InternalError: public std::exception
	{
	public:
		explicit InternalError(const char* message = "<none>")
		{
			m_msg = "An unexpected internal exception was thrown. This is likely a bug! "
				"Please submit a issue at: https://arggit.usask.ca/ARGPackages/SasktranDO. "
				"The following error message was given: " + std::string(message);
		}
		explicit InternalError(const std::string& message = "<none>")
		{
			m_msg = "An unexpected internal exception was thrown. This is likely a bug! "
				"Please submit a issue at: https://arggit.usask.ca/ARGPackages/SasktranDO. "
				"The following error message was given: " + std::string(message);
		}
		virtual ~InternalError() throw () {}
		virtual const char* what() const throw() {
			return m_msg.c_str();
		}
	protected:
		std::string m_msg;
	};

	class InternalRuntimeError: public InternalError
	{
	public:
		explicit InternalRuntimeError(const char* message):
			InternalError(message)
		{
			m_msg = "An internal runtime exception has occured. Likey due to insufficient precision. ERROR MESSAGE:" + std::string(message);
		}
		explicit InternalRuntimeError(const std::string& message):
			InternalError(message)
		{
			m_msg = "An internal runtime exception has occured. Likey due to insufficient precision. ERROR MESSAGE:" + std::string(message);
		}
		virtual ~InternalRuntimeError() throw () {}
		virtual const char* what() const throw() {
			return m_msg.c_str();
		}
	};

	class SASKTRANAtmosphereInterface {
	public:
		SASKTRANAtmosphereInterface() {
			m_brdf = nullptr;
			m_atmospheric_state = new skClimatology_MSIS90;
			m_atmospheric_state->AddRef();
		}

		~SASKTRANAtmosphereInterface() {
			if (m_brdf != nullptr)
			{
				m_brdf->Release();
			}
			if (m_atmospheric_state != nullptr)
			{
				m_atmospheric_state->Release();
			}
		}

		void add_species(const CLIMATOLOGY_HANDLE& species, skClimatology* numberdensityclimatology, skOpticalProperties* particleopticalprops) {
			m_species.emplace_back(SKTRAN_AtmosphericOpticalStateEntry_V21(species));
			m_species.back().Configure(species, numberdensityclimatology, particleopticalprops);
		}

		void set_albedo(double albedo) {
			m_brdf = new SKTRAN_BRDF_Lambertian(albedo);
			m_brdf->AddRef();
		}

		void set_albedo(skBRDF* brdf) {
			m_brdf = brdf;
			m_brdf->AddRef();
		}

		const skBRDF* albedo() {
			return m_brdf;
		}

		void set_atmospheric_state(skClimatology* clim) {
			m_atmospheric_state->Release();
			m_atmospheric_state = clim;
			m_atmospheric_state->AddRef();
		}

		skClimatology* atmospheric_state() {
			return m_atmospheric_state;
		}

		std::vector<SKTRAN_AtmosphericOpticalStateEntry_V21>& species() {
			return m_species;
		}


	private:
		std::vector<SKTRAN_AtmosphericOpticalStateEntry_V21> m_species;
		skBRDF* m_brdf;
		skClimatology* m_atmospheric_state;
	};
}
