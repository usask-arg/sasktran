#pragma once
//#include "sktran_disco/sktran_do.h"

namespace sasktran_disco
{
	// Base storage for the value/derivative pairs, used for terms that contain cross derivatives
	// Unfortunately some terms we want to store are complex
	template <typename T, int CSIZE=-1>
	struct VectorDual {
		Eigen::Vector<T, CSIZE> value;
		Eigen::Matrix<T, -1, CSIZE> deriv;

		VectorDual(size_t numvalue, size_t numderiv) {
			resize(numvalue, numderiv);
		}

		VectorDual() {};

		void resize(size_t numvalue, size_t numderiv) {
            if(numvalue != value.size()) {
                value.resize(numvalue);
            }
            if(numvalue != deriv.cols() || numderiv != deriv.rows()) {
                deriv.resize(numderiv, numvalue);
            }
		}
	};

	// Base storage for value/derivative pairs that do not contain cross derivatives
	// Essentially the same storage but also contains the layer indicies
	template <typename T, int CSIZE=-1>
	struct VectorLayerDual : public VectorDual<T, CSIZE> {
		unsigned int layer_index;
		unsigned int layer_start;

		VectorLayerDual(size_t numvalue, size_t numderiv, unsigned int l, unsigned int lstart) :
			VectorDual<T>(numvalue, numderiv),
			layer_index(l),
			layer_start(lstart)
		{
		}

		VectorLayerDual(const VectorLayerDual& other) {
			resize(other.value.size(), other.deriv.rows(), other.layer_index, other.layer_start);
		}

		VectorLayerDual() {};

		void resize(size_t numvalue, size_t numderiv, unsigned int l, unsigned int lstart) {
			VectorDual<T, CSIZE>::resize(numvalue, numderiv);

			layer_index = l;
			layer_start = lstart;
		}
	};

	template <typename T>
	inline VectorDual<T> operator*=(VectorDual<T>& lhs, double rhs) {
		lhs.value *= rhs;
		lhs.deriv *= rhs;

		return lhs;
	}

	template <typename T>
	inline VectorDual<T> operator*=(VectorDual<T>& lhs, const std::vector<double>& rhs) {
		for (int i = 0; i < lhs.value.size(); ++i) {
			lhs.value(i) *= rhs[i];
			lhs.deriv(Eigen::all, i) *= rhs[i];
		}
		return lhs;
	}

}