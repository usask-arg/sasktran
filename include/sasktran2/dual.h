#pragma once

#include <sasktran2/internal_common.h>

namespace sasktran2 {
    enum dualstorage {
        dense,
        denseRowMajor,
        sparse,
        sparseview
    };


    /** The core storage class that defines variables in SASKTRAN2 where we need to keep track of derivatives.  A standard situation is
     *  we have a vector of values \f$y\f$ with length \f$N_y\f$ and we wish to store derivatives with respect to \f$N_d\f$
     *  quantities.  Then the derivative matrix is of size \f$(N_y, N_d)\f$.  Note that we specifically place the derivative
     *  dimension on the right since often in calculations we have a matrix operating on the left side of the values,
     *  and this naturally extends to the derivative matrix.
     *
     *  Commonly the size of the storage is known at compile time.  This can be the case for scalars, or NSTOKES vectors,
     *  or objects that depend on a templated number of DO streams.  The template parameter CSIZE allows the user to specify
     *  a constant size of values.  For dynamic sized Duals CSIZE can be set to -1.
     *
     *  Another common situation is that the derivative matrix is highly sparse.  The user can specify if the derivative
     *  matrix should be stored sparsely or densely with the DerivStorageE template parameter.
     *
     * @tparam T The underlying storage for the dual, typically double, but could be float or std::complex if necessary
     * @tparam DerivStorageE Storage type for the derivative matrix, either dense or sparse
     * @tparam CSIZE -1 if the value storage is dynamically sized, else set to the constant size of the values
     */
    template<typename T, dualstorage DerivStorageE=dualstorage::dense, int CSIZE=-1>
	class Dual {
	private:
        using DerivStorage = typename std::conditional<DerivStorageE==dualstorage::dense,
                Eigen::Matrix<T, CSIZE, -1>, typename std::conditional<DerivStorageE==dualstorage::denseRowMajor,
                Eigen::Matrix<T, CSIZE, -1, Eigen::RowMajor>,
                typename std::conditional<DerivStorageE==dualstorage::sparse,
                Eigen::SparseMatrix<T>, Eigen::SparseMatrix<T>>::type>::type>::type;

	public:
        Eigen::Vector<T, CSIZE> value; /**< values */
        DerivStorage deriv; /**< Derivative matrix of size (num_values, num_deriv) */

        Dual() {}

        /** Constructs the Dual.  Take special note of the parameter setzero which defaults to false.  By default
         *  the storage for the values and derivatives is not initialized to zero, this offers potential speed improvements
         *  if the user can initialize these directly from other values.
         *
         * @param numvalues Number of values, if CSIZE is not -1 then this must be equal to CSIZE
         * @param numderivatives Number of derivatives
         * @param setzero if true, then the derivative and value storages are initialized to 0.
         */
        Dual(int numvalues, int numderivatives, bool setzero=false) {
            resize(numvalues, numderivatives, setzero);
        }

        /** Resizes the dual storage, only possible if CSIZE is -1
         *
         * @param numvalues Number of values, if CSIZE is not -1 then this must be equal to CSIZE
         * @param numderivatives Number of derivatives
         * @param setzero if true, then the derivative and value storages are initialized to 0.
         */
        void resize(int numvalues, int numderivatives, bool setzero) {
            if constexpr(CSIZE != -1) {
                assert(CSIZE == numvalues);
            }

            if constexpr(DerivStorageE == dualstorage::sparseview) {
                // Don't have to resize, will be provided
            } else {
                deriv.resize(numvalues, numderivatives);
            }

            value.resize(numvalues);

            if(setzero) {
                value.setZero();
            }

            // If we are using dense storage for the derivatives and setzero is requested
            // then we zero out the array
            if constexpr(DerivStorageE == dualstorage::dense) {
                if(setzero) {
                    deriv.setZero();
                }
            }
        }

        /** NOT IMPLEMENTED
         *
         * @tparam T2
         * @tparam E
         * @tparam rsize
         * @param a
         * @return
         */
        template<typename T2, dualstorage E, int rsize>
        Dual<T, DerivStorageE, CSIZE>& operator+=(const Dual<T2, E, rsize>& a) {
            assert((this->value_size() == a.value_size()) || (a.value_size() == 1) );

            return *this;
        }

        /** NOT IMPLEMENTED
         *
         * @tparam T2
         * @tparam E
         * @tparam rsize
         * @param a
         * @return
         */
        template<typename T2, dualstorage E, int rsize>
        Dual<T, DerivStorageE, CSIZE>& operator*=(const Dual<T2, E, rsize>& a) {
            assert((this->value_size() == a.value_size()) || (a.value_size() == 1) );

            return *this;
        }

        /**
         *
         * @return The number of values stored
         */
        int value_size() const { return (int)value.size(); };

        /**
         *
         * @return The number of derivatives stored
         */
        int derivative_size() const { return (int)deriv.cols(); };

    };


    /** A special convenience struct to work with the derivatives of layer optical depth which is a frequently
     *  occuring situation.
     *
     *  The struct stores the optical depth of a layer, as well as an iterator into a sparse matrix that defines
     *  all of the derivatives.
     *
     */
    struct SparseODDualView {
        double od; /**< Optical depth */
        double exp_minus_od; /**< exp(-od) */
        Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator deriv_iter; /**< Iterator to a row of a sparse matrix that defines the derivatives of od */

        /**
         *
         * @param in_od optical depth
         * @param deriv_matrix A full derivative matrix for a set of optical depths
         * @param col The column of the full derivative matrix that this od value applies to
         */
        SparseODDualView(double in_od,
                         double exp_minus_od,
                         const Eigen::SparseMatrix<double, Eigen::RowMajor>& deriv_matrix,
                         int col
                         ) : od(in_od), exp_minus_od(exp_minus_od), deriv_iter(deriv_matrix, col) {};
    };

}

