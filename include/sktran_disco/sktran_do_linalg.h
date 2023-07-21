#pragma once
#include "sktran_disco/sktran_do.h"


namespace sasktran_disco
{
	// Linear algebra namespace.
	namespace la
	{
		// Used to index a submatrix of some larger matrix. Essentially just 
		// calls  mat(i + tl_i, j + tl_j) on any index operation (i, j).
		template<class MatrixClass>
		class MatrixBlock
		{
		public:
			MatrixBlock(MatrixClass& mat, uint tl_i, uint tl_j, uint rows, uint cols):
				M_I0(tl_i),
				M_J0(tl_j),
				M_ROWS(rows),
				M_COLS(cols),
				M_IS_MATRIX(true),
				m_mat(mat)
			{
				// empty
			}
			MatrixBlock(MatrixClass& mat, uint tl_i, uint rows, uint cols):
				M_I0(tl_i),
				M_ROWS(rows),
				M_COLS(cols),
				M_IS_MATRIX(false),
				m_mat(mat)
			{
				// empty
			}
			inline double& operator()(StreamIndex i, SolutionIndex j) {
				return m_mat(M_I0 + i, M_J0 + j);
			}
			inline double operator()(StreamIndex i, SolutionIndex j) const {
				return m_mat(M_I0 + i, M_J0 + j);
			}
			inline double& operator[](StreamIndex i) {
				return m_mat[M_I0 + i];
			}
			inline double operator[](StreamIndex i) const {
				return m_mat[M_I0 + i];
			}
		private:
			const bool M_IS_MATRIX;
			const uint M_I0, M_J0, M_ROWS, M_COLS;
			MatrixClass& m_mat;
		};

		// Used for storing the boundary-value-problem (BVP) matrix. This
		// matrix is block-tridiagonal and stored according to the LAPACK
		// band storage scheme.
        template <int NSTOKES, int CNSTR=-1>
		class BVPMatrix
		{
		public: // Constructors and destructor
			BVPMatrix(uint NSTR, uint NLYR):
				M_NCD(3 * NSTR * NSTOKES / 2 - 1),
				M_NCOLS(NSTR*NSTOKES*NLYR),
				M_NROWS(3 * (3 * NSTR*NSTOKES / 2 - 1) + 1),
				M_NRM1(3 * (3 * NSTR*NSTOKES / 2 - 1)),
				M_NSTR(NSTR),
				M_NLYR(NLYR)
			{
				m_data = std::unique_ptr<double[]>(new double[M_NROWS * M_NCOLS]);
			}

            ~BVPMatrix(){};

		public: // Interface
			inline double& operator()(uint i, uint j) {
				return m_data[2 * M_NCD + i + j*M_NRM1];
			}
			inline double operator()(uint i, uint j) const {
				return m_data[2 * M_NCD + i + j*M_NRM1];
			}
			inline void setZero() const {
				std::fill_n(m_data.get(), M_NCOLS*M_NROWS, 0.0);
			}

			inline void multiplyVector(const Eigen::VectorXd& rhs, Eigen::VectorXd& result) {
				// TODO: Fix this
				assert(result.size() == rhs.size());
				assert(result.size() == (M_NSTR * M_NLYR));

				result.setZero();
			
				// We go through every non-zero element of the data array
				for (int i = 2*M_NCD; i < int(M_NROWS * M_NCOLS); ++i) {
					if (m_data[i] == 0) {
						continue;
					}
					// Internal array is a matrix (M_NROWS, M_NCOLS) convert to this format
					int internal_row = i % M_NROWS;
					int internal_col = i / M_NROWS;

					// Each column in the internal array is a column of the band matrix 
					int external_row = internal_row + internal_col - 2*M_NCD;
					int external_col = internal_col;

					result[external_row] += rhs[external_col] * m_data[i];
				}
			}

			// Used to wrap a block within a BVPMatrix
			class Block
			{
			public: // Constructor
				Block(BVPMatrix& mat, BoundaryIndex b, uint nstr, uint nlyr):
					m_mat(mat)
				{
					m_row_offset = (b == 0) ? 0 : nstr / 2 * NSTOKES + (b - 1)*nstr*NSTOKES;
					m_col_offset = (b == 0) ? 0 : (b == nlyr) ? m_mat.M_NCOLS - nstr*NSTOKES : (b - 1) * nstr*NSTOKES;
				}
			public: // Interface
				inline double& operator()(StreamIndex i, SolutionIndex j) {
					return m_mat(m_row_offset + i, m_col_offset + j);
				}
				inline double operator()(StreamIndex i, SolutionIndex j) const {
					return m_mat(m_row_offset + i, m_col_offset + j);
				}
			private: // Members
				uint m_row_offset;
				uint m_col_offset;
				BVPMatrix& m_mat;
			};

			// Returns a block to the matrix with a top-left element (i0, j0).
			inline MatrixBlock<BVPMatrix> block(uint i0, uint j0, uint rows, uint cols) {
				return MatrixBlock<BVPMatrix>(*this, i0, j0, rows, cols);
			}

			// Get the block corresponding to the given boundary index.
			inline Block getBlock(BoundaryIndex b) {
				return Block(*this, b, M_NSTR, M_NLYR);
			}

			// Returns mutable pointer to underlying data.
			inline double* data() {
				return m_data.get();
			}

			// Returns the number of elements in data
			inline size_t size() {
				return M_NSTR * NSTOKES * M_NLYR;
			}

			// Returns the number of elements above/below the diagonal.
			inline lapack_int NCD() const {
				return M_NCD;
			}
			// Returns the number of columns (band storage).
			inline lapack_int N() const {
				return M_NCOLS;
			}
			// Returns the leading dimension of this matrix (band storage). 
			// Here it is equivalent to the number of rows (**in band 
			// storage**).
			inline lapack_int LD() const {
				return M_NROWS;
			}

		private: // Members
			// Number of elements above/below the diagonal
			const uint M_NCD;

			// The number of rows minus 1. Equivalent to M_NROWS-1
			const uint M_NRM1;

			// Number of rows in band-stored scheme
			const uint M_NROWS;

			// Number of columns in the original matrix (also in band-stored 
			// scheme)
			const uint M_NCOLS;

			const uint M_NSTR;
			const uint M_NLYR;
			std::unique_ptr<double[]> m_data;
		};

		int dgbsv_pentadiagonal(int N, int NRHS, double* AB, double* B, int LDB);

	}

    // Need dense blocks of the BVP matrix for the derivative calculations
    template <int NSTOKES, int CNSTR=-1>
	class BVPMatrixDenseBlock {
	public:
		BVPMatrixDenseBlock(BoundaryIndex b, uint nstr, uint nlyr) {
			m_b = b;
			if (b + 1 == nlyr) {
				m_upper_data.resize(nstr / 2 * NSTOKES, nstr * NSTOKES);
			}
			else
			{
				m_upper_data.resize(nstr * NSTOKES, 2 * nstr * NSTOKES);
			}
			m_upper_row_offset =  nstr / 2 * NSTOKES + (b)*nstr*NSTOKES;
			m_upper_col_offset = (b+1 == nlyr) ? (nlyr*nstr*NSTOKES) - nstr*NSTOKES : (b) * nstr*NSTOKES;
			if (m_b == 0) {
				m_data.resize(nstr * NSTOKES / 2, nstr * NSTOKES);
			}
			else {
				m_data.resize(nstr * NSTOKES, 2 * nstr * NSTOKES);
			}
			m_row_offset = (b == 0) ? 0 : nstr / 2 * NSTOKES + (b - 1)*nstr*NSTOKES;
			m_col_offset = (b == 0) ? 0 : (b == nlyr) ? (nlyr*nstr*NSTOKES) - nstr*NSTOKES : (b - 1) * nstr * NSTOKES;
		}
		Eigen::MatrixXd& upper() {
			return m_upper_data;
		}

		Eigen::MatrixXd& layer() {
			return m_data;
		}

		inline void setZero() {
			m_data.setZero();
			m_upper_data.setZero();
		}

		inline void multiplyVector(const Eigen::VectorXd& rhs, Eigen::VectorXd& result) {
			result.setZero();
			for (uint i = 0; i < m_data.rows(); ++i) {
				for (uint j = 0; j < m_data.cols(); ++j) {
					result[i + m_row_offset] += rhs[j + m_col_offset] * m_data(i, j);
				}
			}

			for (uint i = 0; i < m_upper_data.rows(); ++i) {
				for (uint j = 0; j < m_upper_data.cols(); ++j) {
					result[i + m_upper_row_offset] += rhs[j + m_upper_col_offset] * m_upper_data(i, j);
				}
			}
		}

        inline void assign_rhs_d_bvp(Eigen::MatrixXd& rhs, int colindex, const Eigen::VectorXd& d_b, const Eigen::VectorXd& b) {
            // First assign the column to d_b
            rhs(Eigen::all, colindex) = d_b;

            // Then subtract off the multiplication terms
            for (uint i = 0; i < m_data.rows(); ++i) {
                for (uint j = 0; j < m_data.cols(); ++j) {
                    rhs(i + m_row_offset, colindex) -= b[j + m_col_offset] * m_data(i, j);
                }
            }

            for (uint i = 0; i < m_upper_data.rows(); ++i) {
                for (uint j = 0; j < m_upper_data.cols(); ++j) {
                    rhs(i + m_upper_row_offset, colindex) -= b[j + m_upper_col_offset] * m_upper_data(i, j);
                }
            }
        }

	private:
		Eigen::MatrixXd m_upper_data;
		Eigen::MatrixXd m_data;
		BoundaryIndex m_b;
		uint m_row_offset;
		uint m_col_offset;

		uint m_upper_row_offset;
		uint m_upper_col_offset;
	};




}
