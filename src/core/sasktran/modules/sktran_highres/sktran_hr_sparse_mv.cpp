#include "include/sktran_hr_internals.h"



/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Sparse_mmv		2013-06-27*/
/** Performs sparse matrix by dense vector multiplication.  
 *  \param mvals The nonzero values of the sparse matrix, size=nnz
 *  \param mcolind The column each nonzero value resides in, size=nnz
 *  \param mrowptr Array of size=numrows+1, the first element is 0, then every
 *  element after is the running total of non-zero values in the current row
 *  \param v The dense vector to multiply
 *  \param result Vector to store the result in, must be preallocated
 *  \param numthreads The number of threads to use
 **/
/*---------------------------------------------------------------------------*/

void SKTRAN_HR_Sparse_mmv( std::vector<double>& mvals,
						   std::vector<int>&    mcolind,
						   std::vector<int>&    mrowptr,
						   std::vector<double>& v,
						   std::vector<double>& result,
						   int numthreads)
{
	int curthreads = omp_get_num_threads();
	omp_set_num_threads( numthreads );
	#pragma omp parallel for schedule(dynamic,10)
	for( int rowptridx = 0; rowptridx <(int)mrowptr.size()-1; rowptridx++ )
	{
		for( int colidx = mrowptr[rowptridx]; colidx < mrowptr[rowptridx+1]; colidx++ )
		{
			result[rowptridx] += mvals[colidx] * v[mcolind[colidx]];
		}
	}
	omp_set_num_threads( curthreads );
}
