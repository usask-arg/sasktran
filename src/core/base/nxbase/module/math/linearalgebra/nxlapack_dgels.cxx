#include "nxbase_linearalgebra.h"
#include "sklapack.h"


/*-----------------------------------------------------------------------------
 *					LapackMatrix::DGELS		2007-7-24*/
/** Solves an over determined system using least squares approach. **/
/*---------------------------------------------------------------------------*/

bool LapackMatrix::Dgels( const nxArrayLinear<double>& user_y, nxArrayLinear<double>* x )
{
	blas::f77integer	M      = (int)NumRows();
	blas::f77integer	N      = (int)NumColumns();
	size_t				n      = N;
	blas::f77integer	LDA    = M;
	blas::f77integer	INFO   = -1;
	blas::f77char		TRANS  = 'N';
	blas::f77integer	NRHS   = 1;
	blas::f77integer	LDB    = M;
	blas::f77integer	LWORK  = 16*M;
	bool				ok;
    nx1dArray<double>	work;
	nx1dArray<double>	y;
	LapackMatrix		column;
	size_t				i;

	ok = (M == (int)user_y.size() );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "LapackMatrix::DGELS, This matrix and the y array are not the same size, I have not made the fit");
	}
	else
	{
		ok = ok && work.SetSize( LWORK );				// Allocate a large amount of memory in anticipation of calling DGELS
		ok = ok && x->SetSize(1, &n);
		ok = ok && y.DeepCopy(user_y);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "LapackMatrix::DGELS, Error allocating buffers for analysis");
		}
		else
		{
			LAPACK(dgels) ( &TRANS, &M, &N, &NRHS, UnsafeArrayBasePtr(), &LDA, y.UnsafeArrayBasePtr(), &LDB, work.UnsafeArrayBasePtr(), &LWORK, &INFO );
			ok = (INFO ==0);
			if (ok)
			{
				for (i = 0; i < n; i++)
				{
					x->At(&i) = y.At(i);
				}
			}
		}
	}
	if (!ok) x->erase();
	return ok;
}
