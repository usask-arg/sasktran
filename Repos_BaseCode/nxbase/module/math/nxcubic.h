
/*-----------------------------------------------------------------------------
 *					class nxCubic									2005-3-17*/
/** \ingroup Math_fitting
 *	A class for solving cubic equations for real roots.
 *	Built so we can solve the cubic equation for thermistor curves.
 **/
/*---------------------------------------------------------------------------*/

class nxCubic
{
	public:
		nxBOOL		Solve                     ( double* coeffs, int numcoeffs, nx1dArray<double>* answers );
		nxBOOL		SolveAndChooseBestSolution( double* coeffs, int numcoeffs, double* bestanswer, double minx, double maxx );
};

