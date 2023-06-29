/*@doc DataAnalysis*/

/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/


 double  LSQ_ANALYTIC_CHI_SQUARE   ( const nxdblarr& y, const nxdblarr& yfit, int lim1, int lim2 );

 nxBOOL LSQ_ANALYTIC_MATRIX_INVERT(	const nxdblarr&			 BETA,
													const nx2dArray<double>& ALPHA,
													double					 lamda,
													nxdblarr*				 NEWA,
													const nxdblarr&			 OLDA);

 void   ALPHA_AND_BETA( nx2dArray<double>&	 derivs,
									   nxdblarr*	         BETA,
									   nx2dArray<double>*	 ALPHA,
									   const nxdblarr&		 Y,
									   const nxdblarr&		 yfit,
									   const nxdblarr&		 parameters,
									   int				     lim1,
									   int					 lim2);

/*<--------------------------------------------------------------------------
 *'					LsqAnalyticFit
 *	Performs a least squares fit to a 1-D function that whose value and
 *	gradients can be evaluated analytically (or numerically) given an array
 *	of parameters.
 *
 *	NEEDS UPGRADE WORK
 *
 *	Follows the algorithm given by Bevington.
 *
 *''PARAMETERS:
 *`		ANALYTIC&   function
 *		A class object that supports at least two functions {br}{br}
 *		{br}
 *		a) {code}void ANALYTIC::operator() ( const nxdblarr& parameters, nxdblarr* yfit, int lim1, int lim2 ){\code}{br}
 *		b) {code}void ANALYTIC::CHI_DERIVS ( const nxdblarr& parameters, nx2dArray<double>* derivs, int lim1, int lim2 ) {\code}{br}
 *		<br>
 *		The first function evaluates the function in the range lim1 to lim2 (yfit ranges from 0 to npts)
 *		The second function evaluates the derivatives for each parameter in the range  lim1 to lim2.
 *		derivs is of dimension (nparam, npts )
 *
 *`		const nxdblarr&  y
 *			The data to be fitted to. dimension (npts)
 *
 *`		nxdblarr*   parameters
 *			The array of parameters to be evaluated (These should be initialised
 *			to a close guess) dimension (nparam).  Upon exit these contain the
 *			best fit parameters.
 *
 *`		int lim1
 *			Starting limit for fitting.
 *
 *`		int lim2
 *			End limit for fitting
 *
 *`		int maxloop
 *			Maximum number of times round the system before failing
 *
 *`		nxdblarr* fittedline
 *			If not NULL then is set to the fitted lime in the range lim1 to lim2 only.
 *			Its f fimesnion (npts). All other points set to 0.
 *
 *
 *''	RETURNS:
 *	nxBOOL, nxTRUE if SUCCESS.
 *
 *''	HISTORY:
 * 2002-9-30
 *>------------------------------------------------------------------------*/


template <class ANALYTIC>
nxBOOL LsqAnalyticFit ( ANALYTIC&			function,
						const nxdblarr&		y,
						nxdblarr*			parameters,
						int					lim1,
						int					lim2,
						int					maxloop,
						nxdblarr*		    fittedline)
{

	int					npts      = y.N_Elements();
	int					nparams   = parameters->N_Elements();
	nx2dArray<double>   derivs( nparams, npts );
	nx2dArray<double>	alpha(nparams,nparams);
	nxdblarr		    beta(nparams);
	nxdblarr			olda(nparams);
	nxdblarr			newa(nparams);
	nxdblarr			yfit(npts);
	double				chisq1;
	double				chisqr = 0.0;
	double				change = 0.0;
	nxBOOL				reinvert= nxFALSE;
	int					L         = 0;
	double				lamda     = 0.1;
	double				chilimit  = 0.0001;
	nxBOOL				converged = nxFALSE;
	nxBOOL				ok;// = nxFALSE;

	yfit.SetTo(0.0);															// and clear it all to zero
	olda = *parameters;															// copy over the parameters
	newa = olda;																// and copy the old answers to the new answers
	function(olda, &yfit, lim1, lim2);											// get theoretical curve from
	chisq1 = LSQ_ANALYTIC_CHI_SQUARE( y, yfit, lim1, lim2 );					// get initial chisqr
	do
	{																			// get the alpha and beta arrays using
		function.CHI_DERIVS( olda, &derivs, lim1, lim2);						// Get the derivatives wrt variation in chi
		ALPHA_AND_BETA( derivs, &beta, &alpha, y, yfit, olda, lim1, lim2 );		// Get the gradient matrices around here.
		do																		// Keep increasing lamda
		{																		// invert ALPHA and BETA
			ok = LSQ_ANALYTIC_MATRIX_INVERT(beta, alpha, lamda, &newa, olda );	// solve the with this lamda use new answer to get new fitted
			if (!ok) nxLog::Record(NXLOG_WARNING, "::LsqAnalyticFit, Error inverting derivative array ");
			if (ok)																// if it inverted ok
			{																	// then
				function( newa, &yfit, lim1, lim2  );							// calculate the new function
				chisqr   = LSQ_ANALYTIC_CHI_SQUARE( y, yfit, lim1, lim2 );		// Get chi-square
				change   = (chisqr - chisq1) / chisq1;							// check the change in chi-squared
				reinvert = (change > 0.0);										// Shall I reinvert ALPHA
				if (reinvert) lamda *= 10;										// do I keep looking around chi-squared space.
				L++;															// keep track of the loops.
				ok = L <= maxloop;												// and abort procesisng if it overflows.
			}																	// otherwise matrix invert Failed so abort procesing
		} while (reinvert && ok);												// keep looping until ch-squared decreases.
		if (ok)																	// if we successfully found a lower range
		{																		// then
			converged = (fabs(change) < chilimit);								// check if converged
			olda    = newa;														// pass over new values
			lamda  /= 10.0;														// Decrease lamda
			chisq1  = chisqr;													// reset initial chisqr
		}																		//
	} while (!converged && ok);													// Loop around until convergence. So
																				// exit when converged
	*parameters = newa;
	if (fittedline != NULL) *fittedline = yfit;
	return ok;
}
