#include "nxbase_math.h"

/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/

//---------------------------------------------------------------------------
//+
//NAME:
//						nxArray::Max
//	Returns the minimum value in the array.  If the array is undefined then
//	returns 0.
//-
//---------------------------------------------------------------------------

double nxdblarr::Max( int i1, int i2) const
{
	double	maxval;
	int		i;

   	i1   = LimitBoundOfIndex( i1 );
	i2   = LimitBoundOfIndex( i2 );
	if (size() > 0)
	{
		maxval  = At(i1);
		for ( i= i1+1; i <= i2; i++) maxval = (maxval > At(i))? maxval: At(i);
	}
	else
	{
		maxval = 0.0;
	}
	return maxval;
}

double nxdblarr::Min( int i1, int i2) const
{
	double	maxval;
	int		i;

   	i1   = LimitBoundOfIndex( i1 );
	i2   = LimitBoundOfIndex( i2 );
	if (size() > 0)
	{
		maxval  = At(i1);
		for ( i= i1+1; i <= i2; i++) maxval = (maxval < At(i))? maxval: At(i);
	}
	else
	{
		maxval = 0.0;
	}
	return maxval;
}

//---------------------------------------------------------------------------
//						nxArray::MaxValueBelowThreshold
//	Find the maximum value in the array which is less than the
//  specified threshold.  Useful for quickly sorting the order
//  of arrays.
//
//	Returns the minimum value in the arary within the specified value.
//	maxbin contains the index of the maximum value (below the thershold)
//	if maxbin is less than 0 then no maximum was found and the value 0
//	is returned.
//---------------------------------------------------------------------------

double nxdblarr::MaxValueBelowThreshold(  double threshold, int lim1, int lim2, int* maxbin) 
{ 
	double		value;
	double		maxval = 0.0;
	int		i;
	
   	lim1   = LimitBoundOfIndex( lim1 );					// Make sure Lim1
	lim2   = LimitBoundOfIndex( lim2 );					// and lim2 are within array bounds
	*maxbin = -1;										// default to no maximum found

	if ( N_Elements() > 0)								// if the array is defined
	{													// then
		for (i = lim1; i <= lim2; i++)					// find the first value
		{												// whose value
			maxval = At(i);								// is below the threshold
			if (maxval < threshold)						// have we found it yet
			{											// yup
				*maxbin = i;							// so start the maxbin at this point 
				break;									// and break out of this loop
			}
		}

		if (*maxbin >= 0)										// if we found a point less  tha threshold
		{														// then scan remainder of array for max value
			for (i = *maxbin+1; i <= lim2; i++)
			{
				value = At(i);
				if ((value > maxval) && ( value < threshold ))
				{
					maxval = value;
					*maxbin = i;
				}
			}
		}
	}
	if (*maxbin == -1) maxval = 0.0; //(T)0;
	return maxval;
} 

//---------------------------------------------------------------------------
//+
//NAME:
//						nxdblarr::QuadraticCentredOnPeak
//	Least squares a parabola (quadratic) to an array between limits lim1
//	and lim2.  The code was developed for FPI fringes. The basic idea is
//	to adjust the limits until the peak is centred in the middle of the
//	limits.  This attempts to remove slight biases due to asymetric limits
//	is the least squares analysis.
//-
//---------------------------------------------------------------------------

nxBOOL nxdblarr::QuadraticCentredOnPeak( int     lim1, int lim2,
										 double* peak,
										 double* peakerr,
										 double* intensity,
										 double* A,    double* B,       double* C)
{
	const int	maxtry  = 5;										// maximum nuber of fitting tries 
	double		window;
	int			ntry;												// counts Number of fit tries 
	double		centre;												// Limits centre position 
	nxBOOL		centred; // = nxFALSE;
	nxBOOL		ok;
	double		wind2;												// half the fitting window

	lim1    = LimitBoundOfIndex( lim1 );							// make sure the limits are within bounds
	lim2    = LimitBoundOfIndex( lim2 );							// for both upper and lower limits
	window  = (lim2 - lim1 + 1);									// window for fitting
	wind2   = window / 2.0;											// half the fitting window
	centred = nxFALSE;												// Clear centred on peak flag
	ntry    = 0;													// count number of tries

	while (!centred && (ntry < maxtry))								// iterate upto maxtry times
	{
		ok = LeastSquareQuadratic( A, B, C, peak, intensity, peakerr, lim1, lim2);		// do the fit
		if ( ok && (*peak > lim1) && (*peak < lim2) && (*A < 0.0))	// is the peak within fitting limits.
		{															// yes so check peakk centering of limits.
			centre = lim2 - wind2;									// get centre of limits
			if ((fabs(*peak - centre) < 1.9))						// is peak centred about limits
			{														// Yup 
				centred = nxTRUE;									// Set Centred on peak flag 
				(*peak) += 0.5;										// bin 0 is really X= +0.5  
			}
			else													// peak not centred within limits
			{														// so
				lim1 = int(*peak - wind2 + 0.5);						// setup new lower limit
				lim2 = int(*peak + wind2 + 0.5);						// setup new upper limit 
				lim1 = LimitBoundOfIndex( lim1 );
				lim2 = LimitBoundOfIndex( lim2 );
				ntry++;												// increment No of tries
				*peak = 0;											// place duff No in peak value
			} 
		}
		else
		{
			ntry = maxtry + 1;										// drops out of loop
			*peak = 0;												// set a bad value
		}
	} 
	return centred;
} 

//---------------------------------------------------------------------------
//
//						nxdblarr:LeastSquaresQuadratic
//	Routine to fit a parabola to integer spaced points. The routine	 
//	uses the mid-range of the points as the zero position. The sum
//	sets SUM X, SUM X**3 and SUM Y to zero. The answer is still the
//	same, but the computation is more accurate.
//
//	Parameters
//	----------
//	y(i)	The array containing the data.
//	a	Returned value of a in quadratic eqn.
//	b	Returned value of b in quadratic eqn.
//	c	Returned value of c in quadratic eqn.
//	pk	Peak position along x axis.
//	int	Height of parabola at turning point.
//	pkerr	absolute error on peak position.
//	lim1	lower bin number for fit
//	lim2	upper bin number for fit
//
//	If an error occurs then all parameters, except A, are set to
//	zero. 'A' is set to +1.
//
//---------------------------------------------------------------------------

nxBOOL nxdblarr::LeastSquareQuadratic(double *A, double *B, double *C, double *PK, double *pkHeight, double *PKERR, int LIM1, int LIM2)
{
	double       AD, BD, CD;
	double       DEL;
	double       YMAX;
	double       YX1, YX2;
	double       X2, X4, FN;
	double       ZEROX, ZEROY;
	double       X, Y;
	int          I;
	nxBOOL		ok = nxFALSE;


	*A = 1.0;											/* Error value Coefficient of X**2 */ 
	*B = 0.0;											/* Coefficient of X */ 
	*C = 0.0;											/* Constant */ 
	*PK = 0.0;											/* Peak/Trough position */ 
	*pkHeight = 0.0;									/* Peak trough Value */ 
	*PKERR    = 2000.0;									/* peak error  */ 

	YMAX = 0.0;											/* Max y value in range */ 
	YX1 = 0.0;											/* sum of Y*X */ 
	YX2 = 0.0;											/* sum of Y*X**2 */ 

	X2 = 0.0;											/* sum of X**2 */ 
	X4 = 0.0;											/* sum of X**4 */ 
	FN = LIM2 - LIM1 + 1.0;								/* Number of data points */ 

	ZEROY = 0.0;										/* new Y origin */ 
	if ((FN > 2.0))										/* only fit if enough points */ 
	{ 
		ZEROX = LIM1 + 0.5 * (LIM2 - LIM1);				/* get the new X origin */ 
		for (I = LIM1; I <= LIM2; I++) 
		{ 
			ZEROY += (double)At(I);						/* get average Y coord */ 
		} 
		ZEROY = ZEROY / FN;								/* get New Y origin */ 
		for (I = LIM1; I <= LIM2; I++) 
		{ 
			double yval;

			yval = (double)At(I);
			if (yval > YMAX) YMAX = yval;				/* Find max Y value */ 
			Y   = (yval - ZEROY);						/* New Y coord */ 
			X   = I - ZEROX;							/* new X coord */ 
			YX1 += Y*X;									/* Sum of Y*X */ 
			X   = nxmath::sqr(X);								/* X**2 */ 
			X2  += X;									/* Sum of X**2 */ 
			YX2 += Y*X;									/* sum of Y*X**2 */ 
			X    = nxmath::sqr(X);								/* X**4 */ 
			X4  += X;									/* sum of X**4 */ 
		} 
		DEL = (FN * X4 - nxmath::sqr(X2));						/* partial result */ 
		AD = (FN * YX2) / DEL;							/* A coeff in new origin space */ 
		BD = YX1 / X2;									/* B coeff in new origin space */ 
		CD = -AD* X2;									/* C coeff in new origin space */ 
		if (AD != 0.0)									/* if all data = 0 */ 
		{ 
			*PKERR = (FN / (DEL * nxmath::sqr(AD)) + 1.0 / X2) / (4.0 * nxmath::sqr(AD)) * YMAX; 
			*PKERR = sqrt(*PKERR);						/* get the peak error */ 
		} 
		*C = CD + AD * nxmath::sqr(ZEROX) - BD * ZEROX + ZEROY;	/* transform back to old space */ 
		*B = BD - 2.0 * ZEROX * AD;						/* actual B coeff */ 
		*A = AD;										/* actual A coeff */ 
		ok = nxTRUE;
		if (fabs(*A) > 1.0E-6)							/* find peak if A <> 0.0 */
		{ 
			*PK = -(*B)/(2.0*(*A));						/* Peak/Trough position */ 
			*pkHeight = *C - nxmath::sqr(*B) / (4.0 * (*A));		/* intensity at peak */ 
		} 
	} 
	return ok;
}
//---------------------------------------------------------------------------
//						nxdblarr::AdjustParabolaFittingWindow
//
//   Routine adjust the limits takes the maximum and minimum bins (maxbin and MINBIN) and 
//   calculates a suitable range of limits given by a peak half-height
//   criterion. If the peak is near the edge of the interferogram then
//   it is possible that either one limit or the other does fall below
//   the half-height criterion. In this case goodlimits is set to false
//   and the calling program should seek a new set of limits that will
//   satisfy the half-height criterion.
//---------------------------------------------------------------------------
 
nxBOOL nxdblarr::AdjustParabolaFittingWindow( int maxbin, double minvalue, int *lim1, int *lim2, int maxrad) 
{ 
	double      halfht;
	int         l1, l2;															// Buffers limit values 
	nxBOOL		goodlimits;

	maxbin = LimitBoundOfIndex(maxbin);
	halfht = (double)minvalue + 0.5*((double)At(maxbin) - (double)minvalue);	// get the half height 
 //	l2     = *lim2;															// starting point for upper limit
 //	l1     = *lim1;															// starting point for lower limit
	maxrad = LimitBoundOfIndex( maxrad );

//	while ((l2 >  maxbin) && ( (double)At(l2) < halfht)) l2--;					// Get the upper edge for the parabola fit.
//	while ((l1 <  maxbin) && ( (double)At(l1) < halfht)) l1++;					// Get the lower edge for the parabola fit

	l2 = maxbin;
	l1 = maxbin;
	while ((l2 < maxrad)  && ( (double)At(l2) > halfht)) l2++;					// Get the upper edge for the parabola fit.
	while ((l1 >= 0)      && ( (double)At(l1) > halfht)) l1--;					// Get the lower edge for the parabola fit
	goodlimits = (l1 >= 0) && (l2 < maxrad );
	if (goodlimits)																// if limits are good
	{																			// then pass results 
		*lim1 = l1;																// back to caller 
		*lim2 = l2;																// so he can use them
	}																			// for a parabola fit
	else
	{
	//	TRACE("Failed to get goodlimits, maxbin %d, max %lf, min %lf, halfht %lf, lim1 %d, lim2 %d \n" , (int)maxbin, (double)At(maxbin), (double)minvalue, (double)halfht, (int)*lim1, (int)*lim2);
		*lim1 = 1;
		*lim2 = 0;
	}
	return goodlimits;
}
//--------------------------------------------------------------------------- 
//						nxdblarr::ParabolaPeakFitter
//
//	Routine to try and locate the best peak within the
//	range lim1 to lim2. The best peak is defined as the
//	heighest peak in the range as long as parabola limits
//	can extend half-way down the peak. This is okay as long
//  as the peak is not near the edge of the window.
//
//	If limits cannot be fitted around the first maximum
//	then the second maximum is tried and so on until a good
//	fit limit is found or until data runs out.
//	
//	Once a range of limits are found then a parabola is
//	fitted to the data. The parabola fit may also adjust
//	the fit limits to centre them on the determined peak
//	position.
//
//	Upon exit lim1 and lim2 contain the limits used for
//	the final parabola fit. If an error occurs then lim2
//	is less than lim1.
//	
//---------------------------------------------------------------------------

nxBOOL nxdblarr::ParabolaPeakFitter( double *peak, double *peakerr, double *pkHeight, int* Lim1, int* Lim2, int maxrad) 
{ 
	int         maxbin;					/* maximum bin from sorting arrays */ 
	double		maxvalue;				/* Last maximum value of scan beforehand */ 
	double		minvalue;				/* Current minimum value in the whole array */ 
	double      A, B, C;				/* A,B,C of quadratic equation. Dummies for QUADRO */ 
	nxBOOL      goodlimits;				/* checks if we have good limits */ 
	int			lim1 = *Lim1;
	int			lim2 = *Lim2;
	int         ntry = 0;

	maxvalue = Max(lim1, lim2) + 1;
	minvalue = Min(lim1,  lim2);
	do
	{
		lim1 = *Lim1;
		lim2 = *Lim2;
		maxvalue   = MaxValueBelowThreshold( maxvalue, lim1, lim2, &maxbin       );			/* passed in maxvalue */ 
		goodlimits = AdjustParabolaFittingWindow  ( maxbin, minvalue, &lim1, &lim2, maxrad);		/* about maximum bin */ 
		ntry++;
//		TRACE("parabola limits = %d %d max at %d\n", lim1, lim2, maxbin);
	} while (!goodlimits && (maxvalue > minvalue) && (ntry < (*Lim2 - *Lim1)));																/* keep going until good */ 

	if (goodlimits)
	{
		goodlimits = QuadraticCentredOnPeak( lim1, lim2, peak, peakerr, pkHeight, &A, &B, &C);
//		TRACE( "parabola centre, lim1, lim2, peak, %d %d %lf \n", lim1, lim2, (double)*peak );
		*Lim1 = lim1;
		*Lim2 = lim2;
	}
	else
	{
		*Lim1 = 1;
		*Lim2 = 0;
	//	TRACE("No good limits found in range");
	}
	return goodlimits;
}
