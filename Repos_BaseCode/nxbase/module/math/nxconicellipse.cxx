#include "nxbase_math.h"

#if defined(NX_WINDOWS)
#pragma message ("nxconicellipse.cxx still need major edits to properly compile *** DO IT***. I have kluged in defines to ignore it right now.") 
#endif

#if 0
#include "nxbase_linearalgebra.h"

using namespace nxmath;
/*---------------------------------------------------------------------------
 *						nxConicEllipse::nxConicEllipse
 *-------------------------------------------------------------------------*/

nxConicEllipse::nxConicEllipse()
{
	m_fitrotated = nxTRUE;
	m_A     = 0;		// A coeff from Generalised equation F = Ax^2 + Bxy +Cy^2 +Dx + Ey +F = 0	 x = cos(@) y = sin(@)
	m_B     = 0;		// B coeff from Generalised equation
	m_C     = 0;		// C coeff from Generalised equation
	m_D     = 0;		// D coeff from Generalised equation
	m_E     = 0;		// E coeff from Generalised equation
	m_F     = 0;		// F coeff from Generalised equation
	m_x0    = 0;		// X Origin of ellipse
	m_y0    = 0;		// Y Origin of ellipse
	m_theta = 0;		// Inclination of ellipse "m_a" axis to coordinate X axis
	m_a     = 0;		// Length of semi-major/minor axis closest to coordinate X axis.
	m_b		= 0;		// Length of semi/major/minor axis closest to coordinate Y axis
}


/*---------------------------------------------------------------------------
 *						nxConicEllipse::MapParamToSpace
 *-------------------------------------------------------------------------*/

nxBOOL nxConicEllipse::MapParamToSpace(double xparam, double  yparam, double* xspace, double* yspace )
{
	double xprime, yprime;
	double costheta = cos(m_theta);
	double sintheta = sin(m_theta);

	xprime =   xparam*costheta - yparam*sintheta;
	yprime =   xparam*sintheta + yparam*costheta;
	*xspace = xprime + m_x0;
	*yspace = yprime + m_y0;
	return nxTRUE;
}

/*---------------------------------------------------------------------------
 *						nxConicEllipse::MapSpaceToParam
 * The X parameter space is given by x = a.cos(angle)
 * The Y parameter space is given by y = b.sin(angle)
 *-------------------------------------------------------------------------*/

nxBOOL nxConicEllipse::MapSpaceToParam(double xspace, double  yspace, double* xparam, double* yparam )
{
	double xprime, yprime;
	double costheta = cos(m_theta);
	double sintheta = sin(m_theta);

	xprime = xspace - m_x0;
	yprime = yspace - m_y0;

	*xparam =  xprime*costheta + yprime*sintheta;
	*yparam = -xprime*sintheta + yprime*costheta;
	return nxTRUE;
}

/*---------------------------------------------------------------------------
 *						nxConicEllipse::MapSpaceToParamAngle
 * The X parameter space is given by x = a.cos(angle)
 * The Y parameter space is given by y = b.sin(angle)
 *-------------------------------------------------------------------------*/

nxBOOL nxConicEllipse::MapSpaceToParamAngle(double xspace, double  yspace, double* angle)
{
	double xprime, yprime;

	MapSpaceToParam(xspace, yspace, &xprime, &yprime);
	*angle = atan2d( yprime/m_b, xprime/m_a);
	return nxTRUE;
}


/*---------------------------------------------------------------------------
 *						nxConicEllipse::Clear
 *-------------------------------------------------------------------------*/

void nxConicEllipse::Clear()
{
	m_A = 0.0;
	m_B = 0.0;
	m_C = 0.0;
	m_D = 0.0;
	m_E = 0.0;
	m_F = 0.0;
}

/*---------------------------------------------------------------------------
 *						nxConicEllipse::FromParameters
 *-------------------------------------------------------------------------*/

void nxConicEllipse::FromParameters( double a, double b, double x0, double y0,  double rotationangle )
{
	double cosa;
	double sina;
/*
	double			Q_Storage     [3*3];
	double			TEMP_Storage  [3*3];
	double			QPRIME_Storage[3*3];
	double			R_Storage     [3*3];
	double			T_Storage     [3];
	double			S_Storage     [3];

	LapackMatrix 	   Q(3,3,Q_Storage);
	LapackMatrix	TEMP(3,3,TEMP_Storage);
	LapackMatrix  QPRIME(3,3,QPRIME_Storage);
	LapackMatrix       R(3,3,R_Storage);
	LapackMatrix	   T(3,1,T_Storage);
	LapackMatrix	   S(3,1,S_Storage);
*/
	
	cosa     = cosd(rotationangle);
	sina     = sind(rotationangle);
	m_A      = sqr( cosa/a) + sqr(sina/b);
	m_C      = sqr( sina/a) + sqr(cosa/b);
	m_B      = 2*cosa*sina*(sqr(1.0/a) - sqr(1.0/b));
	m_D      = -m_B*y0 - 2*m_A*x0;
	m_E      = -m_B*x0 - 2*m_C*y0;
	m_F      = m_A*sqr(x0) + m_C*sqr(y0) + m_B*x0*y0 - 1;

	UpdateMmatrix();
}


/*---------------------------------------------------------------------------
 *						nxConicEllipse::UpdateMmatrix
 *-------------------------------------------------------------------------*/

void nxConicEllipse::UpdateMmatrix()
{
	double			COSTHETA;
	double			SINTHETA;
	double			scale;

	scale    = 1.0/(4.0*m_A*m_C -m_B*m_B);
	m_x0     = ( -2.0*m_D*m_C  + m_E * m_B) * scale;			// Solve for X0 and Y0
	m_y0     = ( -2.0*m_A*m_E  + m_D * m_B) * scale;			// Which is the origin of the ellips ein user coordinates
	double AMC = m_A-m_C;
	if (AMC == 0.0) m_theta = 0.0;			// Circle solution
	else									// otherwise
	{
		m_theta = atan( m_B/AMC)/2.0;		// Get the ellipse solution
	}

	COSTHETA = cos(m_theta);
	SINTHETA = sin(m_theta);

	double cost2  = sqr(COSTHETA);
	double sint2  = sqr(SINTHETA);

	m_a = (cost2-sint2)/( m_A*cost2 - m_C*sint2);
	m_b = (sint2-cost2)/( m_A*sint2 - m_C*cost2);
	m_a = sqrt(m_a);
	m_b = sqrt(m_b);
	double f;


	if (m_a >= m_b )
	{ 	
		f = (m_a-m_b)/m_a;
		SetGeoid( m_a, f );
	}
	else
	{
		f = (m_b-m_a)/m_b;
		SetGeoid( m_b, f );
	}
}

/*---------------------------------------------------------------------------
 *						nxConicEllipse::GradientAsAngle
 *-------------------------------------------------------------------------*/

double nxConicEllipse::GradientAsAngle( double x, double y )
{
	double dy = -m_B*y - 2.0*m_A*x - m_D;
	double dx =  m_B*x + 2.0*m_C*y + m_E;
	double theta;

	if (dx == 0) theta = sign(dy)*90.0;
	else         theta = atan( dy/dx )*ONE_RADIAN;
	return theta;
}

/*---------------------------------------------------------------------------
 *						nxConicEllipse::FromDataPoints
 *-------------------------------------------------------------------------*/

nxBOOL nxConicEllipse::FromDataPoints( double*xp, double* yp, int npts )
{
	nxBOOL	ok;

	if (m_fitrotated) ok = FitRotatedEllipse( xp, yp, npts);
	else              ok = FitNonRotatedEllipse( xp, yp, npts );
	return ok;
}

/*---------------------------------------------------------------------------
 *						nxConicEllipse::FromFittedData
 *-------------------------------------------------------------------------*/

nxBOOL nxConicEllipse::FitRotatedEllipse( double*xp, double* yp, int npts )
{
	double						A_Storage[5*5];				// Use local storage for the matrices.
	double						B_Storage[5*1];				// as it helps on the dynamic memory allocation
	double						X_Storage[6*1];				// X_Storage  must be 6 even though X is dimesnionsed as 5.
	double						C_Storage[6*6];
	double						S_Storage[6*6];
	double						E_Storage[6*6];

	LapackFactorLU				A(5,5, A_Storage);			// Matrices used
	LapackMatrix				B(5,1, B_Storage);
	LapackMatrix				X(5,1, X_Storage);
	LapackMatrix				D;							// This memory allocation is dynamic
	LapackMatrix				C(6,6,C_Storage);
	LapackGeneralizedSymEigen	S(6,6,S_Storage);
	LapackMatrix				E(6,6,E_Storage);
	nxArray<double>				U;
	nxArray<double>				eigenvalues;
	int							i;
	double						lambda;
	nxBOOL						ok;
	double						maxval;
	int							maxindex = -1;

	ok = (npts >= 6);
	if (!ok)
	{
		nxLog::Verbose(NXLOG_WARNING, "nxConicEllipse::FromDataPoints, cannot fit conic to %d points.  Need at least 6", (int)npts );
	}
	else
	{
		D.Create( npts, 6, 0.0 );							// Create the Design matrix
		C(3,1) =  2.0;										// Get the constraint matrix
		C(1,3) =  2.0;
		C(2,2) = -1.0;

		nxArray<double>&	x2  = D.Column(1);				// Get the columns from the Design matrix
		nxArray<double>&	xy  = D.Column(2);
		nxArray<double>&	y2  = D.Column(3);
		nxArray<double>&	x   = D.Column(4);
		nxArray<double>&	y   = D.Column(5);
		nxArray<double>&	one = D.Column(6);
		one.SetTo(1);

		for ( i = 0; i < npts; i++)							// Fill the design matrix with values
		{
			double px = xp[i];
			double py = yp[i];
			x [i] = px;
			y [i] = py;
			xy[i] = px*py;
			y2[i] = py*py;
			x2[i] = px*px;
		}

		LapackBLAS::GEMM(1.0, D, 'T', D, ' ', &S, 0.0);		// Calculate the scatter matrix, S = D^T*D
		ok = S.CholeskyFactorize();							// Factorize using SPD matrix
		if (ok)												// It will fail if the data is a perfect conic
		{
			S.SolveEigenvector( C, &E, &eigenvalues );			// Solve for the generalized eigenvector 
			maxval   = 1.0E-20;
			maxindex = -1;
			lambda = 0.0;										// Choose the elliptical solution (there will be one)
			for (i = 0; i < eigenvalues.N_Elements(); i++)		// Now find the most positive solution
			{													// There should be  a few solutions close to zero
				if (eigenvalues[i] > maxval)					// but just one solution that is clearly 
				{												// non zero.
					maxval   = eigenvalues[i];
					maxindex = i+1;
					lambda    = eigenvalues[i];
				}
			}
			ok = (maxindex > 0);
			if (ok) U.AssignToABuffer( (double *)E.Column(maxindex), 6);
		}
		else													// We could do a least squares fit
		{														// This usually means its a perfect ellipse. so just solve for that!
			nxLog::Verbose(NXLOG_WARNING, "nxConicEllipse::FromDataPoints, It looks like a perfect conic. We shall just fit the first 6 points");
			
			int		dn = (npts-1)/4;					// Select 5 points evenly spaced over the data point interval
			int		idx;

			for (i=1; i <= 5; i++)						// Set up the 5 x 5 array
			{
				idx = (i-1)*dn;							
				A.At(i,1) = x2[idx];
				A.At(i,2) = xy[idx];
				A.At(i,3) = y2[idx];
				A.At(i,4) = x[idx];
				A.At(i,5) = y[idx];
				B.At(i,1) = 1;
			}	
			ok = A.LUFactorize();							// Factorize the matrix
			ok = ok && A.LinearSolve( B, &X );				// Solve the 5 equations simultaneously
			if (ok)											// and if that was ok
			{												// then
				X_Storage[6-1] = -1.0;						// Copy the value of "F" to the X storage array	
				U.AssignToABuffer( X_Storage, 6 );			// Now assign the answer array to the X array.
			}
		}

		if (ok)												// If one fit or the other was good then renormalize the conic
		{													// 	so
			double a,b,c,d,e,f,scale, factor;				// internal variable

			a     = U.At(0);								// copy the coefficients
			b     = U.At(1);
			c     = U.At(2);
			d     = U.At(3);
			e     = U.At(4);
			f     = U.At(5);

			factor = (4.0*a*c - b*b);						// get the "ellipse" factor. It mus be positive
			ok = (factor > 0);
			if (!ok)
			{
				nxLog::Verbose(NXLOG_WARNING, "nxConicEllipse::FromDataPoints, The solution is not an ellipse because 4ac-b*b  is less than or equal to 0. (%e)",(double)factor);
			}
			else
			{
				scale = 1.0/factor;									// Calculate a normalizing factor
				double x0 = ( -2.0*d*c  + e * b) * scale;			// Solve for X0 and Y0
				double y0 = ( -2.0*a*e  + d * b) * scale;			// Which is the origin of the ellips ein user coordinates
																	// Now use the fact that F(x0,y0) = -1 for an ellipse
				f += a*x0*x0 + b*x0*y0 + c*y0*y0 + d*x0 + e*y0;		// This can be used to get the proper scaling
				scale = -1.0/f;										// get the scaling
				U  *= scale;
				m_A = U.At(0);
				m_B = U.At(1);
				m_C = U.At(2);
				m_D = U.At(3);
				m_E = U.At(4);
				m_F = U.At(5);
				UpdateMmatrix();
			}
		}
	}
	if (!ok) Clear();
	return ok;
}


/*---------------------------------------------------------------------------
 *						nxConicEllipse::FromFittedData
 *-------------------------------------------------------------------------*/

nxBOOL nxConicEllipse::FitNonRotatedEllipse( double*xp, double* yp, int npts )
{
	double						A_Storage[4*4];				// Use local storage for the matrices.
	double						B_Storage[4*1];				// as it helps on the dynamic memory allocation
	double						X_Storage[4*1];				// X_Storage  must be 6 even though X is dimesnionsed as 5.

	double						C_Storage[5*5];
	double						S_Storage[5*5];
	double						E_Storage[5*5];

	LapackFactorLU				A(4,4, A_Storage);			// Matrices used
	LapackMatrix				B(4,1, B_Storage);
	LapackMatrix				X(4,1, X_Storage);

	LapackMatrix				D;							// This memory allocation is dynamic
	LapackMatrix				C(5,5,C_Storage);
	LapackGeneralizedSymEigen	S(5,5,S_Storage);
	LapackMatrix				E(5,5,E_Storage);
	nxArray<double>				U;
	nxArray<double>				eigenvalues;
	int							i;
	double						lambda;
	nxBOOL						ok;
	double						maxval;
	int							maxindex = -1;

	ok = (npts >= 5);
	if (!ok)
	{
		nxLog::Verbose(NXLOG_WARNING, "nxConicEllipse::FromDataPoints, cannot fit conic to %d points.  Need at least 6", (int)npts );
	}
	else
	{
		D.Create( npts, 5, 0.0 );							// Create the Design matrix
		C(2,1) =  2.0;										// Get the constraint matrix
		C(1,2) =  2.0;

		nxArray<double>&	x2  = D.Column(1);				// Get the columns from the Design matrix
		nxArray<double>&	y2  = D.Column(2);
		nxArray<double>&	x   = D.Column(3);
		nxArray<double>&	y   = D.Column(4);
		nxArray<double>&	one = D.Column(5);
		one.SetTo(1);

		for ( i = 0; i < npts; i++)							// Fill the design matrix with values
		{
			double px = xp[i];
			double py = yp[i];
			x [i] = px;
			y [i] = py;
			y2[i] = py*py;
			x2[i] = px*px;
		}

		LapackBLAS::GEMM(1.0, D, 'T', D, ' ', &S, 0.0);		// Calculate the scatter matrix, S = D^T*D
		ok = S.CholeskyFactorize();							// Factorize using SPD matrix
		if (ok)												// It will fail if the data is a perfect conic
		{
			S.SolveEigenvector( C, &E, &eigenvalues );			// Solve for the generalized eigenvector 
			maxval   = 1.0E-20;
			maxindex = -1;
			lambda = 0.0;										// Choose the elliptical solution (there will be one)
			for (i = 0; i < eigenvalues.N_Elements(); i++)		// Now find the most positive solution
			{													// There should be  a few solutions close to zero
				if (eigenvalues[i] > maxval)					// but just one solution that is clearly 
				{												// non zero.
					maxval   = eigenvalues[i];
					maxindex = i+1;
					lambda    = eigenvalues[i];
				}
			}
			ok = (maxindex > 0);
			if (ok) U.AssignToABuffer( (double *)E.Column(maxindex), 5);
		}
		else													// We could do a least squares fit
		{														// This usually means its a perfect ellipse. so just solve for that!
			nxLog::Verbose(NXLOG_WARNING, "nxConicEllipse::FromDataPoints, It looks like a perfect conic. We shall just fit the first 6 points");
			
			int		dn = (npts-1)/4;					// Select 5 points evenly spaced over the data point interval
			int		idx;

			for (i=1; i <= 4; i++)						// Set up the 5 x 5 array
			{
				idx = (i-1)*dn;							
				A.At(i,1) = x2[idx];
				A.At(i,2) = y2[idx];
				A.At(i,3) = x[idx];
				A.At(i,4) = y[idx];
				B.At(i,1) = 1;
			}	
			ok = A.LUFactorize();							// Factorize the matrix
			ok = ok && A.LinearSolve( B, &X );				// Solve the 5 equations simultaneously
			if (ok)											// and if that was ok
			{												// then
				X_Storage[5-1] = -1.0;						// Copy the value of "F" to the X storage array	
				U.AssignToABuffer( X_Storage, 5 );			// Now assign the answer array to the X array.
			}
		}

		if (ok)												// If one fit or the other was good then renormalize the conic
		{													// 	so
			double a,b,c,d,e,f,scale, factor;				// internal variable

			a     = U.At(0);								// copy the coefficients
			b     = 0;
			c     = U.At(1);
			d     = U.At(2);
			e     = U.At(3);
			f     = U.At(4);

			factor = (4.0*a*c);						// get the "ellipse" factor. It mus be positive
			ok = (factor > 0);
			if (!ok)
			{
				nxLog::Verbose(NXLOG_WARNING, "nxConicEllipse::FromDataPoints, The solution is not an ellipse because 4ac-b*b  is less than or equal to 0. (%e)",(double)factor);
			}
			else
			{
				scale = 1.0/factor;									// Calculate a normalizing factor
				double x0 = ( -2.0*d*c) * scale;			// Solve for X0 and Y0
				double y0 = ( -2.0*a*e) * scale;			// Which is the origin of the ellips ein user coordinates
																	// Now use the fact that F(x0,y0) = -1 for an ellipse
				f += a*x0*x0  + c*y0*y0 + d*x0 + e*y0;		// This can be used to get the proper scaling
				scale = -1.0/f;										// get the scaling
				U  *= scale;
				m_A = U.At(0);
				m_B = 0;
				m_C = U.At(1);
				m_D = U.At(2);
				m_E = U.At(3);
				m_F = U.At(4);
				UpdateMmatrix();
			}
		}
	}
	if (!ok) Clear();
	return ok;
}


/*---------------------------------------------------------------------------
 *						nxConicEllipse::FromFittedData
 *-------------------------------------------------------------------------*/

nxBOOL nxConicEllipse::GetNearestPoint( double xactual, double yactual, double* xnearest, double* ynearest)
{
	nxBOOL	ok;
	double	xp,yp;
	double	longitude;

	ok    = MapSpaceToParam( xactual, yactual, &xp, &yp);		// Map user coords to translated, and rotated coordinate system

	if (m_a < m_b )
	{
		double temp = yp;
		yp = xp;
		xp = -temp;;
	}

	FromGeocentric( nxVector(xp,0,yp));							// Work out the geodetic latitude, longitude and height
	longitude = (xp  >= 0) ? 0: 180.0;							// now if get  proper longitude
	FromGeodetic( GeodeticLatitude(), longitude, 0 );			// and work out position of 0 height (ie on surface of ellipse)
	xp = Location().X();										// get the X location
	yp = Location().Z();										// get the Y location 

	if (m_a < m_b )
	{
		double temp = yp;
		yp = -xp;
		xp = temp;
	}

	ok    = ok && MapParamToSpace( xp, yp, xnearest, ynearest);	// translate and rotate coordinates
	return ok;
}

/*---------------------------------------------------------------------------
 *						nxConicEllipse::IntersectionWithLine
 *	Finds the intersection of the ellipse with the line
 *
 *	ny + mx + k = 0
 *
 *	Vertical line is given by n = 0 , x = -k/m
 *  Horizontal line is given by m = 0, y = -k/n
 *

 *-------------------------------------------------------------------------*/

nxBOOL nxConicEllipse::IntersectionWithLine( double n, double m, double k, double* ux1, double* uy1, double* ux2, double* uy2)
{
	double A,B,C;
	double s;
	double x1,y1,x2,y2;
	nxBOOL ok;
	nxBOOL	isvertical = (n == 0);

	if ( !isvertical)							// if n is non zero
	{											// then convert to
		m = -m/n;								// y = mx+k form
		k = -k/n;
		A = m_A + m*(m_B + m*m_C);
		B = k*m_B + m*(2*k*m_C +m_E) + m_D;
		C = (m_E+m_C*k)*k + m_F;
		s = (B*B - 4*A*C);						// see
		ok = (s >= 0.0);						// if there is a valid solution
		if (ok)									// if there is
		{										// then
			s = sqrt(s);						// solve the quadratic to get
			x1 = (-B + s) / (2.0*A);			// x1			
			x2 = (-B - s) /	(2.0*A);			// and x2
			y1 = m*x1+k;						// now get y1 and y2;
			y2 = m*x2+k;
		}
	}
	else										// otherwise we have vertical line
	{											// so convert to the form
		k = -k/m;								// x = k
		A = m_C;								// get the quadratic 
		B = k*m_B+m_E;							// coefficients
		C = (m_A*k + m_D)*k + m_F;				// and solve
		s = (B*B - 4*A*C);						// see
		ok = (s >= 0.0);						// if there is a valid solution
		if (ok)									// if there is
		{										// then
			s = sqrt(s);						// solve the quadratic to get
			y1 = (-B + s) / (2.0*A);			// x1			
			y2 = (-B - s) /	(2.0*A);			// and x2
			x1 = k;								// now get y1 and y2;
			x2 = k;
		}
	}
	*ux1 = x1;
	*ux2 = x2;
	*uy1 = y1;
	*uy2 = y2;
	return ok;
}

/*---------------------------------------------------------------------------
 *						nxConicEllipse::IntersectionWithEllipse
 *	Calculates the intersection of two ellipses given a point close to the
 *	intersection.
 *
 *	Works by calculating the gradient at the nearby point and finding intersection
 *	of that with other ellipse.  Repeat process until converged.
 *
 *	returns true if it finds an intersection within "accuracy" distance
 *-------------------------------------------------------------------------*/

nxBOOL nxConicEllipse::IntersectionWithEllipse( nxConicEllipse& other, double approxx,double approxy, double accuracy, double* intersectx, double* intersecty)
{
	double nearx, neary;
	double m; 
	double n;
	double ux1, ux2;
	double uy1, uy2;
	double xint, yint;
	double		d1,d2;
	double	k;
	int			i=0;
	nxBOOL		found = nxFALSE;
	nxBOOL		ok;

	do
	{
		GetNearestPoint( approxx, approxy, &nearx, &neary );					// Get the approximate

		m = -( -m_B*neary - 2.0*m_A*nearx - m_D);								// Get the gradient at the closest point
		n =  m_B*nearx + 2.0*m_C*neary + m_E;									// Get the straight line equation
		k = -m*nearx - n*neary;													// of the gradient at that point

		ok = other.IntersectionWithLine( n, m, k, &ux1, &uy1, &ux2, &uy2);		// Get the intersection of g
		if (ok)
		{
			d1 = sqr(nearx - ux1) + sqr(neary-uy1);				// sqr of distance to 1st solution
			d2 = sqr(nearx - ux2) + sqr(neary-uy2);				// sqr of distance to 2nd solution
			if (d1 < d2) { xint = ux1; yint = uy1;}				// Get straight line intersection closest
			else         { xint = ux2; yint = uy2;}				// original point
			GetNearestPoint( xint, yint, &ux1, &uy1);			// Get closest point on other ellipse
			d1 = sqrt( sqr(xint-ux1) + sqr(yint-uy1));			// get distance between points;
			*intersectx = ux1;									// update the closest point
			*intersecty = uy1;									// if
			found = (d1 <= accuracy);							// and see if its within accuracy;
			if (!found)											// if it is then we are done
			{
				approxx = ux1;
				approxy = uy1;
			}
		}
		i++;
	} while ((i < 40) && ok && !found);
	return found;
}
#endif
