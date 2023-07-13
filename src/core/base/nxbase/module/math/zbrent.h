/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/

/*----------------------------------------------------------------------------
 *				zbrent2														*/
/**	\ingroup MATH_functions
 *	A templated function that finds the zero of a function in the range x1 and x2 using the
 *	Van Wijngaarden-Dekker-Brent method. Adapted from from Numerical recipes. Ch 9.3.  THe user
 *	is responsible for ensuring there is only one zero crossing point in the specified x range.
 *
 *	\par template <class EVALUATE>
 *		The EVALUATE class is a function object that has a member function \c double \c operator \c()(double x).
 *		It can also be a regular function
 *
 *	\param func
 *		The templated class that evaluates the function at location x with
 *		EVALUATE::operator()(double x)
 *
 *	\param x1
 *		The lower bounding point of the zero.
 *
 *	\param x2
 *		The upper bounding point of the zero. The function evaluated at x1 and x2 must
 *		be of diferent sign.
 *
 *	\param tol
 *		The tolerance of the final solution. I.e. how close to zero is the function
 *
 *	\param zerovalue
 *		Returns the value of x at the zero value. If function fails because beginning
 *		and end are on same side of zero then zerovalue has the value of the function at the
 *		startpoint. If status equals 2 then the zerovalue may not have quite converged to specified
 *		convergence criteria.
 *
 *	\param status
 *		returns -0 if fit is good, 1 if the function does not cross zero in the specified range
 *		and 2 if the fit did not converge.
 *
 *	\param maxiteration
 *		maximum number of iterations.
 *
 *	\return 
 *		true if the fit converged. 
 *
 *>----------------------------------------------------------------------------*/

template <class EVALUATE> bool zbrent2( EVALUATE func, double x1, double x2, double tol, double* zerovalue, int* status, int maxiteration)
{
	int    iter;
	int    ITMAX = maxiteration;
	double EPS = 3.0e-12;

	double a=x1;
	double b=x2;
	double c=0.0;
	double d=0.0;
	double e=0.0;
	double min1,min2;

	double fa=func(a);
	double fb=func(b);

	double fc,p,q,r,s,tol1,xm;
	bool	ok;

	*status = 0;
    ok = (fb*fa < 0.0);
	if (!ok)
	{
		ok = (fb*fa == 0.0);
		if (ok)
		{
			*zerovalue = (fa == 0.0) ? a : b;
			status = 0;
			return ok;
		}
		else
		{
			*zerovalue = fa;
			*status    = 1;
			return ok; 
		}
	}
    fc=fb;
    for (iter=1; iter <= ITMAX; iter++)
	{
		if (fb*fc > 0.0)
		{
			c  = a;
			fc = fa;
			e  = d = b-a;
		}
		if (fabs(fc) < fabs(fb))
		{
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0)
		{
			*zerovalue = b;
			return ok;
        }
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb))
		{
			s=fb/fa;
			if (a == c)
			{
				p=2.0*xm*s;
				q=1.0-s;
			}
			else
			{
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0)  q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2))
			{
				e=d;
				d=p/q;
			}
			else
			{
				d=xm;
				e=d;
			}
		}
		else
		{
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1) b += d;
		else                b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
		fb=func(b);
	}
//	nxLog::Record(NXLOG_WARNING,"::zbrent, Maximum number of iterations exceeded");
	*zerovalue = b;
	*status = 2;
	return false;
}

/*----------------------------------------------------------------------------
 *				zbrent														*/
/**	\ingroup MATH_functions
 *	A templated function that finds the zero of a function in the range x1 and x2 using the
 *	Van Wijngaarden-Dekker-Brent method. Adapted from from Numerical recipes. Ch 9.3.  THe user
 *	is responsible for ensuring there is only one zero crossing point in the specified x range.
 *
 *	\par template <class EVALUATE>
 *		The EVALUATE class is a function object that has a member function \c double \c operator \c()(double x).
 *		It can also be a regular function
 *
 *	\param func
 *		The templated class that evaluates the function at location x with
 *		EVALUATE::operator()(double x)
 *
 *	\param x1
 *		The lower bounding point of the zero.
 *
 *	\param x2
 *		The upper bounding point of the zero. The function evaluated at x1 and x2 must
 *		be of diferent sign.
 *
 *	\param tol
 *		The tolerance of the final solution. I.e. how close to zero is the function
 *
 *	\return 
 *		the location of the zero crossing point in x.  
 *
 *>----------------------------------------------------------------------------*/

template <class EVALUATE> double zbrent( EVALUATE func, double x1, double x2, double tol)
{
	int    iter;
	int    ITMAX = 800;
	double EPS = 3.0e-12;

	double a=x1;
	double b=x2;
	double c=0.0;
	double d=0.0;
	double e=0.0;
	double min1,min2;

	double fa=func(a);
	double fb=func(b);

	double fc,p,q,r,s,tol1,xm;


    if (fb*fa > 0.0) return 0;     // Zero must be bracketed by plus and minus.
    fc=fb;
    for (iter=1; iter <= ITMAX; iter++)
	{
		if (fb*fc > 0.0)
		{
			c  = a;
			fc = fa;
			e  = d = b-a;
		}
		if (fabs(fc) < fabs(fb))
		{
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0)
		{
			return b;
        }
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb))
		{
			s=fb/fa;
			if (a == c)
			{
				p=2.0*xm*s;
				q=1.0-s;
			}
			else
			{
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0)  q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2))
			{
				e=d;
				d=p/q;
			}
			else
			{
				d=xm;
				e=d;
			}
		}
		else
		{
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1) b += d;
		else                b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
		fb=func(b);
	}
	nxLog::Record(NXLOG_WARNING,"::zbrent, Maximum number of iterations exceeded");
	return 0.0;
}
