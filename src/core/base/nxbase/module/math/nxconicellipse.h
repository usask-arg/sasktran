#include "../science/geodesy/nxgeodetic.h"

/*--------------------------------------------------------------------------
 *						class nxConicEllipse							*/
/** \ingroup Math_fitting
 *	Class developed that least squares fits ellipses to data points. The
 *	ellipse fitting is based upon the papers:-
 *	-  Ellipse-Specific Direct Least-Square Fitting:\n
 *	   Maurizio Pilu, Andrew W. Fitzgibbon, Robert B. Fisher.\n
 *	   http://vision.dai.ed.ac.uk/maurizp/ElliFitDemo/demo.html.\n
 *	-  Direct Least squares fitting of ellipses:\n
 *     A.W. Fitzgibbon, M. Pilu and R.B. Fisher.\n
 *     In International Conference on Pattern Recognition. August 1996.
 *
 *	\par Overview
 *	The generalised conic equation is given by \f$ F = Ax^{2} + Bxy +Cy^{2} +Dx + Ey +F = 0 \f$
 *	and the ellipse condition \f$ 4AC-B*B > 0 \f$ or \f$ (\frac{x'}{a})^{2}   + (\frac{y'}{b})^{2}  = 1 \f$
 *	where: \f$x',y'\f$ have been translated to \f$(x_{0},y_{0})\f$ and rotated by angle \f$\theta\f$
 *	then:
 *	\f[ A = (\frac{\cos(\theta)}{a})^{2} + (\frac{\sin(\theta)}{b})^{2} \f]
 *	\f[ B = 2\cos(\theta)\sin(\theta)(\frac{1}{a^{2}}-\frac{1}{b^{2}})  \f]
 *	\f[ C = (\frac{\sin(\theta)}{a})^{2} + (\frac{\cos(\theta)}{b})^{2} \f]
 *	\f[ D = -2Ax_{0} - y_{0}B \f]
 *	\f[ E = -2Cy_{0} - x_{0}B \f]
 *	\f[ F = Ax_{0}^{2} + By_{0}^{2} + Cy_{0}^{2} -1 \f]
 *
 *  \par Nearest Point determination:
 *	I used the nxGeodetic routines to do this at it seemed to be the best
 *	technique.  It finds the geodetic latitude of each point which corresponds
 *	to the perpendicular distance.
 *	Another, quicker, technique would use the parametrized coordinates
 *	to work out the parametric coordinate.  xparam = a.cos(u), yparam = b.sin(u),
 *	u = atan2d( yparam, xparam ).  The parametric coordinate would then be
 *	the geocentric latitude.
**/

/*     COS(@)^2    SIN(@)^2
 * A = --------  + --------
 *        a^2         b^2
 *
 *
 * B = 2.COS(@).SIN(@).( 1/a^2 - 1/b^2)
 *
 *     SIN(@)^2    COS(@)^2
 * C = --------  + --------
 *        a^2         b^2
 *
 * D = -2.A.x0 - y0.B
 *
 * E = -2.C.y0 - x0.B
 *
 * F = A.x0^2 + B.y0^2 + C.y0^2 -1
 *
 * {/code}
 *
 *
 *>------------------------------------------------------------------------*/

class nxConicEllipse : protected nxGeodetic

{
	private:
		nxBOOL			m_fitrotated;
		double			m_A;		// A coeff from Generalised equation F = Ax^2 + Bxy +Cy^2 +Dx + Ey +F = 0	 x = cos(@) y = sin(@)
		double			m_B;		// B coeff from Generalised equation
		double			m_C;		// C coeff from Generalised equation
		double			m_D;		// D coeff from Generalised equation
		double			m_E;		// E coeff from Generalised equation
		double			m_F;		// F coeff from Generalised equation
		double			m_x0;		// X Origin of ellipse
		double			m_y0;		// Y Origin of ellipse
		double			m_theta;	// Inclination of ellipse "m_a" axis to coordinate X axis
		double			m_a;		// Length of semi-major/minor axis closest to coordinate X axis.
		double			m_b;		// Length of semi/major/minor axis closest to coordinate Y axis

	private:
		void			UpdateMmatrix();	// Updates Ellipse parameters from Generalised conic equation.
		void			Clear();			// Clears all of the ellipse parameters and conic equation parameters.
		nxBOOL			FitRotatedEllipse   ( double*xp, double* yp, int npts );
		nxBOOL			FitNonRotatedEllipse( double*xp, double* yp, int npts );

	public:
						nxConicEllipse();
		void			SetNonRotated() { m_fitrotated = nxFALSE;}
		void			FromParameters   ( double a, double b, double x0, double y0,  double rotationangle );
		nxBOOL			FromDataPoints   ( double*xp, double* yp, int npts );
		nxBOOL			MapParamToSpace  ( double xparam, double  yparam, double* xspace, double* yspace );
		nxBOOL			MapSpaceToParamAngle(double xspace, double  yspace, double* angle);
		nxBOOL			MapSpaceToParam  ( double xspace, double  yspace, double* xparam, double* yparam );
		nxBOOL			GetNearestPoint  ( double xactual, double yactual, double* xnearest, double* ynearest);
		double			Theta  () const	 { return m_theta*nxmath::ONE_RADIAN;}
		double			ALength() const	 { return m_a;}
		double			BLength() const	 { return m_b;}
		double			XOrig  () const  { return m_x0;}
		double			YOrig  () const	 { return m_y0;}
		double			GradientAsAngle(double x, double y);
		nxBOOL			IntersectionWithLine( double n, double m, double k, double* ux1, double* uy1, double* ux2, double *uy2);
		nxBOOL			IntersectionWithEllipse( nxConicEllipse& other, double approxx,double approxy, double accuracy, double* intersectx, double* intersecty);
};

