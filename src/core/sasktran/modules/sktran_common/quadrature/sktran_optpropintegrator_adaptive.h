//#include "sktran_common_internals.h"


class SKTRAN_SolarTransmission_Base;





/*-----------------------------------------------------------------------------
 *					SKTRAN_OpticalPropertiesIntegrator_Adaptive		2014-2-7*/
/** @ingroup odintegrate
 *	Adaptive ray integrator.
 *
 *  \par Motivation
 *  The default ray tracing mode for the SASKTRAN engine involves finding the 
 *  intersections with a set of spherical shells, using these intersection points
 *  the ray is divided into a discrete set of cells.  The advantage of this method
 *  is primarily its simplicity, however there it can also be numerically justified
 *  by noting that (usually) the optical properties is a function only of altitude.
 *  Unfortunately, tracing through spherical shells produces cells that are very large
 *  near the tangent point, where the radiance contribution is often the most important.
 *  Here we outline a new adaptive integration method that aims to solve this problem.

 *  \par Radiance Line Integral For a Cell
 *  The radiance along a ray is the attenuated line integral of the source function,
 *  namely, 
 *  \f[
 *     I(s_0,\hat{p}) = \int_{s_0}^{s_1} \exp{(-\tau(s_0,s))} J(s_0,\hat{p}) ds
 *   \f]
 *  Numerically this integral is done by calculating the radiance contribution of
 *  each cell and summing the results.  Similarly the radiance contribution for one
 *  cell can be written,
 *  \f[
 *     I_{cell} = \exp (-\tau(cellstart)) \int_0^{\Delta s} \exp{(-\tau(0,s)} J(s) ds
 *  \f]
 *	Here \f$s\f$ runs from 0 to \f$\Delta s\f$, the length of the cell.  This formula
 *  evaluates the source function along the cell and attenuates it back to the start
 *  of the cell, and then this result is attenuated back to the observer.
 *
 *  \par Old SASKTRANV21 Method
 *  The source term and extinction are sampled at the midpoint of the cell, and are
 *  assumed to be constant within the cell.  The integral can then be done exactly
 *  using \f$\tau(0,s) = k s\f$ and \f$J(s) = J\f$, then we get the expression
 *  \f[
 *     I_{cell} = \exp(-\tau(cellstart)) (1 - \exp(-k\cdot \Delta s)) / k
 *  \f]
 *
 *  \par Problems with the Old Method
 *  - The source term is not constant along each cell.  Near the tangent point, the source
 *    term varies by a very large amount.
 *  - The optical depth is not linear across the cell, the actual optical depth of each cell
 *    can be quite different than \f$k \cdot \Delta s\f$
 *
 *  \par New Method For Calculating the Source Term
 *  The source term is sampled at the beginning, mid point, and end point of the cell.
 *  Then we can write, \f$J(s) = a + b\cdot s + c \cdot s^2\f$. The integral over the
 *  cell can still be performed analytically.  Furthermore, rather than sampling \f$k\f$
 *  for the cell, we use the already known optical depth of the cell, \f$ k = \sigma / \Delta s\f$.
 *  
 *  \par Adaptive Ray-Tracing
 *  If a cell is found to have an optical depth larger than the maximum allowed cell
 *  optical depth, the cell is split in two.  This combats the problem of having cells
 *  where the optical depth is highly non-linear.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_OpticalPropertiesIntegrator_Adaptive : public SKTRAN_OpticalPropertiesIntegrator_Straight
{
	protected:
		double											m_maxopticaldepth;
		double											m_maxextinctiongradient;
		double											m_maxrayopticaldepthtosplit;

	public:
														SKTRAN_OpticalPropertiesIntegrator_Adaptive();
		virtual bool                                    CalculateRayScalarTransmission_withMinContainer( SKTRAN_RayOptical_Base* ray, double* transmission, bool totaltransmissiononly, bool usecachedtransmission  ) const override;
		virtual bool									CalculateRayScalarTransmission	( SKTRAN_RayOptical_Base*			ray,
																					  double*							transmisison,
																					  bool								totaltransmissiononly,
																					  bool								usecachedtransmission
																					) const override;
//		virtual bool									CalculateRayTransmission_old( SKTRAN_RayOptical_Base* ray, double* transmission, bool totaltransmissiononly, bool usecachedtransmission  ) const;
		
		void									        SetMaxOpticalDepthOfCell    ( double opticaldepth ) { m_maxopticaldepth = opticaldepth; }
		void									        SetMaxExtinctionGradientOfCell( double gradient ) { m_maxextinctiongradient = gradient; }
		void											SetMaxRayOpticalDepthToSplit(double opticaldepth) { m_maxrayopticaldepthtosplit = opticaldepth; }
};


