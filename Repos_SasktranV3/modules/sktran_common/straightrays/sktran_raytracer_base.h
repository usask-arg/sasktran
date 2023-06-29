//#pragma once
//
//#include "sktran_common_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracer_Base		2013-05-22*/
/** @ingroup rays
 *	Class which provides the basic interface for ray tracing.  Derived 
 *  classes must be able to trace rays through the atmosphere, and provide
 *  the required information to a given base optical ray.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_RayTracer_Base
{
	private:
		std::shared_ptr<const SKTRAN_CoordinateTransform_V2>		m_coords;

	public:
																	SKTRAN_RayTracer_Base		(std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords) : m_coords(coords) {};
		virtual													   ~SKTRAN_RayTracer_Base		()	{};
		std::shared_ptr<const SKTRAN_CoordinateTransform_V2>		CoordsObject				() const { return m_coords;}
		virtual bool ConfigureOptical(SKTRAN_AtmosphericOpticalState_V21 *opticalstate, double wavelen_nm, GEODETIC_INSTANT referencepoint) { return true; }; // Not necessary to implement, used for ray tracers that need optical information like curved rays
};



