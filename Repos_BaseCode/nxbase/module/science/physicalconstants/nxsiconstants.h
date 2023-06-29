#if !defined(nxSIPhysiscalConstants_H)
#define nxSIPhysiscalConstants_H

/** \ingroup Science
 *	\namespace nxSI	
 *	Namespace for storing physical constants in the SI system.  A similar class
 *	with similar variables also exists for the cgs system. This allows users to use
 *	the same constant symbols in their code but choose an appropriate namespace
 *
 *	\code
 
 using namespace nxSI;						// Use physical constants from the SI namespace

 void main()
 {
	double	distance;
	double	lighttime;

	distance  = ASTRONOMICALUNIT;			// Get one astronomical unit in SI units
	lighttime = distance/CLIGHT;			// Get the time for light to travel one astronomical unit
	printf("1 AU light time = %f seconds\n", (double)lighttime);
}
	\endcode
 **/

namespace nxSI
{
	extern       double HCOK;				//!< hc/k, a ratio frequently used in spectroscopy in m*K
	extern const double KBOLTZMAN;			//!< Boltzmanns constant in J/K
	extern const double CLIGHT;				//!< Speed of light in m/s
	extern const double HPLANCK;			//!< Planck's constant in J*s
	extern const double AMU;				//!< One Atomic Mass unit in kg
	extern const double ASTRONOMICALUNIT;	//!< one astronomical unit in m.
	extern const double SOL;				//!< Speed of light in m/s
};
#endif
