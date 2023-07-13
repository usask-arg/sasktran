#if !defined(nxcgsPhysiscalConstants_H)
#define nxcgsPhysiscalConstants_H

/** \ingroup Science
 *	\namespace nxcgs	
 *	Namespace for storing physical constants in the CGS system.  A similar class
 *	with similar variables also exists for the SI system. This allows users to use
 *	the same constant symbols in their code but choose an appropriate namespace
 *
 *	\code
 
 using namespace nxcgs;						// Use physical constants from the CGS namespace

 void main()
 {
	double	distance;
	double	lighttime;

	distance  = ASTRONOMICALUNIT;			// Get one astronomical unit in cgs units
	lighttime = distance/CLIGHT;			// Get the time for light to travel one astronomical unit
	printf("1 AU light time = %f seconds\n", (double)lighttime);
}
	\endcode
 **/

namespace nxcgs
{
	extern       double HCOK;				//!< hc/k in cm*K, used in spectroscopy
	extern const double KBOLTZMAN;			//!< Boltzmanns constant in erg/K
	extern const double CLIGHT;				//!< Speed of light in cm/s
	extern const double HPLANCK;			//!< Plancks constant in erg*s
	extern const double AMU;				//!< One Atomic Mass unit in grams
	extern const double ASTRONOMICALUNIT;	//!< One astronomical unit in cm
	extern const double SOL;				//!< Speed of light in cm/s
};
#endif
