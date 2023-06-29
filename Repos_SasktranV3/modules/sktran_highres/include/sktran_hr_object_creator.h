//#include "sktran_hr_internals.h"

// here some basic functions are provided which create derived
// instances of base objects that may be used in multiple parts of
// the engine, e.g. rays.  All functions are prefixed with SKTRAN_HR
// indicating their global scope

// This is done to easily switch between different ray types
// within the engine without having to create complicated objects
// that would be required to exist within many objects
//
//bool SKTRAN_HR_CreateGeoRay( SKTRAN_RayGeometry_Base** raygeo );
//bool SKTRAN_HR_CreateOptRay( SKTRAN_RayOptical_Base**  rayopt );
