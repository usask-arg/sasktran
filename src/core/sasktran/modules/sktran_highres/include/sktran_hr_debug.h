#pragma once

#define SKTRAN_HR_ADD_SOLAR_DISK_TO_LOS 0

// Directory some debug output files are put in
#define SKTRAN_HR_DEBUG_DIR "C:/development-work/hr_debug/"
//#define SKTRAN_HR_DEBUG_DIR "C:/ARGsoftware/Repos_SasktranV21_MonteCarlo/output/2015/diffusePoints/"
//#define SKTRAN_HR_DEBUG_DIR "C:/ARGsoftware/Repos_SasktranV21_MonteCarlo/output/2014/trash/dpScalar/"

// Dumps outgoing radiances for the diffuse points
#define SKTRAN_HR_DUMP_DIFFUSE 0

// Dumps scattering matrix for every diffuse point
#define SKTRAN_HR_DUMP_SCATTER_MATRIX 0

// Dumps path for curved rays
#define SKTRAN_HR_DUMP_CURVED_PATH 0

// Dumps information about line of sight rays, cell lengths, source function, optical depth, etc.
#define SKTRAN_HR_DUMP_LOS_DIFFUSE 0

// Outputs timing information to the console
#define SKTRAN_HR_VERBOSE_TIMING 0

// Uncouples diffuse profiles, i.e., each profile computes its own diffuse field separate from the others
#define SKTRAN_HR_DECOUPLE_DIFFUSE_PROFILES 0

void DumpRay( SKTRAN_RayOptical_Base* ray, const SKTRAN_TableOpticalProperties_Base& opttable, const SKTRAN_Source_Term& source );
