#include <skopticalproperties21.h>
#include <sources/sasktranif_opticalimpl/skbrdf_stubs.h>


bool skBRDF::CheckCosines(double &MU_in, double &MU_out, double &COSDPHI, nxString functionname) const
{
	nxString message;

	// This is the cosine of 89 degrees. We use it to avoid an infinity at 90 degrees
	const double minmu = 0.01745240643728351281941897851632;	

	// force incoming zenith angle to be from 0 to 89
	MU_in = std::max(minmu, std::min(1.0, MU_in));

	// force outgoing zenith angle to be from 0 to 89
	MU_out = std::max(minmu, std::min(1.0, MU_out));

	// force COSDPHI to be within -1 and 1
	if (COSDPHI > 1.0) {
		message = nxString(", the given cosine of the relative azimuth angle is greater than 1. It has been set to 1");
		nxLog::Record(NXLOG_WARNING, functionname + message);
	}
	if (COSDPHI < -1.0) {
		message = nxString(", the given cosine of the relative azimuth angle is less than -1. It has been set to 1");
		nxLog::Record(NXLOG_WARNING, functionname + message);
	}
	COSDPHI = std::max(-1.0, std::min(1.0, COSDPHI));

	return true;
}