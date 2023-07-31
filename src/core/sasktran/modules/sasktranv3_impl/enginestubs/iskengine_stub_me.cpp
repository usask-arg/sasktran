#include "../dllimplementation/stdafx.h"
#include "modules/sasktranv3_impl/sktranif_impl_helperclasses.h"
#include <boost/regex.hpp>



bool ISKEngine_Stub_ME::AddLineOfSight(double mjd, const nxVector& observer, const nxVector& lookvector, int* losindex) {

	return true;
}

bool ISKEngine_Stub_ME::AddSpecies(const CLIMATOLOGY_HANDLE& species, ISKClimatology_Stub* climatology, ISKOpticalProperty_Stub* opticalproperty) {
	return true;
}

bool ISKEngine_Stub_ME::AddEmission(const EMISSION_HANDLE& species, ISKEmission_Stub* emission) {

	return true;

}

bool ISKEngine_Stub_ME::SetAtmosphericState(ISKClimatology_Stub* climatology) {

	return true;


}

bool ISKEngine_Stub_ME::SetAlbedo(double albedo) {

	return true;


}

bool ISKEngine_Stub_ME::SetBRDF(ISKBrdf_Stub* brdf) {

	return true;


}

bool ISKEngine_Stub_ME::SetPolarizationMode(int polarizationmode) {

	return true;


}

bool ISKEngine_Stub_ME::SetWavelengths(const double* wavelen, int numwavelen) {

	return true;


}

bool ISKEngine_Stub_ME::InitializeModel() {
	// 

	return true;

}


bool ISKEngine_Stub_ME::CalculateRadiance(const double** radiance, int* numwavelens, int* numlinesofsight) {
	return true;


}

bool ISKEngine_Stub_ME::CalculateStokesVector(const ISKStokesVector** radiancep, int* numwavelens, int* numlinesofsight) {
	return true;


}

bool ISKEngine_Stub_ME::GetWeightingFunctions(const double** wf, int* numwavel, int* numlinesofsight, int* numwf) {
	return true;

}
