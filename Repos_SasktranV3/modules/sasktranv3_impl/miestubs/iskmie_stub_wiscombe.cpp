#include "../dllimplementation/stdafx.h"
#include "modules/sasktranv3_impl/sktranif_impl_helperclasses.h"

ISKMie_Stub_Wiscombe::ISKMie_Stub_Wiscombe() {
	// Set default scatter angles

	MakeVectorSetFunctions();
	MakeScalarSetFunctions();
}

void ISKMie_Stub_Wiscombe::MakeScalarSetFunctions() {
	AddSetScalarFunction("ipolzn",
		[&, this](double value)
		{
			return m_mie.Set_IPolzn(round(value));
		}
	);
	AddSetScalarFunction("maxlegendre",
		[&, this](double value)
		{
			return m_mie.Set_MaxLegendreMoment(round(value));
		}
	);
}

void ISKMie_Stub_Wiscombe::MakeVectorSetFunctions() {
	AddSetVectorFunction("scatteringangles",
		[&, this](const double* value, int n)
		{
			nx1dArray<double> angles(n, const_cast<double*>(value));

			return m_mie.Set_AnyScatterAngles(angles);
		}
	);
}

bool ISKMie_Stub_Wiscombe::Calculate(double lambda, double radius, double refrac_real, double refrac_imag) {
	bool ok = true;

	ok = ok && m_mie.Set_Wavelength(lambda);
	ok = ok && m_mie.Set_Radius(radius);
	ok = ok && m_mie.Set_RefractiveIndex(refrac_real, refrac_imag);

	return ok;
}

double ISKMie_Stub_Wiscombe::Qext() {
	return m_mie.Qext();
}

double ISKMie_Stub_Wiscombe::Qsca() {
	return m_mie.Qsca();
}

double ISKMie_Stub_Wiscombe::Qabs() {
	return m_mie.Qabs();
}

double ISKMie_Stub_Wiscombe::Cext() {
	return m_mie.Cext();
}

double ISKMie_Stub_Wiscombe::Csca() {
	return m_mie.Csca();
}

double ISKMie_Stub_Wiscombe::Cabs() {
	return m_mie.Cabs();
}

nx1dArray< std::complex<double> >* ISKMie_Stub_Wiscombe::S1() {
	return m_mie.S1();
}

nx1dArray< std::complex<double> >* ISKMie_Stub_Wiscombe::S2() {
	return m_mie.S2();
}

nx2dArray<double>* ISKMie_Stub_Wiscombe::PMom() {
	return m_mie.PMom();
}

std::complex<double> ISKMie_Stub_Wiscombe::SForward() {
	return m_mie.SForward();
}

std::complex<double> ISKMie_Stub_Wiscombe::SBackward() {
	return m_mie.SBackward();
}

std::complex<double> ISKMie_Stub_Wiscombe::TForward(int i) {
	return m_mie.TForward(i);
}

std::complex<double> ISKMie_Stub_Wiscombe::TBackward(int i) {
	return m_mie.TBackward(i);
}

double ISKMie_Stub_Wiscombe::Spike() {
	return m_mie.Spike();
}