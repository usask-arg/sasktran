#include "sasktranif_internals.h"


ISKMie::ISKMie(const char* name)
{
	SasktranIF_ClassFactoryLocator	classfactory;

	classfactory.CreateISKMie(name, &m_mie, DllNamePtr());
}


ISKMie::~ISKMie()
{
	if (m_mie != NULL) m_mie->Release();
}

bool ISKMie::SetPropertyScalar(const char* propertyname, double value) {
	if (m_mie != NULL) return m_mie->SetPropertyScalar(propertyname, value);
}

bool ISKMie::SetPropertyArray(const char* propertyname, const double* value, int numpoints) {
	if (m_mie != NULL) return m_mie->SetPropertyArray(propertyname, value, numpoints);
}

bool ISKMie::SetPropertyObject(const char* propertyname, ISKModuleBase* object) {
	if (m_mie != NULL) return m_mie->SetPropertyObject(propertyname, (nxUnknown*) object);
}

bool ISKMie::SetPropertyString(const char* propertyname, const char* str) {
	if (m_mie != NULL) return m_mie->SetPropertyString(propertyname, str);
}


bool ISKMie::Calculate(double lambda, double radius, double refrac_real, double refrac_imag) {
	return m_mie->Calculate(lambda, radius, refrac_real, refrac_imag);
}

double ISKMie::Qext() {
	return m_mie->Qext();
}

double ISKMie::Qsca() {
	return m_mie->Qsca();
}

double ISKMie::Qabs() {
	return m_mie->Qabs();
}

double ISKMie::Cext() {
	return m_mie->Cext();
}

double ISKMie::Csca() {
	return m_mie->Csca();
}

double ISKMie::Cabs() {
	return m_mie->Cabs();
}

void ISKMie::S1(std::complex<double>** s, int* numpoints) {
    *s = m_mie->S1()->UnsafeArrayBasePtr();

    *numpoints = (int)m_mie->S1()->size();
}

void ISKMie::S2(std::complex<double>** s, int* numpoints) {
    *s = m_mie->S2()->UnsafeArrayBasePtr();

	*numpoints = (int)m_mie->S2()->size();
}

void ISKMie::PMom(double** pmom, int* numlegendre) {
	const nx2dArray<double>* p = m_mie->PMom();

    *numlegendre = p->XSize();

	*pmom = m_mie->PMom()->UnsafeArrayBasePtr();
}

void ISKMie::Sforward(double* real, double* imag) {
	std::complex<double> s = m_mie->SForward();
	*real = s.real();
	*imag = s.imag();
}
void ISKMie::SBackward(double* real, double* imag) {
	std::complex<double> s = m_mie->SBackward();
	*real = s.real();
	*imag = s.imag();
}
void ISKMie::TForward(int i, double* real, double* imag) {
	std::complex<double> T = m_mie->TForward(i);
	*real = T.real();
	*imag = T.imag();
}
void ISKMie::TBackward(int i, double* real, double* imag) {
	std::complex<double> T = m_mie->TBackward(i);
	*real = T.real();
	*imag = T.imag();
}

double ISKMie::Spike() {
	return m_mie->Spike();
}