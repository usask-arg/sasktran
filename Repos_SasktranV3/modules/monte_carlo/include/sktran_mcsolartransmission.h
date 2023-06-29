#pragma once
#include "sktran_montecarlo_internals.h"

class SKTRAN_TableOpticalProperties_MCBase;

class SKTRAN_SolarTransmission_NoTable_reuseRays_MC : public SKTRAN_SolarTransmission_NoTable_reuseRays
{
	protected:
		virtual bool				ScatterCoefficient( const SKTRAN_SourceTermQueryObject_Base& qobj, const SKTRAN_MCPhoton_Base* photon, const std::vector<double>& outgoingwavelengths, std::vector<double>& incomingwavelengths, std::vector<double>& scattCoeffs ) const;
		virtual bool				MonteCarlo_SourceAtPoint( const SKTRAN_SourceTermQueryObject_Base & qobj, const SKTRAN_MCPhoton_Base* photon, const std::vector<double> wavelengths, const std::vector<double>& scattCoeffs, std::vector<SKTRAN_MCPhoton_RadInfo>& sources ) const;

	public:
									SKTRAN_SolarTransmission_NoTable_reuseRays_MC( ) {}
									~SKTRAN_SolarTransmission_NoTable_reuseRays_MC( ) {}

		virtual bool				MonteCarlo_SingleScatteredRadianceAtPoint( const SKTRAN_SourceTermQueryObject_Base & qobj, SKTRAN_MCPhoton_Base* photon ) const;
		virtual bool				MonteCarlo_GroundScatteredRadianceAtPoint( const SKTRAN_SourceTermQueryObject_Base & qobj, SKTRAN_MCPhoton_Base* photon ) const;
};

class SKTRAN_SolarTransmission_NoTable_reuseRays_SolarSpectrum_MC : public SKTRAN_SolarTransmission_NoTable_reuseRays_MC
{
protected:
	std::vector<double>							m_solarSpectrum;
	SKTRAN_GridDefWavelength*					m_wavelengthgrid;

protected:

	virtual bool				MonteCarlo_SourceAtPoint(const SKTRAN_SourceTermQueryObject_Base & qobj, const SKTRAN_MCPhoton_Base* photon, const std::vector<double> wavelengths, const std::vector<double>& scattCoeffs, std::vector<SKTRAN_MCPhoton_RadInfo>& sources) const;

public:
	SKTRAN_SolarTransmission_NoTable_reuseRays_SolarSpectrum_MC() {}
	~SKTRAN_SolarTransmission_NoTable_reuseRays_SolarSpectrum_MC() {}

	bool						SetSolarSpectrum(const std::vector<double>& solarSpectrum);
	bool						SetWavelengthGrid(SKTRAN_GridDefWavelength* wavelengthgrid) { m_wavelengthgrid = wavelengthgrid; if (m_wavelengthgrid != nullptr) m_wavelengthgrid->AddRef(); return true; }

};


class SKTRAN_SolarTransmission_Inelastic_MC : public SKTRAN_SolarTransmission_NoTable_reuseRays_SolarSpectrum_MC
{
protected:
	const SKTRAN_TableOpticalProperties_MCBase*	m_mcOptProps;

protected:
	virtual bool	FillTable_ClassSpecific() override;
	virtual bool	InelasticScatterCoefficient( const SKTRAN_SourceTermQueryObject_Base& qobj, const SKTRAN_MCPhoton_Base* photon, const double& randNum, const std::vector<double>& outgoingwavelengths, std::vector<double>& incomingwavelengths, std::vector<double>& scattCoeffs ) const;

public:
	SKTRAN_SolarTransmission_Inelastic_MC();
	~SKTRAN_SolarTransmission_Inelastic_MC();

	virtual bool	MonteCarlo_SingleScatteredRadianceAtPoint( const SKTRAN_SourceTermQueryObject_Base& qobj, double& radiance ) const override { return false; }
	virtual bool	MonteCarlo_GroundScatteredRadianceAtPoint( const SKTRAN_SourceTermQueryObject_Base& qobj, double& radiance ) const override { return false; }
	//virtual bool    MonteCarlo_SingleScatteredRadianceAtPoint( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, double& radiance ) const override;
	//virtual bool    MonteCarlo_GroundScatteredRadianceAtPoint( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, double& radiance ) const override;
	virtual bool			MonteCarlo_SingleScatteredRadianceAtPoint(const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_MCPhoton_Base* photon) const override;
};


class SKTRAN_SolarTransmission_Ring_MC : public SKTRAN_SolarTransmission_Inelastic_MC
{

public:
	SKTRAN_SolarTransmission_Ring_MC() {}
	~SKTRAN_SolarTransmission_Ring_MC() {}

	virtual bool			MonteCarlo_SingleScatteredRadianceAtPoint(const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_MCPhoton_Base* photon) const override;
	virtual bool			MonteCarlo_GroundScatteredRadianceAtPoint(const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_MCPhoton_Base* photon) const override;

};