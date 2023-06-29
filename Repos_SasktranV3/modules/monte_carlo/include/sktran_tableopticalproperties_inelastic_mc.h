#pragma once

#include <omp.h>
#include "sktran_montecarlo_internals.h"

class SKTRAN_TableOpticalProperties_3D_UnitSphere_MC;


class SKTRAN_TableOpticalProperties_Inelastic_Base
{
	public:
															SKTRAN_TableOpticalProperties_Inelastic_Base( ) {}
		virtual												~SKTRAN_TableOpticalProperties_Inelastic_Base( ) {}

		virtual bool										MakeThreadsafeFor( size_t numThreads ) = 0;
		virtual bool										SetPolarizationProperties( std::unique_ptr< SKTRAN_PolarizationProperties_Base >& polprops ) = 0;
		virtual bool										ConfigureGeometry( const SKTRAN_TableOpticalProperties_Base* baseTable ) = 0;
		virtual bool										ConfigureOptical( double wavelen, SKTRAN_AtmosphericOpticalState_V21& opticalstate ) = 0;

		virtual double										InelasticExtinctionPerCM( const double& outgoingwavelength, const HELIODETIC_POINT& point ) const = 0;
		virtual bool										GetScatteringMatrixCM2( double wavelength, const HELIODETIC_POINT& point, double cosAngle, SKTRAN_ScatMat_MIMSNC&  pmatrix ) const = 0;
		virtual bool										GetScatteringCoefficientCM2( double wavelength, const HELIODETIC_POINT& point, double cosangle, SKTRAN_PhaseMatrixScalar* scatcoeff ) const = 0;
		virtual bool										GetIncomingWavelength( const double& outgoingwavelength, const HELIODETIC_POINT& point, const double& randNum, double& incomingwavelength) const = 0;
		virtual bool										GetIncomingWavelength( const std::vector<double>& outgoingwavelength, const size_t& primary, const HELIODETIC_POINT& point, const double& randNum, std::vector<double>& incomingwavelength, std::vector<double>& xsratios ) const = 0;
		virtual bool										GetCosScatteringAngle( const double& outgoingwavelength, const HELIODETIC_POINT& point, const double& randNum, double &cosScatAngle, SKTRAN_ScatMat_MIMSNC* pmatrix = nullptr ) const = 0;
		virtual bool										GetWavelengthRange( double& minwavelength, double& maxwavelength ) const = 0;
};

class SKTRAN_TableOpticalProperties_Inelastic_DoNothing : public SKTRAN_TableOpticalProperties_Inelastic_Base
{
	public:
		virtual bool										MakeThreadsafeFor(size_t numThreads) override { return true; }
		virtual bool										SetPolarizationProperties(std::unique_ptr< SKTRAN_PolarizationProperties_Base >& polprops) override { return true; }
		virtual bool										ConfigureGeometry(const SKTRAN_TableOpticalProperties_Base* baseTable) override { return true; }
		virtual bool										ConfigureOptical(double wavelen, SKTRAN_AtmosphericOpticalState_V21& opticalstate) override { return true; }

		virtual double										InelasticExtinctionPerCM(const double& outgoingwavelength, const HELIODETIC_POINT& point) const override { return 0.0; }
		virtual bool										GetScatteringMatrixCM2(double wavelength, const HELIODETIC_POINT& point, double cosAngle, SKTRAN_ScatMat_MIMSNC&  pmatrix) const override { return false; }
		virtual bool										GetScatteringCoefficientCM2(double wavelength, const HELIODETIC_POINT& point, double cosangle, SKTRAN_PhaseMatrixScalar* scatcoeff) const override { return false; }
		virtual bool										GetIncomingWavelength(const double& outgoingwavelength, const HELIODETIC_POINT& point, const double& randNum, double& incomingwavelength) const { return false; }
		virtual bool										GetIncomingWavelength(const std::vector<double>& outgoingwavelength, const size_t& primary, const HELIODETIC_POINT& point, const double& randNum, std::vector<double>& incomingwavelength, std::vector<double>& xsratios ) const { return false; }
		virtual bool										GetCosScatteringAngle(const double& outgoingwavelength, const HELIODETIC_POINT& point, const double& randNum, double &cosScatAngle, SKTRAN_ScatMat_MIMSNC* pmatrix = nullptr) const override { return false; }
		virtual bool										GetWavelengthRange(double& minwavelength, double& maxwavelength) const override;
};

class SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere : public SKTRAN_TableOpticalProperties_Inelastic_Base
{
	protected:
		const SKTRAN_TableOpticalProperties_Base*				m_baseTable;
		const SKTRAN_TableOpticalProperties_3D_UnitSphere_MC*	m_baseTableMC;
		std::vector< std::vector< std::vector<double> > >		m_inelasticextinction;
		std::vector<double>										m_scattercdf;
		skClimatology*		 									m_clim;
		skOpticalProperties_Inelastic*							m_optProp;
		mutable std::vector<double>								m_wnlookup;
		mutable std::vector<double>								m_cdflookup;
		std::unique_ptr< SKTRAN_PolarizationProperties_Base >	m_inelasticscatprops;

		double				m_wavelen;
		double				m_mjd;
		size_t				m_wavelindex;
		bool				m_firsttime;

		mutable omp_lock_t	m_optPropLock;

	private: // near-copies of private functions from SKTRAN_TableOpticalProperties_3D_UnitSphere (if the const functions could be made protected we wouldn't have to copy them here)
		virtual double			InterpTable(const std::vector< std::vector< std::vector<double> > >& table, double wavelength, const HELIODETIC_POINT& loc) const;
		bool					GetUniquePointWeights( double wavelength, const HELIODETIC_POINT& point, double cosAngle, SKTRAN_GridIndex gridindices[24], double gridweights[24], size_t& numNonZero ) const;
		SKTRAN_GridIndex        TableSubToInd(size_t wavidx, SKTRAN_GridIndex locidx, SKTRAN_GridIndex altidx, SKTRAN_GridIndex angidx) const;

		virtual bool			FillTablesAtIndex(size_t altidx, size_t locidx, SKTRAN_AtmosphericOpticalState_V21& opticalstate);
		virtual bool			FillTablesAtIndexMultiWavel(size_t altidx, size_t locidx, SKTRAN_AtmosphericOpticalState_V21& opticalstate);

	private:
		bool					GetInelasticExtinction( SKTRAN_AtmosphericOpticalState_V21& opticalstate, double wavelength, double* inelxs ) const;
		bool					GetInelasticScatteringMatrix(SKTRAN_AtmosphericOpticalState_V21& opticalstate, double wavelength, double cosscatteringangle, skRTPhaseMatrix* kinel) const;
		bool					makeScatterCdf( SKTRAN_AtmosphericOpticalState_V21& opticalstate );
		void					ReleaseResources();

	public:
								SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere();
		virtual					~SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere();
								
		virtual bool			MakeThreadsafeFor( size_t numThreads ) override { return true; }
		virtual bool			SetPolarizationProperties( std::unique_ptr< SKTRAN_PolarizationProperties_Base >& polprops ) override;
		virtual bool			ConfigureGeometry( const SKTRAN_TableOpticalProperties_Base* baseTable ) override;
		virtual bool			ConfigureOptical( double wavelen, SKTRAN_AtmosphericOpticalState_V21& opticalstate ) override;

		virtual double			InelasticExtinctionPerCM( const double& wavelength, const HELIODETIC_POINT& point ) const override;
		virtual bool			GetScatteringMatrixCM2( double wavelength, const HELIODETIC_POINT& point, double cosAngle, SKTRAN_ScatMat_MIMSNC&  pmatrix ) const;
		virtual bool			GetScatteringCoefficientCM2( double wavelength, const HELIODETIC_POINT& point, double cosangle, SKTRAN_PhaseMatrixScalar* scatcoeff ) const override;

		virtual bool			GetIncomingWavelength(const double& outgoingwavelength, const HELIODETIC_POINT& point, const double& randNum, double& incomingwavelength) const;
		virtual bool			GetIncomingWavelength(const std::vector<double>& outgoingwavelength, const size_t& primary, const HELIODETIC_POINT& point, const double& randNum, std::vector<double>& incomingwavelength, std::vector<double>& xsratios ) const;
		virtual bool			GetCosScatteringAngle( const double& outgoingwavelength, const HELIODETIC_POINT& point, const double& randNum, double &cosScatAngle, SKTRAN_ScatMat_MIMSNC* pmatrix = nullptr ) const override;
		virtual bool			GetWavelengthRange(double& minwavelength, double& maxwavelength) const override;

};
