#pragma once

#include "sktran_montecarlo_internals.h"

#include <omp.h>

#include "boost/filesystem.hpp"	// used in first hackey version of scatter cloud -- get rid of this in final version


class SKTRAN_Engine_DebugLog
{
	private:
		double m_pathDistance;

	public:
		SKTRAN_Engine_DebugLog  ( ) { SetDebugDistance( -1.0 );}
		void   SetDebugDistance ( double d ) { m_pathDistance = d;}
		double GetDebugDistance ( ) const { return m_pathDistance;}		
};

class SKTRAN_SimultaneousWavelengthManager
{
	private:
		bool				m_active; // true if calculating multiple wavelengths simultaneously
		bool				m_complete; // true if the simultaneous calculation is complete
		std::vector<double> m_wavelengths;
		double				m_primaryWavelength;
		size_t				m_primaryWavelengthIndex;
		double				m_currentWavelength;
		size_t				m_currentWavelengthIndex;
	public:
		SKTRAN_SimultaneousWavelengthManager() {}
		~SKTRAN_SimultaneousWavelengthManager() {}

		void Configure(std::vector<double> wavelengths, double primaryWavelength, size_t primaryWavelengthIndex)
		{ 
			m_wavelengths = wavelengths; 
			m_primaryWavelength = primaryWavelength; 
			m_active = wavelengths.size() > 0;
			m_primaryWavelength = primaryWavelength;
			m_primaryWavelengthIndex = primaryWavelengthIndex;
			m_complete = false;
		}
		bool						SetCurrentWavelength(double wavelength) { 
			bool ok = true;
			m_currentWavelength = wavelength;
			if (m_active)
			{
				auto it = std::find(m_wavelengths.begin(), m_wavelengths.end(), wavelength);
				ok = ok && it != m_wavelengths.end();
				m_currentWavelengthIndex = it - m_wavelengths.begin(); // will be an invalid index if wavelength is not found
			}
			else
			{
				m_currentWavelengthIndex = 0;
				m_primaryWavelength = wavelength;
				m_primaryWavelengthIndex = 0;
			}
			return ok;
		}
		void						SetComplete(bool complete) { m_complete = complete; }

		double						GetCurrentWavelength()		const { return m_currentWavelength; }
		size_t						GetCurrentWavelengthIndex() const { return m_currentWavelengthIndex; }
		double						GetPrimaryWavelength()		const { return m_primaryWavelength; }
		size_t						GetPrimaryWavelengthIndex() const { return m_primaryWavelengthIndex; }
		const std::vector<double>&	GetWavelengths()			const { return m_wavelengths; }
		bool						Active()					const { return m_active; }
		bool						Complete()					const { return m_complete; }
		size_t						GetNumWavelengths()			const { if (m_active) return m_wavelengths.size(); else return 1; }
};

/*-----------------------------------------------------------------------------
 *					SKTRAN_Engine_MC_V21			2011-07-08*/
/** This is the class used for calculating LOS spectral radiance using Monte
 *  Carlo calculations (a backward method, as opposed to the SASKTRAN forward
 *  method).  It is initialized the same way as #SKTRANSO_Engine, and calculates
 *  the radiance for an array of lines of sight for a specified wavelength.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_Engine_MC_V21 : public SKTRAN_Engine_Base
{		
	private:
		std::unique_ptr<const SKTRAN_MCPhoton_Base> m_photonTemplate;
		std::vector<std::unique_ptr<SKTRAN_MCPhoton_Base>> m_higherOrderPhotons;         //!< Higher order photons' pdfs are not cached -- they change after each scatter
		std::vector<std::unique_ptr<SKTRAN_MCPhoton_Base>> m_firstOrderPhotons;      //!< First order scattered photons' pdfs are cached -- same ray coming out of observer every time
		SKTRAN_TableOpticalProperties_Base*         m_opticalpropertiestable; //!< A table which contains all optical properties for the atmosphere
		
		SKTRAN_LineOfSightArray_V21					m_linesofsight;
		double										m_prevWavelen;
		nx2dArray<double>							m_radiances;				// the photon weights for all los's
		nx2dArray<double>							m_variances;				// Total signal variance (in I) for each LOS
		nx2dArray<double>							m_secondaryMeasurements;
		nx2dArray<double>							m_secondaryMeasurementVariances;

		nx2dArray<skRTStokesVector>                 m_svecs;
		vector<size_t>								m_numPhotons;
		double                                      m_targetStd;
		double										m_minimumRelativeWeight;
		size_t										m_maxNumScatters;
		double										m_scatterPosResolution;		// Scatter positions are guaranteed to be more accurate than this [meters]
		size_t										m_chunkSize;

		std::vector< SKTRAN_RNG >                   m_randGens;
		bool										m_optStateAltered;
		//std::unique_ptr< SKTRAN_Sun_Base >          m_sun;
		std::shared_ptr< SKTRAN_MCScatterOperator_Base > m_scatterop;
		//SKTRAN_SourceTermIntegrator_Order2          m_integratorOrder2;
		
		std::unique_ptr< SKTRAN_ConfigurationManager_MC > m_mcconfig;
		int											m_externalNumThreads;
		int											m_internalNumThreads;
		
		bool curveAllRays, curveLOSRays;


		SKTRAN_OpticalPropertiesIntegrator_Base*    m_opticalpropsintegrator;

//		SKTRAN_RayTracer_Shells*					m_raytracer_shells;
//		SKTRAN_RayTracer_Shells_Curved*				m_raytracer_shells_curved;
		
		std::shared_ptr<SKTRAN_RayFactory_Base>		m_rayFactory_los;
		std::shared_ptr<SKTRAN_RayFactory_Base>		m_rayFactory_secondary;
		std::shared_ptr<SKTRAN_RayFactory_Base>		m_rayFactory_solar;

		std::string									m_saveDir;
        SKTRAN_EmissionTable_Base*                  m_emissionTable;
		SKTRAN_SolarTransmission_Base*				m_solarTransmissionTable;
		//SKTRAN_ThermalEmission_2D*					m_ThermalEmission;

		SKTRAN_TableOpticalProperties_MCBase*		m_mcOptTable;

		SKTRAN_PhotonLog_Base*                m_aveKernel;

		std::vector<std::vector<double> >			m_sigmaks;
		std::vector<std::vector<double> >			m_sigmafs;
		std::vector<HELIODETIC_UNITVECTOR>			m_dpIncs;
		double m_saveAlb;

		std::vector<double>                                     m_minFractionHigherOrder;
		
		SKTRAN_Engine_DebugLog srdLog;

        nx2dArray<double> m_lastRunTiming;

		nx2dArray<double>										m_airMassFactors;			// x = los idx | y = layer idx
		nx2dArray<double>										m_airMassFactorVariances;

		std::unique_ptr<SKTRAN_MCAirMassFactorCalculator_Base>	m_amfcalculator;
		SKTRAN_TableOpticalProperties_Base*						m_amfopticalpropertiestable;
		SKTRAN_OpticalPropertiesIntegrator_Base*				m_amfopticalpropsintegrator;

		SKTRAN_SimultaneousWavelengthManager					m_simwl;
		std::unique_ptr<SKTRAN_OptimalScatterSequenceManager_Base>	m_seqManager;

	private:
		bool                            ConfigureModel_SetThreads( const SKTRAN_Specifications_MC* mcspecs, size_t numthreads );

	protected:
		virtual	void					InitializeEngine				( );
		virtual bool					CreateOpticalPropertyTables		( SKTRAN_Specifications_MC* mcspecs);
		bool							InitializeRandomNumberGenerators( const std::uint32_t seed = 0);
		bool							MonteCarloMultiScatter			( size_t losIdx);	// Integrate by sampling uncorrelated elements of path space
		//void							CreateDiffuseIncomingPoints     ( );
		bool							TraceMotherRays					( );
		void							DeleteMotherRays				( );

		bool                            ChangeScatterDirection      ( const SKTRAN_MCScatterOperator_Base* scatterop, const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* mcphoton, size_t threadid, size_t orderOfScatter ); // Should be made const
		bool                            ChooseScatterPoint          ( SKTRAN_MCPhoton_Base const* incomingPhoton, SKTRAN_MCPhoton_Base* mcphoton, HELIODETIC_POINT& scatterPoint, size_t order ) const;
		bool                            MultipleScatterContribution (	SKTRAN_MCScatterOperator_Base*		scatterop, 
																		SKTRAN_MCPhoton_Base const* const		motherPhoton,
			                                                            SKTRAN_MCPhoton_Base*                    resultDest, 
																		HELIODETIC_VECTOR&                  prevScatterPoint, 
																		size_t scatterSequenceIndex, 
																		int&                                orderOfScatter, 
																		size_t                              threadid
																    );
		
		virtual bool					CalculateOpticalPropertiesTable( double								 wavelen,
																		 SKTRAN_AtmosphericOpticalState_V21* opticalstate,
																		 bool								 userupdateclimatology );
		virtual bool					CalculateAMFOpticalPropertiesTable ( double								 wavelen,
																			 SKTRAN_AtmosphericOpticalState_V21* opticalstate,
																			 bool								 userupdateclimatology );	
		void							ReleaseResources();
		

	public:
										SKTRAN_Engine_MC_V21 	();
		virtual							~SKTRAN_Engine_MC_V21	();
		GEODETIC_INSTANT                ReferencePoint          () const;
		bool							SetRayFactory_LOS		( std::shared_ptr<SKTRAN_RayFactory_Base>	rayFactory_los);
		bool							SetRayFactory_SOLAR		( std::shared_ptr<SKTRAN_RayFactory_Base>	rayFactory_solar);
		bool							SetRayFactory_SECONDARY	( std::shared_ptr<SKTRAN_RayFactory_Base>	rayFactory_secondary);
		void							SetOptStateAltered		(const bool& value)			{ m_optStateAltered = value;}
		bool                            ConfigureMinorResolutionUpdates ( SKTRAN_SpecsUser_Base& userspecs_base );
		virtual bool					ConfigureModel			(      SKTRAN_SpecsUser_Base&			modelspecifications,
																 const SKTRAN_LineOfSightArray_V21&		linesofsight,
																 size_t									numthreads = 0);

		virtual bool					CalculateRadiance		(std::vector<SKTRAN_StokesScalar>*		losradiance,
																 double									wavelen,
																 size_t									numordersofscatter,
																 SKTRAN_AtmosphericOpticalState_V21*	opticalstate,
																 std::vector<skRTStokesVector>*         losvector = nullptr,
																 bool									updateclimatology = false,
																 SKTRAN_DiagnosticInterface*			diag = NULL) override;

		bool							GetStokesVectors	    ( nx1dArray<skRTStokesVector>& dest ) const;
		bool                            GetMeasurementVariance  ( size_t losidx, double& variance ) const;
		bool                            GetAirMassFactors       ( size_t losidx, std::vector<double>& amf ) const;
		bool                            GetAirMassFactorVariance( size_t losidx, std::vector<double>& amfvar) const;
		bool							GetSecondaryMeasurement ( size_t losidx, double& measurement) const;
		bool							GetSecondaryVariance	( size_t losidx, double& variance) const;

		//bool							GenerateRandomNumber	( SKTRAN_RNG rng );

		const SKTRAN_CoordinateTransform_V2* Coordinates		() const {return m_mcconfig->CoordinateSystemPtr();}
		// Debugging functions mostly
		void							SetSaveDir				( std::string dir) {m_saveDir = dir; boost::filesystem::create_directory(m_saveDir); }
		void							SetAlbForSave			(double d) {m_saveAlb = d;}

		const SKTRAN_LineOfSightArray_V21 LinesOfSight() { return m_linesofsight; }

        bool GetTimingData ( nx2dArray<double>& t ) const;

		friend bool SKTRAN_Specifications_MC::CreateRayTracers            ( SKTRAN_Engine_MC_V21* engine ) const;
		friend bool SKTRAN_Specifications_MC::SetRayTracers               ( SKTRAN_Engine_MC_V21* engine, std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords) const;	// This is ugly, should probably be done internally to specs
};		
