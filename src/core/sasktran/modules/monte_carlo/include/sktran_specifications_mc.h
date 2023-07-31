#pragma once

#include "sktran_montecarlo_internals.h"


class SKTRAN_Specifications_MC : public SKTRAN_Specifications_Base
{
	public:
		enum class SunType           {point, randomDisc};
		enum class SolarTableType    {doNothing, noTable, dim2, dim3};
        enum class EmissionTableType {doNothing, dim1};
		//enum class ThermalType       {none, dim2};
		enum class RayTracerType     {shell, curved, generic};
		enum class OptPropIntType    {straight, adaptive, constant};
		enum class EngineType        {sktran, mc, perturbed, tdLatLon, polarized};
		enum class OptTableType      {dim1, dim2_plane, dim3_delaunay, dim1_constant};
		enum class LogType           {none, obsPlane, stDev, radAlongLOS, photAlongLOS, scatPtOnLOS};
		enum class SymmetryType      {none, horiz};
		enum class PolType           {none, pol, pv1};
		enum class AtmoHasDelta      {unset, yes, no};
		enum class ProfileType		 {rayTracing, airMassFactor, solarTable, opticalProperties};
		enum class ScatterType		 {elastic, inelastic, manualInelastic, both, manualBoth};
		enum class WavelengthType	 {single, simultaneous};
		enum class SecondaryOutput	 {none, lengthAMF, opticalDepthAMF, elasticRaman, ringSpectrum, fillingInParameter};

	private: // Types
		EngineType        m_mcEngineType;
		SunType           m_sunType;
        EmissionTableType m_emissionTableType;
		SolarTableType    m_solarTableType;
		//ThermalType       m_thermalTableType;
		RayTracerType     m_solarRayTracerType;
		RayTracerType     m_LOSRayTracerType;
		RayTracerType     m_MSRayTracerType;
		OptPropIntType    m_optPropsIntType;
		OptTableType      m_optTableType;
		LogType           m_photonLogType;
		SymmetryType      m_symmetryType;
		PolType           m_polType;
		AtmoHasDelta      m_atmosphereHasDelta;
		ScatterType		  m_scatterType;
		WavelengthType	  m_wavelengthType;
		SecondaryOutput m_secondary;

	private: // Grid resolutions, etc.
		std::vector<size_t>                         m_numPhotonsPerLOS;         // Maximum number of ray paths to trace per LOS 
		size_t										m_numOptPropAlts;			// ignored if m_manualopticalheights.size() > 0
		std::vector<double>							m_manualopticalheights;
		size_t										m_numRayTracingAlts;		// ignored if m_manualRayTracingHeights.size() > 0
		std::vector<double>							m_manualRayTracingHeights;
		std::vector<double>							m_manualAMFHeights;
        double                                      m_solarTableAltDelta;		// ignored if m_manualSolarTableHeights.size() > 0
		std::vector<double>							m_manualSolarTableHeights;
		size_t										m_numSolTableCosSza;
		size_t                                      m_thread_chunkSize;
		std::uint32_t                               m_rngSeed;                  // Set to zero for runtime-generated "random" seed
		std::vector<double>                         m_precisionMC;              // Desired stdev/measurement for each LOS
		double										m_scattAngleResolution;
		double										m_TOAHeight;
		double										m_curvedRayStepSize;
		double										m_adaptivemaxopticaldepth;
		double                                      m_adaptiveminratio;
		double                                      m_minRelScatterWeight;      // Ray paths that contribute less than this, as a fraction of the radiance already measured along the path, are terminated
		double                                      m_scatterPosRes;            // Scatter positions are accurate to within (0.5,1)*m_scatterPosRes
		double                                      m_sineSunApexAngle;         // The half-apperture of the (unobstructed) sun 
		std::vector<double>                         m_minFractionHigherOrder;   // The probability of scattering to higher order, ignoring VDOS
        double                                      m_groundShiftAltitude;
		bool                                        m_threads_allowDynamic;
        bool                                        m_hasBeenFinalized;
		CLIMATOLOGY_HANDLE							m_amfSpeciesHandle;
		double										m_primaryWavelength;
		size_t										m_primaryWavelengthIndex;
		std::vector<double>							m_radiancewavelengths;
		std::vector<double>							m_opticalpropertieswavelengths;
		std::vector<double>							m_solarSpectrum;
		std::vector<size_t>							m_maxRamanOrders;
		std::string									m_exportFileName;
		size_t										m_chunkSize;
		bool										m_nadirreferencepointonground;
		std::vector<double>							m_referencepoint;							

		
	private: // Set before/when specs are finalized
		nxVector												m_sunposition;
		std::shared_ptr<SKTRAN_GridDefRayTracingShells_V21>		m_raytracingshells;
		std::shared_ptr<SKTRAN_GridDefAirMassFactorShells>		m_amfshells;
		SKTRAN_GridDefScatterAngle_V21*							m_scatteranglegrid;
		SKTRAN_GridDefWavelength*								m_wavelengthgrid;
		SKTRAN_GridDefOpticalPropertiesRadii_V21*				m_opticalpropradii;
		SKTRAN_UnitSphere_V2*									m_opticalunitsphere;

	private:
		void										ReleaseGrids();
		
        bool                                        ConfigureSLonGrid                     ( SKTRAN_GridDefSLON_V21& slongrid, std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords ) const;
		bool										CreateEmissionTable_DoNothing         ( SKTRAN_EmissionTable_Base** target ) const;
		bool										CreateEmissionTable_1DTable           ( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords, SKTRAN_EmissionTable_Base** target ) const;
        bool                                        CreateSolarTransmissionTable_DoNothing(SKTRAN_SolarTransmission_Base** target) const;
		bool										CreateSolarTransmissionTable_NoTable  ( SKTRAN_SolarTransmission_Base**             target ) const;
		bool										CreateSolarTransmissionTable_2DTable  ( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords, SKTRAN_SolarTransmission_Base** target ) const;
        bool										CreateSolarTransmissionTable_3DTable  ( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords, SKTRAN_SolarTransmission_Base** target ) const;
		bool										CreateSolarTransmissionTable_Inelastic_NoTable( SKTRAN_SolarTransmission_Base** target ) const;
		bool										CreateSolarTransmissionTable_Ring_NoTable( SKTRAN_SolarTransmission_Base** target ) const;
		//bool                                        CreateThermalEmissionTable_none       ( SKTRAN_ThermalEmission_2D**                 target ) const;
		//bool                                        CreateThermalEmissionTable_2DTable    ( SKTRAN_ThermalEmission_2D**                 target ) const;
		bool                                        CreateOpticalPropsIntegrator_Straight ( SKTRAN_OpticalPropertiesIntegrator_Base**   target ) const;
		bool                                        CreateOpticalPropsIntegrator_Adaptive ( SKTRAN_OpticalPropertiesIntegrator_Base**   target ) const; 
		bool                                        CreateOpticalPropsIntegrator_Constant ( SKTRAN_OpticalPropertiesIntegrator_Base**   target ) const; 
		bool                                        CreateSun_Point                       ( std::unique_ptr<SKTRAN_Sun_Base>& target ) const;
		bool                                        CreateSun_RandomDisc                  ( std::unique_ptr<SKTRAN_Sun_Base>& target, const std::vector<SKTRAN_RNG>& rngs, size_t numthreads  ) const;
//		bool                                        CreateRayTracer_Shells                ( SKTRAN_RayTracer_Shells**                   target ) const;
//		bool                                        CreateRayTracer_Shells_Curved         ( SKTRAN_RayTracer_Shells_Curved**            target ) const;
		bool                                        CreatePolarizationObject              ( std::unique_ptr< SKTRAN_PolarizationProperties_Base >& target ) const;
		bool										CreateOpticalPropertiesUnitSphere     ( const SKTRAN_CoordinateTransform_V2& coords );
		bool										CreateInelasticTable				  ( std::shared_ptr< SKTRAN_TableOpticalProperties_Inelastic_Base >& target ) const;

		bool										CreateDelaunaySphere( const SKTRAN_CoordinateTransform_V2& coords );
		nxVector									CalcRotatedVector		( const nxVector& tangent, double zenith, double azimuth ) const;
		bool										AddInfoToGenericRayTracer( SKTRAN_RayTracer_Straight_Generic& raytracer, const SKTRAN_CoordinateTransform_V2& coords, const SKTRAN_GridDefRayTracingShells_V21* shells ) const;
        
		bool										GetProfileAlts						  ( SKTRAN_Specifications_MC::ProfileType profileType, std::vector<double>& profile, bool& uniform) const;

	public:
													SKTRAN_Specifications_MC     ( );
												   ~SKTRAN_Specifications_MC     ( );
		void										ReleaseResources             ( );
		bool										Allocate                     ( );
		bool										ConfigureDefaults            ( );
		bool										FinalizeSpecs                ( );
		bool                                        CreateSun                    ( std::unique_ptr<SKTRAN_Sun_Base>& target, const std::vector<SKTRAN_RNG>& rngs, size_t numthreads   ) const;
		bool										CreateConfigurationManager   ( std::unique_ptr<SKTRAN_ConfigurationManager_MC>& target ) const;
		bool										CreateEmissionTable          ( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords, SKTRAN_EmissionTable_Base** target ) const;
		bool										CreateSolarTransmissionTable ( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords, SKTRAN_SolarTransmission_Base** target, size_t numthreads  ) const;
		//bool                                        CreateThermalEmissionTable   ( SKTRAN_ThermalEmission_2D**                target ) const;
		bool										CreateOpticalPropertyTables  ( SKTRAN_TableOpticalProperties_Base** optprop, SKTRAN_TableOpticalProperties_MCBase** mcoptprop, const SKTRAN_CoordinateTransform_V2& coords);
		bool                                        CreateOpticalPropsIntegrator ( SKTRAN_OpticalPropertiesIntegrator_Base**  target ) const;
		bool                                        CreateAveragingKernel        ( SKTRAN_PhotonLog_Base**              target ) const;
		bool                                        CreateScatterOperator        ( std::shared_ptr<SKTRAN_MCScatterOperator_Base>& target ) const;
		bool                                        CreateRayTracers             ( SKTRAN_Engine_MC_V21*                      engine ) const;
		bool										SetRayTracers                ( SKTRAN_Engine_MC_V21*                      engine,  	std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords) const;
		bool										ConfigureSimultaneousWavelengths( SKTRAN_SimultaneousWavelengthManager& simwl) const;
		bool										CreateOptimalScatterSequenceManager( std::unique_ptr<SKTRAN_OptimalScatterSequenceManager_Base>& target );
		bool										CreatePhotonTemplate		 ( std::unique_ptr<const SKTRAN_MCPhoton_Base>& target ) const;

		bool										CreateAirMassFactorCalculator				( std::unique_ptr<SKTRAN_MCAirMassFactorCalculator_Base>& amfcalc ) const;
		bool										CreateAirMassFactorOpticalPropertiesTable	( SKTRAN_TableOpticalProperties_Base* optprop, SKTRAN_TableOpticalProperties_Base** amfoptprop ) const;
		bool										CreateAirMassFactorOpticalPropsIntegrator	( SKTRAN_OpticalPropertiesIntegrator_Base* integrator, SKTRAN_TableOpticalProperties_Base* amfoptprop, SKTRAN_OpticalPropertiesIntegrator_Base** amfintegrator ) const;
		bool										SetAirMassFactorRayTracers					( std::unique_ptr<SKTRAN_MCAirMassFactorCalculator_Base>&, std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords ) const;

		bool										SetSunGeographicPosition (const nxVector& sun);
		bool                                        SetPrecisionMC           (double d);
		bool                                        SetPrecisionMC           ( const std::vector<double>& dvec);
		bool										SetScattAngleResolution  (double d);
		bool										SetCurvedRayStepSize     (double d);
		bool										SetAdaptOptDepthMax      (double d);
		bool										SetAdaptOptDepthMinRatio (double d);
		bool										SetTOAHeight             (double t);
		bool                                        SetMinimumRelPathWeight  (double d); // Allows value zero
		bool                                        SetScatterPositionRes    (double d);
		bool                                        SetSineSunApexAngle      (double d);
		bool                                        SetMinFractionHigherOrder(double d);
		bool                                        SetMinFractionHigherOrder(const std::vector<double>& d);
        bool                                        SetSolarTableAltDelta    (double d);
		bool										SetManualSolTableHeights (const std::vector<double>& heights) { m_manualSolarTableHeights = heights; return true; }
        bool                                        SetGroundShiftAlt        (double d);
		bool                                        SetAllowDynamicThreads   ( bool b) { m_threads_allowDynamic = b; return true;}
		bool                                        SetNumPhotonsPerLOS      (size_t n);
		bool                                        SetNumPhotonsPerLOS      (const std::vector<size_t>& n);
		bool										SetNumRayTracingAlts	 (size_t n) { m_numRayTracingAlts    = n; return true;}
		bool										SetManualRayTracingHeights(const std::vector<double>& heights) { m_manualRayTracingHeights = heights; return true; }
		bool										SetNumOptPropAlts		 (size_t n) { m_numOptPropAlts       = n; return true;}
		bool										SetManualOpticalHeights  (const std::vector<double>& heights) { m_manualopticalheights = heights; return true; }
		bool										SetManualAMFHeights      (const std::vector<double>& heights) { m_manualAMFHeights = heights; return true; }
		bool										SetAMFSpeciesHandle      (CLIMATOLOGY_HANDLE& handle) { m_amfSpeciesHandle = handle; return true; }
		bool										SetNumSolTableCosSza     (size_t n) { m_numSolTableCosSza    = n; return true;}
		bool                                        SetRngSeed               (std::uint32_t n=0) { m_rngSeed   = n; return true;}
		bool										SetOpticalPropertiesWavelengths(const std::vector<double>& wavelengths) { m_opticalpropertieswavelengths = wavelengths; return true;}
		bool										SetRadianceWavelengths	 (const std::vector<double>& wavelengths) { m_radiancewavelengths = wavelengths; return true;}
		bool										SetPrimaryWavelength	 (double d);
		bool										SetMaxRamanOrders		 (const std::vector<size_t>& maxOrders) { m_maxRamanOrders = maxOrders; return true; }
		bool										SetStatisticsExport		 (std::string filename)						  { m_exportFileName = filename; return true; }
		bool										SetChunkSize			 (size_t n) { m_chunkSize = n; return true; }
		bool                                        SetSunType               (SKTRAN_Specifications_MC::SunType        t) { m_sunType            = t; return true;}
		//bool                                        SetThermalEmissionType   (SKTRAN_Specifications_MC::ThermalType    t) { m_thermalTableType   = t; return true;}
        bool                                        SetEmissionTableType     (SKTRAN_Specifications_MC::EmissionTableType t) { m_emissionTableType = t; return true;}
		bool										SetSolarTableType        (SKTRAN_Specifications_MC::SolarTableType t) { m_solarTableType     = t; return true;}
		bool										SetSolarRayTracerType    (SKTRAN_Specifications_MC::RayTracerType  t) { m_solarRayTracerType = t; return true;}
		bool										SetLOSRayTracerType      (SKTRAN_Specifications_MC::RayTracerType  t) { m_LOSRayTracerType   = t; return true;}
		bool										SetMSRayTracerType       (SKTRAN_Specifications_MC::RayTracerType  t) { m_MSRayTracerType    = t; return true;}
		bool                                        SetOptPropIntType        (SKTRAN_Specifications_MC::OptPropIntType t) { m_optPropsIntType    = t; return true;}
		bool                                        SetOptTableType          (SKTRAN_Specifications_MC::OptTableType   t) { m_optTableType       = t; return true;}
		bool                                        SetKernelType            (SKTRAN_Specifications_MC::LogType        t) { m_photonLogType      = t; return true;}
		bool                                        SetSymmetryType          (SKTRAN_Specifications_MC::SymmetryType   t) { m_symmetryType       = t; return true;}
		bool                                        SetPolType               (SKTRAN_Specifications_MC::PolType        t) { m_polType            = t; return true;}
		bool                                        SetAtmosphereHasDelta    (SKTRAN_Specifications_MC::AtmoHasDelta   t) { m_atmosphereHasDelta = t; return true;}
		bool										SetScatterType			 (SKTRAN_Specifications_MC::ScatterType    t) { m_scatterType        = t; return true;}
		bool										SetWavelengthType		 (SKTRAN_Specifications_MC::WavelengthType t) { m_wavelengthType     = t; return true;}
		bool										SetSecondaryOutput		 (SKTRAN_Specifications_MC::SecondaryOutput t) { m_secondary     = t; return true;}
		bool										SetNadirReferencePointOnGround(bool t) { m_nadirreferencepointonground = t; return true; }
		bool										SetReferencePoint(double latitude, double longitude, double height_metres, double mjd) { m_referencepoint = { latitude, longitude, height_metres, mjd }; return true; }

		bool										GetSun                   (nxVector* target) const;
		bool                                        HasBeenFinalized         () const { return m_hasBeenFinalized;        }
		bool                                        GetAllowDynamicThreads   () const { return m_threads_allowDynamic;    }
		double                                      GetPrecisionMC           () const { return m_precisionMC[0];          }
		double										GetScatterAngleResolution() const { return m_scattAngleResolution;    }
		double										GetCurvedRayStepSize     () const { return m_curvedRayStepSize;       }
		double										GetAdaptOptDepthMax      () const { return m_adaptivemaxopticaldepth; }
		double                                      GetAdaptOptDepthMinRatio () const { return m_adaptiveminratio;        }
		double										GetTOAHeight             () const { return m_TOAHeight;           }
		double                                      GetMinimumRelPathWeight  () const { return m_minRelScatterWeight;     }
		double                                      GetScatterPositionRes    () const { return m_scatterPosRes;           }
		double                                      GetSineSunApexAngle      () const { return m_sineSunApexAngle;        }
		const std::vector<double>&                  GetMinFractionHigherOrder() const { return m_minFractionHigherOrder;  }
        double                                      GetSolarTableAltDelta    () const { return m_solarTableAltDelta;      }
		bool                                        GetNumPhotonsPerLOS      ( std::vector<size_t>& vec ) const { vec.resize(m_numPhotonsPerLOS.size()); std::copy(m_numPhotonsPerLOS.begin(), m_numPhotonsPerLOS.end(), vec.begin()); return true;}
        size_t										GetNumRayTracingAlts     () const { return m_numRayTracingAlts;       }
		size_t										GetNumAMFCells           () const { if (m_amfshells == nullptr) return 0; else return m_amfshells->NumCells() - m_amfshells->ExtendedToGround() - m_amfshells->ExtendedToTOA();   }
		size_t										GetNumOptPropAlts        () const { return m_numOptPropAlts;          }
		size_t										GetNumSolTableCosSza     () const { return m_numSolTableCosSza;       }
		size_t                                      GetChunkSize             () const { return m_thread_chunkSize;        }
		std::uint32_t								GetRngSeed               () const { return m_rngSeed;                 }
		const std::vector<size_t>&					GetMaxRamanOrders		 () const { return m_maxRamanOrders;		  }	
		std::string									GetStatisticsExport		 () { return m_exportFileName; }
		size_t 										GetChunkSize			 () { return m_chunkSize; }
		SKTRAN_Specifications_MC::SunType           GetSunType               () const { return m_sunType;                 } 
		//SKTRAN_Specifications_MC::ThermalType       GetThermalEmissionType   () const { return m_thermalTableType;        }
        SKTRAN_Specifications_MC::EmissionTableType GetEmissionTableType     () const { return m_emissionTableType;       }
		SKTRAN_Specifications_MC::SolarTableType    GetSolarTableType		 () const { return m_solarTableType;          }
		SKTRAN_Specifications_MC::RayTracerType     GetSolarRayTracerType    () const { return m_solarRayTracerType;      }
		SKTRAN_Specifications_MC::RayTracerType     GetLOSRayTracerType      () const { return m_LOSRayTracerType;        }
		SKTRAN_Specifications_MC::RayTracerType     GetMSRayTracerType       () const { return m_MSRayTracerType;         }
		SKTRAN_Specifications_MC::OptPropIntType    GetOprPropIntType        () const { return m_optPropsIntType;         }
		SKTRAN_Specifications_MC::OptTableType      GetOptTableType          () const { return m_optTableType;            }
		SKTRAN_Specifications_MC::LogType           GetKernelType            () const { return m_photonLogType;           }
		SKTRAN_Specifications_MC::SymmetryType      GetSymmetryType          () const { return m_symmetryType;            }
		SKTRAN_Specifications_MC::PolType           GetPolType               () const { return m_polType;                 }
		SKTRAN_Specifications_MC::AtmoHasDelta      GetAtmosphereHasDelta    () const { return m_atmosphereHasDelta;      }
        double                                      GetGroundShiftAlt        () const { return m_groundShiftAltitude;     }
		CLIMATOLOGY_HANDLE							GetAMFSpeciesHandle      () const { return m_amfSpeciesHandle;        }
		SKTRAN_Specifications_MC::ScatterType		GetScatterType			 () const { return m_scatterType;             }
		SKTRAN_Specifications_MC::WavelengthType	GetWavelengthType		 () const { return m_wavelengthType;		  }	
		SKTRAN_Specifications_MC::SecondaryOutput	GetSecondaryOutput		 () const { return m_secondary;		  }	
		bool										GetNadirReferencePointOnGround() const { return m_nadirreferencepointonground; }
		const std::vector<double>&					GetReferencePoint		 () const { return m_referencepoint; }

		std::shared_ptr<SKTRAN_GridDefRayTracingShells_V21>	  GetRayTracingShells  ( ) const { return m_raytracingshells;  }
		std::shared_ptr<SKTRAN_GridDefAirMassFactorShells>    GetAMFRayTracingShells() const { return m_amfshells;		   }
		const SKTRAN_GridDefScatterAngle_V21*                 GetScatterAngleGrid  ( ) const { return m_scatteranglegrid;  }
		const SKTRAN_GridDefOpticalPropertiesRadii_V21*       GetOpticalPropRadii  ( ) const { return m_opticalpropradii;  }
		const SKTRAN_UnitSphere_V2*                           GetOpticalUnitSphere ( ) const { return m_opticalunitsphere; }
		const SKTRAN_GridDefWavelength*						  GetWavelengthGridOpticalProperties ( ) const { return m_wavelengthgrid;    }


			// this is probably not the best place to pass in the solar spectrum
		bool					SetSolarSpectrum( const std::vector<double> solarSpectrum ) { m_solarSpectrum = solarSpectrum; return true; }
		const std::vector<double>& GetSolarSpectrum() const { return m_solarSpectrum; }


		 
};

