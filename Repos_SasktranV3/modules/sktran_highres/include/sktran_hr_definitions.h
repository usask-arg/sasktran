//#include "sktran_hr_internals.h"

class SKTRAN_HR_Diffuse_Table_Base;

enum	SKTRAN_HR_RayTracer_Type { SKTRAN_HR_RayTracer_Shells,
								   SKTRAN_HR_RayTracer_Curved, 
								   SKTRAN_HR_RayTracer_Curved_NoCurve, 
								   SKTRAN_HR_RayTracer_Straight_Generic };

enum    SKTRAN_HR_Integrator_Type	{ SKTRAN_HR_IntegratorType_Straight,
									  SKTRAN_HR_IntegratorType_Adaptive, 
									  SKTRAN_HR_IntegratorType_Constant };

enum	SKTRAN_HR_OpticalPropertiesTableType { SKTRAN_HR_OpticalPropertiesTableType_1d, 
											   SKTRAN_HR_OpticalPropertiesTableType_3D_UnitSphere, 
											   SKTRAN_HR_OpticalPropertiesTableType_LOSPlane, 
											   SKTRAN_HR_OpticalPropertiesTableType_SZA,
											   SKTRAN_HR_OpticalPropertiesTableType_1d_ConstantLayers };

enum    SKTRAN_HR_wf_Mode { SKTRAN_HR_wf_Mode_1d_uniform,
							SKTRAN_HR_wf_Mode_1d_los,
							SKTRAN_HR_wf_Mode_2d,
							SKTRAN_HR_wf_Mode_None };

enum SKTRAN_HR_wf_precision { SKTRAN_HR_wf_precision_onlylos,
							  SKTRAN_HR_wf_precision_all };

enum SKTRAN_HR_wf_aerosol_Mode { SKTRAN_HR_wf_aerosol_Mode_numberdensity,
	SKTRAN_HR_wf_aerosol_Mode_moderadius,
	SKTRAN_HR_wf_aerosol_Mode_modewidth
};

enum SKTRAN_HR_wf_Species_Mode { SKTRAN_HR_wf_Species_numberdensity,
    SKTRAN_HR_wf_Species_LogNormal_ModeRadius,
    SKTRAN_HR_wf_Species_LogNormal_ModeWidth};

enum    SKTRAN_HR_DiffuseProfilePlacementType { SKTRAN_HR_DiffuseProfilePlacement_LinearSZA,
	SKTRAN_HR_DiffuseProfilePlacement_LinearSZAForceTP,
	SKTRAN_HR_DiffuseProfilePlacement_LinearLOS,
	SKTRAN_HR_DiffuseProfilePlacement_SmartSZA,
};


enum class SKTRAN_HR_PolHOType { unpolarized, constOut };
enum class SKTRAN_HR_AtmosphereHasDelta { unset, yes, no };
enum class SKTRAN_HR_DiffuseMatrixStorageMethod { scalar, scatter_cache, scatter_interpolate, phase };
enum    SKTRAN_HR_DiffuseIncomingSphereType { SKTRAN_HR_DiffuseIncomingType_Default, SKTRAN_HR_DiffuseIncomingType_HardCode, SKTRAN_HR_DiffuseIncomingType_v21Shifted };

typedef float SKTRAN_HR_WEIGHT_TYPE;
#define	SKTRAN_HR_DBL_TO_WEIGHT(x)	(SKTRAN_HR_WEIGHT_TYPE)((x)) 

typedef unsigned int SKTRAN_HR_INDEX_TYPE;
#define SKTRAN_HR_SIZET_TO_INDEX(x) (SKTRAN_HR_INDEX_TYPE)((x))

typedef std::unique_ptr<SKTRAN_RayTracer_Base>                   RayTracerPtr;
typedef std::shared_ptr<SKTRAN_RayFactory_Base>                  RayFactoryPtr;
typedef std::unique_ptr<SKTRAN_SolarTransmission_Base>           SolarTablePtr;
typedef std::unique_ptr<SKTRAN_EmissionTable_Base>               EmissionTablePtr;
typedef std::unique_ptr<SKTRAN_HR_Diffuse_Table_Base>            DiffuseTablePtr;
typedef std::unique_ptr<SKTRAN_OpticalPropertiesIntegrator_Base> OptIntegratorPtr;
typedef std::unique_ptr<SKTRAN_SourceTermIntegrator_Base>        SrcIntegratorPtr;
typedef std::unique_ptr<SKTRAN_TableOpticalProperties_Base>      OpticalTablePtr;