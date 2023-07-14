/**
 * SASKTRAN TIR Type Definitions
 * 2018-09-13
 */

class SKTRAN_TIR_TableOpticalProperties;
class SKTRAN_TIR_Integrator;

static const double SPEED_OF_LIGHT = 299792458.0;		// speed of light       [m/s]
static const double PLANCK = 6.62607015e-34;			// Planck's constant    [(m^2 kg) / (s)]
static const double BOLTZMANN = 1.380649e-23;			// Boltzmann's constant [(m^2 kg) / (s^2 K)]

enum class RayTracerTypeTIR { straight, curved };
enum class SourceTermOrderTIR { order0, order2 };
enum class OpticalPropertiesIntegratorTypeTIR { straight, adaptive };
enum class AtmosphereDimensionTIR { dim1, dim2 };
enum class MultithreadedDimensionTIR { wavelength, lineofsight };
enum class LayerExtinctionTypeTIR { constant, linearwithheight };
enum class SpeciesWFUnitTIR { numberdensity, vmr };

typedef std::unique_ptr<SKTRAN_RayTracer_Base> RayTracerPtrTIR;
typedef std::shared_ptr<SKTRAN_RayFactory_Base> RayFactoryPtrTIR;
typedef std::unique_ptr<SKTRAN_TIR_Integrator> IntegratorPtrTIR;
typedef std::unique_ptr<SKTRAN_TIR_TableOpticalProperties> OpticalTablePtrTIR;
