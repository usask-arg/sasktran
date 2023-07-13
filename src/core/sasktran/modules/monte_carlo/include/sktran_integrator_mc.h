#include "sktran_montecarlo_internals.h"


class SKTRAN_OpticalPropertiesIntegrator_Straight_MC : public SKTRAN_OpticalPropertiesIntegrator_Straight
{
	protected:
		// wrappers for interpolating in wavelength on the optical property table
		virtual double			TotalExtinctionPerCM(const SKTRAN_RayOptical_Base* ray, const HELIODETIC_POINT& point) const override;
		virtual bool			GetEffectiveExtinctionPerCMWithHeight1(const SKTRAN_RayOptical_Base* ray, const HELIODETIC_POINT& startpoint, HELIODETIC_POINT& endpoint, double* sigma0, double* sigma1) const override;
		virtual bool			GetEffectiveExtinctionPerCMWithHeight1(const SKTRAN_RayOptical_Base* ray, size_t startPtIndex, double* sigma0, double* sigma1) const override;

	public:
		SKTRAN_OpticalPropertiesIntegrator_Straight_MC() {}
		~SKTRAN_OpticalPropertiesIntegrator_Straight_MC() {}
};

class SKTRAN_OpticalPropertiesIntegrator_Adaptive_MC : public SKTRAN_OpticalPropertiesIntegrator_Adaptive
{
	// untested
	protected:
		virtual double			TotalExtinctionPerCM(const SKTRAN_RayOptical_Base* ray, const HELIODETIC_POINT& point) const override;
		virtual bool			GetEffectiveExtinctionPerCMWithHeight1(const SKTRAN_RayOptical_Base* ray, const HELIODETIC_POINT& startpoint, HELIODETIC_POINT& endpoint, double* sigma0, double* sigma1) const override;
		virtual bool			GetEffectiveExtinctionPerCMWithHeight1(const SKTRAN_RayOptical_Base* ray, size_t startPtIndex, double* sigma0, double* sigma1) const override;

	public:
		SKTRAN_OpticalPropertiesIntegrator_Adaptive_MC() {}
		~SKTRAN_OpticalPropertiesIntegrator_Adaptive_MC() {}
};

class SKTRAN_OpticalPropertiesIntegrator_ConstantLayers_MC : public SKTRAN_OpticalPropertiesIntegrator_ConstantLayers
{
	protected:
		// wrappers for interpolating in wavelength on the optical property table
		virtual double			TotalExtinctionPerCM(const SKTRAN_RayOptical_Base* ray, const HELIODETIC_POINT& point) const override;
		virtual bool			GetEffectiveExtinctionPerCMWithHeight1(const SKTRAN_RayOptical_Base* ray, const HELIODETIC_POINT& startpoint, HELIODETIC_POINT& endpoint, double* sigma0, double* sigma1) const override;
		virtual bool			GetEffectiveExtinctionPerCMWithHeight1(const SKTRAN_RayOptical_Base* ray, size_t startPtIndex, double* sigma0, double* sigma1) const override;

	public:
		SKTRAN_OpticalPropertiesIntegrator_ConstantLayers_MC() {}
		~SKTRAN_OpticalPropertiesIntegrator_ConstantLayers_MC() {}
};


