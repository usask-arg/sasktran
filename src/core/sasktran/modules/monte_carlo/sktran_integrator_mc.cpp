#include "include/sktran_montecarlo_internals.h"

double SKTRAN_OpticalPropertiesIntegrator_Straight_MC::TotalExtinctionPerCM(const SKTRAN_RayOptical_Base * ray, const HELIODETIC_POINT & point) const
{
	return m_opticalprops->TotalExtinctionPerCM(ray->GetWavelength(), point);
}

bool SKTRAN_OpticalPropertiesIntegrator_Straight_MC::GetEffectiveExtinctionPerCMWithHeight1(const SKTRAN_RayOptical_Base * ray, const HELIODETIC_POINT & startpoint, HELIODETIC_POINT & endpoint, double * sigma0, double * sigma1) const
{
	return m_opticalprops->GetEffectiveExtinctionPerCMWithHeight1(ray->GetWavelength(), startpoint, endpoint, sigma0, sigma1);
}

bool SKTRAN_OpticalPropertiesIntegrator_Straight_MC::GetEffectiveExtinctionPerCMWithHeight1(const SKTRAN_RayOptical_Base * ray, size_t startPtIndex, double * sigma0, double * sigma1) const
{
	return m_opticalprops->GetEffectiveExtinctionPerCMWithHeight1(ray->GetWavelength(), ray->Storage(), startPtIndex, sigma0, sigma1);
}

double SKTRAN_OpticalPropertiesIntegrator_Adaptive_MC::TotalExtinctionPerCM(const SKTRAN_RayOptical_Base * ray, const HELIODETIC_POINT & point) const
{
	return m_opticalprops->TotalExtinctionPerCM(ray->GetWavelength(), point);
}

bool SKTRAN_OpticalPropertiesIntegrator_Adaptive_MC::GetEffectiveExtinctionPerCMWithHeight1(const SKTRAN_RayOptical_Base * ray, const HELIODETIC_POINT & startpoint, HELIODETIC_POINT & endpoint, double * sigma0, double * sigma1) const
{
	return m_opticalprops->GetEffectiveExtinctionPerCMWithHeight1(ray->GetWavelength(), startpoint, endpoint, sigma0, sigma1);
}

bool SKTRAN_OpticalPropertiesIntegrator_Adaptive_MC::GetEffectiveExtinctionPerCMWithHeight1(const SKTRAN_RayOptical_Base * ray, size_t startPtIndex, double * sigma0, double * sigma1) const
{
	return m_opticalprops->GetEffectiveExtinctionPerCMWithHeight1(ray->GetWavelength(), ray->Storage(), startPtIndex, sigma0, sigma1);
}

double SKTRAN_OpticalPropertiesIntegrator_ConstantLayers_MC::TotalExtinctionPerCM(const SKTRAN_RayOptical_Base * ray, const HELIODETIC_POINT & point) const
{
	return m_opticalprops->TotalExtinctionPerCM(ray->GetWavelength(), point);
}

bool SKTRAN_OpticalPropertiesIntegrator_ConstantLayers_MC::GetEffectiveExtinctionPerCMWithHeight1(const SKTRAN_RayOptical_Base * ray, const HELIODETIC_POINT & startpoint, HELIODETIC_POINT & endpoint, double * sigma0, double * sigma1) const
{
	return m_opticalprops->GetEffectiveExtinctionPerCMWithHeight1(ray->GetWavelength(), startpoint, endpoint, sigma0, sigma1);
}

bool SKTRAN_OpticalPropertiesIntegrator_ConstantLayers_MC::GetEffectiveExtinctionPerCMWithHeight1(const SKTRAN_RayOptical_Base * ray, size_t startPtIndex, double * sigma0, double * sigma1) const
{
	return m_opticalprops->GetEffectiveExtinctionPerCMWithHeight1(ray->GetWavelength(), ray->Storage(), startPtIndex, sigma0, sigma1);
}