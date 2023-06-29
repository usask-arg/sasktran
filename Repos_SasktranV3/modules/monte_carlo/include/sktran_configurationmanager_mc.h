#pragma once

#include "sktran_montecarlo_internals.h"

class SKTRAN_ConfigurationManager_MC
{
	private:
		SKTRAN_RayTracingRegionManager							m_raymanager;
		std::shared_ptr<const SKTRAN_CoordinateTransform_V2>	m_coordinatesystem;
	
	public:
																SKTRAN_ConfigurationManager_MC	();
															   ~SKTRAN_ConfigurationManager_MC	();
		void													ReleaseResources				();
		bool													ConfigureCoordinateTransform	( const nxVector& sun, const SKTRAN_LineOfSightArray_V21& linesofsight, double surfaceHeight, double toaHeight, const std::vector<double>& reference_point, bool nadirreferencepointonground );
		const SKTRAN_CoordinateTransform_V2*					CoordinateSystemPtr				() const { return m_coordinatesystem.get();}
		std::shared_ptr<const SKTRAN_CoordinateTransform_V2>	CoordinateSystemObject			() const { return m_coordinatesystem;}
		std::shared_ptr<const SKTRAN_CoordinateTransform_V2>&	CoordinateSystemObjectVar		()		 { return m_coordinatesystem;}
        const SKTRAN_RayTracingRegionManager&					RayManager						() const { return m_raymanager;}

};

