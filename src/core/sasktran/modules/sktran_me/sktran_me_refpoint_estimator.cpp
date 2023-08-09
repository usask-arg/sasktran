#include "include/sktran_me_internals.h"

namespace sktran_me {
    ReferencePointEstimator::ReferencePointEstimator() {
        m_raymanager = std::make_unique<SKTRAN_RayTracingRegionManager>();
    }

    void ReferencePointEstimator::estimate_reference_point(const SKTRAN_LineOfSightArray_V21 &linesofsight,
                                                           GEODETIC_INSTANT &refpt) {
        m_raymanager->UpdateUndefinedParametersFromLinesOfSight(linesofsight);

        // Populate the latitude/longitude
        m_raymanager->GetReferencePoint(&refpt.latitude, &refpt.longitude);

        // Calculate the mean time
        refpt.mjd = linesofsight.MeanMJD();

        // Refpt we always set to the geoid surface
        refpt.heightm = 0.0;
    }
}