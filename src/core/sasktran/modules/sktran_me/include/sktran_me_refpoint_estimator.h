#pragma once

namespace sktran_me {
    class ReferencePointEstimator {
    private:
        std::unique_ptr<SKTRAN_RayTracingRegionManager> m_raymanager;
    public:
        ReferencePointEstimator();
        ~ReferencePointEstimator() = default;

        void estimate_reference_point(const SKTRAN_LineOfSightArray_V21& linesofsight, GEODETIC_INSTANT& refpt);
    };
}