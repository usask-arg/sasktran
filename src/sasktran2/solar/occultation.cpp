#include <sasktran2/solartransmission.h>
#include <sasktran2/dual.h>
#include <sasktran2/config.h>
#include <sasktran2/atmosphere/atmosphere.h>

namespace sasktran2::solartransmission {
    template <int NSTOKES>
    void OccultationSource<NSTOKES>::initialize_config(const sasktran2::Config &config) {

    }

    template <int NSTOKES>
    void OccultationSource<NSTOKES>::initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay> &los_rays) {

    }

    template <int NSTOKES>
    void OccultationSource<NSTOKES>::integrated_source(int wavelidx, int losidx, int layeridx, int threadidx,
                                                       const sasktran2::raytracing::SphericalLayer &layer,
                                                       const sasktran2::SparseODDualView& shell_od,
                                                       sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES> &source) const {

    }

    template <int NSTOKES>
    void OccultationSource<NSTOKES>::calculate(int wavelidx, int threadidx) {

    }

    template <int NSTOKES>
    void OccultationSource<NSTOKES>::initialize_atmosphere(
            const sasktran2::atmosphere::Atmosphere<NSTOKES> &atmosphere) {

    }

    template <int NSTOKES>
    void OccultationSource<NSTOKES>::end_of_ray_source(int wavelidx, int losidx, int threadidx,
                                                       sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES> &source) const {
        if constexpr(NSTOKES == 1) {
            source.value.array() += 1;
        } else {
            source.value(0) += 1;
        }
    }

    template class OccultationSource<1>;
    template class OccultationSource<3>;


}