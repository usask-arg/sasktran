#include <sasktran2/config.h>

namespace sasktran2 {
    Config::Config() :
        m_nstokes(1),
        m_nthreads(1),
        m_ndostreams(16),
        m_enable_wfs(true),
        m_apply_delta_scaling(true),
        m_wf_precision(WeightingFunctionPrecision::full),
        m_nsinglescatter_moments(16),
        m_ndosza(1),
        m_ndosphericaliterations(0),
        m_hr_nincoming(110),
        m_hr_noutgoing(110),
        m_hr_nspherical_iterations(50),
        m_hr_num_incoming_points(-1),
        m_initialize_hr_with_do_solution(false)
    {
        set_multiple_scatter_source(MultipleScatterSource::none);
        set_single_scatter_source(SingleScatterSource::exact);
        set_occultation_source(OccultationSource::none);
    }


}