#include <sasktran2/atmosphere/atmosphere.h>

namespace sasktran2::atmosphere {

	template<int NSTOKES>
	Atmosphere<NSTOKES>::Atmosphere(AtmosphereGridStorageFull<NSTOKES>&& storage, Surface&& surface, bool calculate_derivatives):
		m_storage(storage),
		m_surface(surface),
        m_calculate_derivatives(calculate_derivatives),
        m_storage_accessor(&m_storage)
	{

	}

    template<int NSTOKES>
    int Atmosphere<NSTOKES>::num_deriv() const {
        if(m_calculate_derivatives) {
            return (int)(m_storage.total_extinction.rows() * (2 + m_storage.phase[0].num_deriv()) + m_surface.num_deriv());
        } else {
            return 0;
        }
    }

    template<int NSTOKES>
    void Atmosphere<NSTOKES>::apply_delta_m_scaling(int order) {
        if(order >= m_storage.phase[0].max_stored_legendre()) {
            BOOST_LOG_TRIVIAL(warning) << "Trying to delta scale without the correct number of legendre orders";
            return;
        }

        if(m_storage_delta_scaled == nullptr) {
            int nwavel = (int)m_storage.phase.size();
            int nlocation = (int)m_storage.total_extinction.rows();
            int nleg = m_storage.phase[0].max_stored_legendre();

            m_storage_delta_scaled = std::make_unique<AtmosphereGridStorageFull<NSTOKES>>(nwavel, nlocation, nleg);
        }

        m_storage_accessor = m_storage_delta_scaled.get();

        int stokes_stack;
        if(NSTOKES == 1) {
            stokes_stack = 1;
        } else if (NSTOKES == 3) {
            stokes_stack = 4;
        }

        for(int w = 0; w < m_storage.phase.size(); ++w) {
            // Scale k = (1-wf) k
            (*m_storage_delta_scaled).total_extinction(Eigen::all, w).array() *=
                    (1 - m_storage.ssa(Eigen::all, w).array() *
                         m_storage.phase[w].storage()(order * stokes_stack, Eigen::all).array());

            // scale w = (1-f)/(1-wf)
            (*m_storage_delta_scaled).ssa(Eigen::all, w).array() *=
                    (1 - m_storage.phase[w].storage()(order * stokes_stack, Eigen::all).array()) /
                    (1 - m_storage.ssa(Eigen::all, w).array() *
                         m_storage.phase[w].storage()(order * stokes_stack, Eigen::all).array());

            // We want wP to be scaled by 1/(1-wf), but we have already scaled it by (1-f) / (1-wf), which implies
            // we need to scale the phase moments by 1/(1-f)
            (*m_storage_delta_scaled).phase[w].storage().array() /= (1 - m_storage.phase[w].storage()(order * stokes_stack,
                                                                                                   Eigen::all).array());
        }
    }


    template class Atmosphere<1>;
    template class Atmosphere<3>;


}