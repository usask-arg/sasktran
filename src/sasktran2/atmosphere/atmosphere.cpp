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

        int stokes_stack;
        if(NSTOKES == 1) {
            stokes_stack = 1;
        } else if (NSTOKES == 3) {
            stokes_stack = 4;
        }

        m_storage.applied_f_location = order * stokes_stack;
        m_storage.applied_f_order = order;

        // First copy the requested order into the f storage
        for(int w = 0; w < m_storage.phase.size(); ++w) {
            m_storage.f(Eigen::all, w) = m_storage.phase[w].storage()(m_storage.applied_f_location, Eigen::all) / (2*order + 1);
        }

        // Now scale k* = (1-wf) k
        m_storage.total_extinction.array() = m_storage.total_extinction.array() * (1 - m_storage.ssa.array() * m_storage.f.array());

        // And then w* = (1 - f) / (1-wf) * w
        m_storage.ssa.array() = (1 - m_storage.f.array()) / (1 - m_storage.ssa.array() * m_storage.f.array()) * m_storage.ssa.array();

        // Lastly we need to scale the legendre coefficients, b = b / (1-f)
        // Note that this is the scaling for single scatter TMS correction, we still have to subtract f/1-f for multiple scatter
        for(int w = 0; w < m_storage.phase.size(); ++w) {
            m_storage.phase[w].storage().array().rowwise() /= (1 - m_storage.f(Eigen::all, w).transpose().array());

            // Also have to scale the derivatives
            for(int d = 0; d < m_storage.phase[w].num_deriv(); ++d) {
                // Set the f derivative
                m_storage.phase[w].f_deriv_storage(d).array() = m_storage.phase[w].deriv_storage(d)(m_storage.applied_f_location, Eigen::all).array() / (2*order + 1);

                // db* = db / 1-f + (b*) * df / (1-f)
                m_storage.phase[w].deriv_storage(d).array() += m_storage.phase[w].storage().array().rowwise() * (m_storage.phase[w].f_deriv_storage(d).array().transpose());
                m_storage.phase[w].deriv_storage(d).array().rowwise() /= (1 - m_storage.f(Eigen::all, w).transpose().array());
            }
        }
    }


    template class Atmosphere<1>;
    template class Atmosphere<3>;


}