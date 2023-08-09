#pragma once

#include "../geometry.h"
#include <sasktran2/atmosphere/grid_storage.h>
#include <sasktran2/atmosphere/constituent.h>


namespace sasktran2::atmosphere {
    /** Class to represent the surface within the atmosphere.  Right now we just have a lambertian surface, plan to 
     *  upgrade this later
     *
     */
    class Surface {
    private:
        Eigen::VectorXd m_albedo;
    public:
        Surface() {}

        Eigen::VectorXd& albedo() { return m_albedo;  };
        const Eigen::VectorXd& albedo() const { return m_albedo; }

        int num_deriv() const { return 0; }
    };


    /** Essentially void base class for the Atmosphere to remove the NSTOKES parameter for SWIG.
     *
     */
    class AtmosphereInterface {
    public:
        virtual ~AtmosphereInterface() {}
    };

    /**
     *
     * @tparam NSTOKES
     */

    template<int NSTOKES>
    class Atmosphere : public AtmosphereInterface {
    private:
        AtmosphereGridStorageFull<NSTOKES> m_storage;
        Surface m_surface;
        bool m_calculate_derivatives;

        AtmosphereGridStorageFull<NSTOKES>* m_storage_accessor;
    public:
        // Direct construction
        Atmosphere(AtmosphereGridStorageFull<NSTOKES>&& storage,
                   Surface&& surface,
                   bool calculate_derivatives=false);

        virtual ~Atmosphere() {}

        void apply_delta_m_scaling(int order);

        // Construction from Constituent
        // Atmosphere(const std::vector<Constituent>& constituents, Surface&& surface);

        const AtmosphereGridStorageFull<NSTOKES>& storage() const { return m_storage; };
        AtmosphereGridStorageFull<NSTOKES>& storage() { return *m_storage_accessor; }

        Surface& surface() { return m_surface;  }
        const Surface& surface() const { return m_surface; }

        int ssa_deriv_start_index() const { return (int)m_storage.total_extinction.rows();}
        int scat_deriv_start_index() const { return (int)m_storage.total_extinction.rows() * 2;}
        int surface_deriv_start_index() const { return scat_deriv_start_index() + m_storage.phase[0].num_deriv(); }

        int num_wavel() const { return (int)m_storage.total_extinction.cols();}
        int num_deriv() const;
    };
}