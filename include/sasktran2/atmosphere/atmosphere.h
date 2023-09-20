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

        int num_deriv() const { return 1; }
    };


    /** Essentially void base class for the Atmosphere to remove the NSTOKES parameter for SWIG.
     *
     */
    class AtmosphereInterface {
    public:
        virtual ~AtmosphereInterface() {}
    };

    /** Stores all of the atmosphere information for SASKTRAN2.  Essentially this is extinction/single scatter albedo/phase information
     *  on a grid that matches the global geometry object, as well as surface parameters.  Eventually terms needed
     *  for emission sources likely will be added here as well.
     *
     * @tparam NSTOKES
     */
    template<int NSTOKES>
    class Atmosphere : public AtmosphereInterface {
    private:
        AtmosphereGridStorageFull<NSTOKES> m_storage; /** The internal storage object */
        Surface m_surface; /** The surface */
        bool m_calculate_derivatives; /** True if we are going to be calculating derivatives */
    public:
        /** Directly constructs the atmosphere from it's base objects
         *
         * @param storage
         * @param surface
         * @param calculate_derivatives
         */
        Atmosphere(AtmosphereGridStorageFull<NSTOKES>&& storage,
                   Surface&& surface,
                   bool calculate_derivatives=false);

        virtual ~Atmosphere() {}

        /** Applies delta_m scaling of a specific order to the internal storage object, overwriting it.
         *  Note this is a "half" delta-m scaling.  We scale the extinction/ssa by the regular scaling factors,
         *  and then scale the phase function by 1-f.  This completes the TMS single scatter correction.  For multiple
         *  scatter it is still necessary to further scale the legendre coefficients.
         *
         * @param order
         */
        void apply_delta_m_scaling(int order);

        // Construction from Constituent
        // Atmosphere(const std::vector<Constituent>& constituents, Surface&& surface);

        const AtmosphereGridStorageFull<NSTOKES>& storage() const { return m_storage; };
        AtmosphereGridStorageFull<NSTOKES>& storage() { return m_storage; }
        int num_wavel() const { return (int)m_storage.total_extinction.cols();}

        Surface& surface() { return m_surface;  }
        const Surface& surface() const { return m_surface; }

        // TODO: refactor the below functions into a derivative handler class of some kind
        int ssa_deriv_start_index() const { return (int)m_storage.total_extinction.rows();}
        int scat_deriv_start_index() const { return (int)m_storage.total_extinction.rows() * 2;}
        int surface_deriv_start_index() const { return scat_deriv_start_index() + m_storage.phase[0].num_deriv()*m_storage.total_extinction.rows(); }

        int num_deriv() const;
        int num_scattering_deriv_groups() const { return m_storage.phase[0].num_deriv(); }
    };
}