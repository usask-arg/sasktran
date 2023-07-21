#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_geometrylayerarray.h"

template <int NSTOKES, int CNSTR>
sasktran_disco::GeometryLayerArray<NSTOKES, CNSTR>::GeometryLayerArray(const PersistentConfiguration<NSTOKES, CNSTR>& config,
                                                                       const sasktran_disco_lowlevel::Atmosphere& atmosphere
) : m_config(config),
GeometryLayerArrayROP<NSTOKES>(config)
{
    // Initialize the chapman factor storage
    m_chapman_factors.resize(this->M_NLYR, this->M_NLYR);
    m_chapman_factors.setZero();

    // Calculate the ceiling heights and floor heights of each layer
    m_ceiling_h.resize(this->M_NLYR);
    m_floor_h.resize(this->M_NLYR);

    for (int p = 0; p < (int)this->M_NLYR; p++) {
        m_ceiling_h(p) = atmosphere.layerboundaryaltitude[p];

        if (p == this->M_NLYR - 1) {
            m_floor_h(p) = 0.0;
        }
        else {
            m_floor_h(p) = atmosphere.layerboundaryaltitude[p + 1];
        }
    }

    // Now we have the layer locations, we can calculate the chapman factors
    calculate_chapman_factors(atmosphere.earthradius);
}

template <int NSTOKES, int CNSTR>
sasktran_disco::GeometryLayerArray<NSTOKES, CNSTR>::GeometryLayerArray(const PersistentConfiguration<NSTOKES, CNSTR>& config,
                                                                       const sasktran2::Geometry1D& geometry
) : m_config(config),
    GeometryLayerArrayROP<NSTOKES>(config)
{
    // Initialize the chapman factor storage
    m_chapman_factors.resize(this->M_NLYR, this->M_NLYR);
    m_chapman_factors.setZero();

    // Calculate the ceiling heights and floor heights of each layer
    m_ceiling_h.resize(this->M_NLYR);
    m_floor_h.resize(this->M_NLYR);

    for (int p = 0; p < (int)this->M_NLYR; p++) {
        m_ceiling_h(p) = geometry.altitude_grid().grid().reverse()(p);

        if (p == this->M_NLYR - 1) {
            m_floor_h(p) = 0.0;
        }
        else {
            m_floor_h(p) = geometry.altitude_grid().grid().reverse()(p+1);
        }
    }

    // Now we have the layer locations, we can calculate the chapman factors
    calculate_chapman_factors(geometry.coordinates().earth_radius());
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::GeometryLayerArray<NSTOKES, CNSTR>::calculate_chapman_factors(double earth_rad) {
    m_chapman_factors.setZero();
    if (m_config.use_pseudo_spherical()) {
        // Calculate the chapman factors for the layer
        double sinthetasq = 1 - this->M_CSZ * this->M_CSZ;

        for(sasktran_disco::LayerIndex p = 0; p < this->M_NLYR; ++p) {
            double rp = earth_rad + m_floor_h(p);

            for (sasktran_disco::LayerIndex q = 0; q <= p; ++q) {
                double rfloor = earth_rad + m_floor_h(q);
                double rceil = earth_rad + m_ceiling_h(q);

                m_chapman_factors(p, q) = (sqrt(rceil * rceil - rp * rp*sinthetasq) - sqrt(rfloor * rfloor - rp * rp*sinthetasq)) / (rceil - rfloor);
            }
        }
    }
    else {
        m_chapman_factors.setConstant(1 / this->M_CSZ);
    }
}

SASKTRAN_DISCO_INSTANTIATE_TEMPLATE(sasktran_disco::GeometryLayerArray);