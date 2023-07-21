#pragma once
#include "sktran_disco/sktran_do.h"
#include <sasktran2/geometry.h>

namespace sasktran_disco {
    template <int NSTOKES, int CNSTR=-1>
    using GeometryLayerArrayROP = ReadOnlyProperties<
    BasicProperties<NSTOKES>,
    SolarProperties<NSTOKES>,
    UserSpecProperties,
    TestProperties
    >;

    // Class which calculates and stores the wavelength independent aspects of constructing the OpticalLayerArray
    template <int NSTOKES, int CNSTR=-1>
    class GeometryLayerArray : public GeometryLayerArrayROP<NSTOKES>, public AzimuthDependencyCascade {
    protected:
        const PersistentConfiguration<NSTOKES, CNSTR>& m_config;
        Eigen::MatrixXd								   m_chapman_factors;
        Eigen::MatrixXd                                m_optical_interpolator;
        Eigen::VectorXd                                m_floor_h;
        Eigen::VectorXd                                m_ceiling_h;

        GeometryLayerArray(const PersistentConfiguration<NSTOKES, CNSTR>& config) : GeometryLayerArrayROP<NSTOKES>(config), m_config(config) {};


        void calculate_chapman_factors(double earth_rad);
    public:
        GeometryLayerArray(const PersistentConfiguration<NSTOKES, CNSTR>& config,
            const sasktran_disco_lowlevel::Atmosphere& atmosphere
        );

        GeometryLayerArray(const PersistentConfiguration<NSTOKES, CNSTR>& config,
                           const sasktran2::Geometry1D& geometry
                           );

        const Eigen::MatrixXd& chapman_factors() const { return m_chapman_factors; }
        const Eigen::MatrixXd& interpolating_matrix() const { return m_optical_interpolator; }

        const Eigen::VectorXd& layer_floor() const { return m_floor_h; }
        const Eigen::VectorXd& layer_ceiling() const { return m_ceiling_h;}
    };
}