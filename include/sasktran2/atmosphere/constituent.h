#pragma once

#include <sasktran2/internal_common.h>
#include <sasktran2/grids.h>
#include <sasktran2/geometry.h>
#include <sasktran2/atmosphere/grid_storage.h>

namespace sasktran2::atmosphere {
    class Constituent {
    private:

    public:
        virtual std::unique_ptr<AtmosphereGridStorage> populate_storage(const Geometry& geometry) const = 0;
    };


    class UserDefined1DExtinctionConstituent : public Constituent {
    private:
        std::unique_ptr<sasktran2::grids::AltitudeGrid> m_altitude_grid;
        Eigen::MatrixXd m_extinction;

    public:
        UserDefined1DExtinctionConstituent(const Eigen::VectorXd& altitudes, const Eigen::MatrixXd& extinction);

        std::unique_ptr<AtmosphereGridStorage> populate_storage(const Geometry& geometry) const;

    };

}