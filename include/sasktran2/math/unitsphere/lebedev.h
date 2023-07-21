#pragma once

#include <sasktran2/internal_common.h>


namespace sasktran2::math::unitsphere::lebedev {
    void get_lebedev_data(int npoints, Eigen::MatrixXd& result);
}
