#pragma once

#include <sasktran2/math/trig.h>

namespace sasktran2::math {
    /** Calculates the angle between two vectors using the dot product.
     *
     * @param v1 First vector
     * @param v2 Second vector
     * @return Angle between the two vectors in degrees
     */
    inline double angle_degrees_between(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) {
        double cos_angle = (v1.dot(v2)) / (v1.norm() * v2.norm());
        if(cos_angle < -1) cos_angle = -1.0;
        if(cos_angle > 1) cos_angle = 1.0;

        return acosd(cos_angle);
    }

    /** Calculates the component of a vector that is perpindicular to another vector.
     *
     * @param v1 First vector
     * @param v2 Second vector
     * @return The components of v1 perpindicular to v2
     */
    inline Eigen::Vector3d component_perpindicular_to(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) {
        double zcomp;

        Eigen::Vector3d zaxis = v2 / v2.norm();
        zcomp = v1.dot(zaxis);

        Eigen::Vector3d xaxis = v1 - (zaxis * zcomp);

        return xaxis;
    }
}