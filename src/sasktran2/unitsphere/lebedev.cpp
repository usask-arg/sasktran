#include <sasktran2/math/unitsphere.h>
#include <sasktran2/math/unitsphere/lebedev.h>

namespace sasktran2::math {
    LebedevSphere::LebedevSphere(int npoints) {
        Eigen::MatrixXd result;
        unitsphere::lebedev::get_lebedev_data(npoints, result);

        Eigen::AngleAxis<double> rot(0.1, Eigen::Vector3d(0, 1, 0));

        const_cast<Eigen::MatrixXd&>(m_xyz) = rot.matrix() * result(Eigen::seq(0, 2), Eigen::all);
        const_cast<Eigen::VectorXd&>(m_weights) = result(3, Eigen::all);
    }

    template<int N>
    void LebedevSphere::integrate_on_grid(const Eigen::Matrix<double, N, -1>& values, Eigen::Vector<double, N>& result) {
        result.setZero();
        for (int i = 0; i < m_weights.size(); i++) {
            result += m_weights(i) * values.col(i);
        }
        result *= 4*EIGEN_PI;
    }

    Eigen::Vector3d LebedevSphere::get_quad_position(int index) const {
        return Eigen::Vector3d(m_xyz(0, index), m_xyz(1, index), m_xyz(2, index));
    }

    void LebedevSphere::interpolate(const Eigen::Vector3d& direction,
                                    std::vector<std::pair<int, double>> &index_weights, int &num_interp) const {
        num_interp = 3;
        // Set storage
        index_weights.resize(num_interp);

        // Store distances in the weights at first to save memory, initialize to a large value
        for(int i = 0; i < num_interp; ++i) {
            index_weights[i].second = 9999;
        }


        // TODO: other interpolation schemes?
        for(int i = 0; i < m_weights.size(); ++i) {
            double sqdist = (m_xyz(Eigen::all, i) - direction).squaredNorm();



            // Insert sorted
            for(int j = 0; j < num_interp; ++j) {
                if(sqdist < index_weights[j].second) {
                    // Have to bump every other weight upwards
                    auto temp = index_weights[j];
                    for(int k = j+1; k < num_interp; ++k) {
                        std::swap(temp, index_weights[k]);
                    }
                    // Then assign it
                    index_weights[j].first = i;
                    index_weights[j].second = sqdist;

                    sqdist = 99999;
                }

            }
        }


        // Now we have the correct indices and all square distances, assign the interpolation weights
        // Start by calculating the total inverse
        double total_distance = 0;
        for(int i = 0; i < num_interp; ++i) {
            // first check for a special case, if any of the distances are ~0 we just want to return that one back
            if(index_weights[i].second < 1e-8) {
                for(int j = 0; j < num_interp; ++j) {
                    index_weights[j].second = 0;
                }
                index_weights[i].first = i;
                index_weights[i].second = 1;
                return;
            }

            total_distance += 1 / sqrt(index_weights[i].second);
        }

        // Else we weight each one by inverse distance
        for(int i = 0; i < num_interp; ++i) {
            index_weights[i].second = (1 / sqrt(index_weights[i].second)) / total_distance;
        }
        /*
        for(int i = 0; i < 3; ++i) {
            BOOST_LOG_TRIVIAL(warning) << "dir:" << m_xyz(Eigen::all, index_weights[i].first).transpose() << " sqdist: " << index_weights[i].second;

        }
         */
    }

    template void LebedevSphere::integrate_on_grid<1>(const Eigen::Matrix<double, 1, -1>& values,
                                                   Eigen::Vector<double, 1>& result);
    template void LebedevSphere::integrate_on_grid<3>(const Eigen::Matrix<double, 3, -1>& values,
                                                   Eigen::Vector<double, 3>& result);

}