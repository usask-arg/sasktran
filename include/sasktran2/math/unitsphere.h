#pragma once


#include <sasktran2/internal_common.h>


namespace sasktran2::math {
    /** Base interface class for UnitSpheres inside Sasktran2
     *
     */
    class UnitSphere {
    private:
    public:
        virtual ~UnitSphere() {};

        virtual int num_points() const = 0;

        virtual Eigen::Vector3d get_quad_position(int index) const = 0;

        virtual double quadrature_weight(int i) const = 0;

        virtual void interpolate(const Eigen::Vector3d& direction,
                                 std::vector<std::pair<int, double>>& index_weights,
                                 int& num_interp
                                 ) const = 0;
    };

    /**
     *
     */
    class UnitSphereGround : public UnitSphere {
    private:
        std::unique_ptr<const UnitSphere> m_full_sphere;
        const Eigen::Vector3d m_location;

        std::vector<int> m_contributing_map;
        std::vector<int> m_reverse_contributing_map;
        std::vector<bool> m_is_full_sphere_looking_up;

    public:
        UnitSphereGround(std::unique_ptr<const UnitSphere>&& sphere, const Eigen::Vector3d location);

        int num_points() const override;

        Eigen::Vector3d get_quad_position(int index) const override;

        double quadrature_weight(int i) const override;

        void interpolate(const Eigen::Vector3d& direction,
                                 std::vector<std::pair<int, double>>& index_weights,
                                 int& num_interp
        ) const override;
    };

    /**
     * @class LebedevSphere
     *
     * @brief A class that represents a unit sphere with Lebedev quadrature points.
     *
     * The class provides a convenient way to access the cartesian coordinates and weights of the Lebedev quadrature points.
     * The cartesian coordinates are stored in the `m_xyz` member variable and the weights are stored in the `m_weights` member variable.
     */
    class LebedevSphere : public UnitSphere {
    private:
        const Eigen::MatrixXd m_xyz;
        const Eigen::VectorXd m_weights;
    public:
        /**
         * @brief Constructs a LebedevSphere with a specified number of points.
         *
         * The number of points must be one of the following: [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890].
         *
         * @param npoints The number of quadrature points to be used in constructing the sphere.
         */
        LebedevSphere(int npoints);

        /**
         *
         * @return The number of points in the unit sphere quadrature
         */
        int num_points() const override { return (int)m_weights.size(); }

        /**
         * Performs integration over the unit sphere for a set of values specified on the grid ponts.  The values
         * can be vector valued, i.e., N values for each grid point.  In this case the result is a vector of length
         * N.
         *
         * @tparam N Length of result
         * @param values Matrix of values to integrate over.  These are specified on the internal m_xyz quadrature points
         * @param result The integrated result of size N
         */
        template<int N>
        void integrate_on_grid(const Eigen::Matrix<double, N, -1>& values, Eigen::Vector<double, N>& result);

        /**
        * @brief Returns the quadrature position for a given index.
        *
        * @param index The index of the desired quadrature position.
        * @return Eigen::Vector3d The quadrature position as a 3-dimensional vector.
        */
        Eigen::Vector3d get_quad_position(int index) const;

        /** Returns the quadrature weight at index i
         *
         * @param i Index
         * @return Weight
         */
        double quadrature_weight(int i) const override { return m_weights[i]; }

        void interpolate(const Eigen::Vector3d& direction,
                         std::vector<std::pair<int, double>>& index_weights,
                         int& num_interp
        ) const override;
    };
}

