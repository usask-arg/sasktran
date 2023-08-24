#pragma once


#include <sasktran2/internal_common.h>

namespace sasktran2 {
    class Coordinates;
    struct Location;
}

namespace sasktran2::grids {

    /** Interpolation mode between layer boundaries.  Note that currently 'constant' is not fully implemented.
     *
     * 'linear' performs linear interpolation between layer boundaries
     *
     * 'constant' assumes the the layer consists of optical properties equal to 0.5*lower + 0.5*upper
     *
     */
    enum interpolation {
        shell,
        linear
    };

    /** Keeps track of the spacing between grid points within a grid.
     *
     * 'constant' indicates that all grid points are evenly spaced.  If this is the case interpolation speed can be
     * greatly improved
     *
     * 'variable' indicates that the spacing may be variable between grid points.  If this is the case then usually a
     * binary search has to be performed to find interpolation indicies.
     *
     */
    enum gridspacing {
        constant,
        variable
    };

    /** The policy for points that fall outside of the grid boundaries.
     *
     * 'extend' Extends the last value beyond the grid points
     *
     * 'setzero' Sets the values to be 0 outside of the grid.
     *
     */
    enum outofbounds {
        extend,
        setzero
    };

    /** Generic class that implements all of the standard functionality for a one-dimensional grid to interpolate
     *  over.
     */
    class Grid {
    private:
        const interpolation m_interp_method;
        const gridspacing m_grid_spacing;
        const outofbounds m_out_of_bounds_mode;
        const Eigen::VectorXd m_grid_values;

        double m_dx; /**< Set when m_grid_spacing is constant, set to the distance between grid points */
        double m_x0; /**< Set when m_grid spacing is constant, set to the first value of the grid */

        /** Internal method to interpolate when the grid spacing is constant.
         *
         * @param x Value to interpolate to
         * @param index Returned indexes
         * @param weight Returned weights
         * @param num_contributing Number of constributing interpolation points
         */
        void interpolate_constant_spacing(
                double x,
                std::array<int, 2> & index,
                std::array<double, 2>& weight,
                int& num_contributing
                ) const;

        /** Internal method to interpolate when the grid spacing is variable.
         *
         * @param x Value to interpolate to
         * @param index Returned indexes
         * @param weight Returned weights
         * @param num_contributing Number of constributing interpolation points
         */
        void interpolate_varying_spacing(
                double x,
                std::array<int, 2> & index,
                std::array<double, 2>& weight,
                int& num_contributing
                ) const;

    public:
        /** Initializes the one-dimensional grid
         *
         * @param grid_values The vector of points defining the grid, it is moved into the class upon creation
         * @param spacing The spacing mode, either constant or varying.  Set to constant if the grid is uniformly spaced
         * @param out_of_bounds_mode Out of bounds mode, either extend or setzero.  Extend will always use the last/first value for out of bounds values, setzero will set it to 0.
         * @param interp Interpolation mode, either linear or constant.
         */
        Grid(Eigen::VectorXd&& grid_values,
             gridspacing spacing,
             outofbounds out_of_bounds_mode,
             interpolation interp
             );

        /** Calculates the interpolation indicies and weights.  The number of contributing indices is stored in
         * num_contributing, the maximum number of contributing points is 2.  Note that it may be 0 if the interpolated
         * value is always 0.
         *
         * @param x Value to interpolate to
         * @param index The indicies in the grid
         * @param weight The corresponding weights to the indicies
         * @param num_contributing The number of contributing points to the interpolation
         */
        void calculate_interpolation_weights(double x,
                                             std::array<int, 2> & index,
                                             std::array<double, 2>& weight,
                                             int& num_contributing
                                             ) const ;

        /**
         *
         * @return  The internal grid
         */
        const Eigen::VectorXd& grid() const { return m_grid_values; }

    };

    /** An AltitudeGrid is essentially the same as a Grid except we give it a new name to be specific
     *
     */
    class AltitudeGrid : public Grid {
    public:
        /** Initializes the one-dimensional grid
        *
        * @param grid_values The vector of points defining the grid, it is moved into the class upon creation
        * @param spacing The spacing mode, either constant or varying.  Set to constant if the grid is uniformly spaced
        * @param out_of_bounds_mode Out of bounds mode, either extend or setzero.  Extend will always use the last/first value for out of bounds values, setzero will set it to 0.
        * @param interp Interpolation mode, either linear or constant.
        */
        AltitudeGrid(Eigen::VectorXd&& grid_values,
        gridspacing spacing,
        outofbounds out_of_bounds_mode,
        interpolation interp
        ) : Grid(std::forward<Eigen::VectorXd&&>(grid_values), spacing, out_of_bounds_mode, interp) {}
    };


    /** A convenience class which combines multiple grids to provide interpolation scheme for source functions.
     *  Note this includes the location component of the source function but not the directional component.
     *
     *  Examples are, in a spherically symmetric atmosphere the source function depends only on solar zenith angle
     *  and altitude.
     *
     *  In a plane parallel horizontally homogenous atmosphere, the source function depends only on altitude.
     *
     *  All source functions will have an altitude component, so we include the altitude dependence in the base class
     *
     */
    class SourceLocationInterpolator {
    protected:
        const AltitudeGrid m_altitude_grid;
    public:
        SourceLocationInterpolator(AltitudeGrid&& altitude_grid);

        virtual ~SourceLocationInterpolator() {};

        virtual int num_interior_points() const = 0;

        virtual int num_ground_points() const = 0;

        virtual Eigen::Vector3d grid_location(const sasktran2::Coordinates& coords, int location_index) const = 0;
        virtual Eigen::Vector3d ground_location(const sasktran2::Coordinates& coords, int ground_index) const = 0;
        virtual void interior_interpolation_weights(const sasktran2::Coordinates& coords, const sasktran2::Location& location, std::vector<std::pair<int, double>>& weights, int& num_interp) = 0;
        virtual void ground_interpolation_weights(const sasktran2::Coordinates& coords, const sasktran2::Location& location, std::vector<std::pair<int, double>>& weights, int& num_interp) const = 0;
    };


    class AltitudeSZASourceLocationInterpolator : public SourceLocationInterpolator {
    private:
        const Grid m_cos_sza_grid;

        int interior_linear_index(int alt_index, int sza_index);
        int ground_linear_index(int sza_index) const;
    public:
        AltitudeSZASourceLocationInterpolator(AltitudeGrid&& altitude_grid,
                                              Grid&& cos_sza_grid);

        int num_interior_points() const override;
        int num_ground_points() const override;

        Eigen::Vector3d grid_location(const sasktran2::Coordinates& coords, int location_index) const override;
        Eigen::Vector3d ground_location(const sasktran2::Coordinates& coords, int ground_index) const override;

        void interior_interpolation_weights(const sasktran2::Coordinates& coords, const sasktran2::Location& location, std::vector<std::pair<int, double>>& weights, int& num_interp) override;
        void ground_interpolation_weights(const sasktran2::Coordinates& coords, const sasktran2::Location& location, std::vector<std::pair<int, double>>& weights, int& num_interp) const override;
    };


}



