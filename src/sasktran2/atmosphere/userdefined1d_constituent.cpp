#include <sasktran2/atmosphere/constituent.h>

namespace sasktran2::atmosphere {
	UserDefined1DExtinctionConstituent::UserDefined1DExtinctionConstituent(const Eigen::VectorXd& altitudes, const Eigen::MatrixXd& extinction) {
		// Construct an Internal Grid from the measurements

		// Copy the altitude grid to a new vector
		Eigen::VectorXd alt_copy = altitudes;

		// Make sure the sizes are good
		assert(altitudes.size() > 1);
		assert(altitudes.size() == extinction.rows());

		// TODO: Determine if the grid is uniformly spaced

		// and move it into the altitude grid
		m_altitude_grid = std::make_unique<sasktran2::grids::AltitudeGrid>(std::move(alt_copy),
			sasktran2::grids::gridspacing::variable,
			sasktran2::grids::outofbounds::setzero,
			sasktran2::grids::interpolation::linear);
		
		m_extinction = extinction;
	}

	std::unique_ptr<AtmosphereGridStorage> UserDefined1DExtinctionConstituent::populate_storage(const Geometry& geometry) const {
		Eigen::MatrixXd extinction;

		for (int i = 0; i < geometry.size(); ++i) {

		}

		return std::make_unique<AtmosphereGridStorageAbsorber>(extinction);
	}
}