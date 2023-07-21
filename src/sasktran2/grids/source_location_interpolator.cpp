#include <sasktran2/grids.h>

namespace sasktran2::grids {

    SourceLocationInterpolator::SourceLocationInterpolator(AltitudeGrid &&altitude_grid)
    : m_altitude_grid(altitude_grid) {

    }

}