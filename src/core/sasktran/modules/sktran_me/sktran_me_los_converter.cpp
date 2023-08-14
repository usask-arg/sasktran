#include "include/sktran_me_internals.h"

namespace sktran_me {
    void GeometryConstructor::construct_viewing_geometry(const sasktran2::Geometry1D& geometry,
                                                         const SKTRAN_LineOfSightArray_V21 &linesofsight,
                                                         nxGeodetic& geodetic,
                                                         std::unique_ptr<sasktran2::viewinggeometry::ViewingGeometryContainer> &viewing_geo) const {
        const SKTRAN_LineOfSightEntry_V2* entry;
        for(int i = 0; i < linesofsight.NumRays(); ++i) {
            linesofsight.GetRay(i, &entry);

            switch (entry->DefaultViewingType(geodetic, geometry.altitude_grid().grid()(Eigen::last))) {
                case SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_UNDEFINED:
                    not_implemented_error("SKTRAN_VIEWING_TYPE_UNDEFINED");
                    break;
                case SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_INSPACE_LOOKINGATSPACE:
                    not_implemented_error("SKTRAN_VIEWING_TYPE_INSPACE_LOOKINGATSPACE");
                    break;
                case SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_INSPACE_LOOKINGATLIMB:
                    add_limb_viewing_los(geometry, geodetic, entry, viewing_geo);
                    break;
                case SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_INSPACE_LOOKINGATNADIR:
                    not_implemented_error("SKTRAN_VIEWING_TYPE_INSPACE_LOOKINGATNADIR");
                    break;
                case SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_INATMOS_LOOKINGATSPACE:
                    not_implemented_error("SKTRAN_VIEWING_TYPE_INATMOS_LOOKINGATSPACE");
                    break;
                case SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_INATMOS_LOOKINGATLIMB:
                    not_implemented_error("SKTRAN_VIEWING_TYPE_INATMOS_LOOKINGATLIMB");
                    break;
                case SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_INATMOS_LOOKINGATNADIR:
                    not_implemented_error("SKTRAN_VIEWING_TYPE_INATMOS_LOOKINGATNADIR");
                    break;
                case SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_NEARGROUND:
                    not_implemented_error("SKTRAN_VIEWING_TYPE_NEARGROUND");
                    break;
                case SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_ENDOFLIST:
                    not_implemented_error("SKTRAN_VIEWING_TYPE_ENDOFLIST");
                    break;
            }

        }

    }

    void GeometryConstructor::add_limb_viewing_los(const sasktran2::Geometry1D& geometry,
                                                   nxGeodetic &geodetic,
                                                   const SKTRAN_LineOfSightEntry_V2* los,
                                                   std::unique_ptr<sasktran2::viewinggeometry::ViewingGeometryContainer> &viewing_geo) const {
        // Get the tangent point of the LOS
        geodetic.FromTangentPointLocation(los->Observer(), los->Look());
        double tangent_altitude = geodetic.Height();

        Eigen::Vector3d tangent_point(geodetic.Location().X(), geodetic.Location().Y(), geodetic.Location().Z());
        Eigen::Vector3d look(los->Look().X(), los->Look().Y(), los->Look().Z());

        sasktran2::Location loc;
        nxVector west, south, up;
        geodetic.GetGeodeticWestSouthUp(&west, &south, &up);

        loc.position = Eigen::Vector3d(up.X(), up.Y(), up.Z());

        double csz;
        double rel_az;

        sasktran2::raytracing::calculate_csz_saz(m_geographic_sun, loc, look, csz, rel_az);



        viewing_geo->observer_rays().emplace_back(
                std::make_unique<sasktran2::viewinggeometry::TangentAltitudeSolar>(tangent_altitude,
                                                                              rel_az,
                                                                              geometry.altitude_grid().grid()(Eigen::last) + 1000,
                                                                              csz)
        );


    }

    void GeometryConstructor::not_implemented_error(const std::string& type) const {
        BOOST_LOG_TRIVIAL(error) << "Adding a line of sight type: " << type << " that is not implemented by sktran_me";
    }

    void GeometryConstructor::construct_from_config(const sasktran2::Config &config,
                                                    const Eigen::VectorXd& altitude_grid,
                                                    nxGeodetic &geodetic,
                                                    const SKTRAN_LineOfSightArray_V21 &linesofsight,
                                                    std::unique_ptr<sasktran2::viewinggeometry::ViewingGeometryContainer> &viewing_geo,
                                                    std::unique_ptr<sasktran2::Geometry1D> &geometry,
                                                    const nxVector& sun_unit
                                                    ) {
        m_geographic_sun = Eigen::Vector3d(sun_unit.X(), sun_unit.Y(), sun_unit.Z());

        m_refpt_estimator.estimate_reference_point(linesofsight, m_refpt);

        // Use the geodetic to construct the local coordinate system
        geodetic.FromGeodetic(m_refpt.latitude, m_refpt.longitude, m_refpt.heightm);

        double earth_radius = geodetic.Location().Magnitude();
        nxVector west, south, up;
        geodetic.GetGeodeticWestSouthUp(&west, &south, &up);

        // Get the cos_sza at the reference point, set SAA=0 mostly for stokes basis reasons
        // This means that rel_az is contained entirely in the LOS azimuth
        double cos_sza = sun_unit & up;
        sasktran2::Coordinates coords(cos_sza, 0, earth_radius);

        Eigen::VectorXd grid_values = altitude_grid;

        sasktran2::grids::AltitudeGrid grid = sasktran2::grids::AltitudeGrid(std::move(grid_values), sasktran2::grids::gridspacing::constant,
                                                                             sasktran2::grids::outofbounds::extend,
                                                                             sasktran2::grids::interpolation::linear);

        geometry = std::make_unique<sasktran2::Geometry1D>(std::move(coords), std::move(grid));

        viewing_geo = std::make_unique<sasktran2::viewinggeometry::ViewingGeometryContainer>();
        construct_viewing_geometry(*geometry, linesofsight, geodetic, viewing_geo);
    }
}