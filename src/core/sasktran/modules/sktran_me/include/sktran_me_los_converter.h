#pragma once

namespace sktran_me {
    class GeometryConstructor {
    private:
        GEODETIC_INSTANT m_refpt;
        ReferencePointEstimator m_refpt_estimator;

        void add_limb_viewing_los(const sasktran2::Geometry1D& geometry,
                                  nxGeodetic& geodetic,
                                  const SKTRAN_LineOfSightEntry_V2* los,
                                  std::unique_ptr<sasktran2::viewinggeometry::ViewingGeometryContainer> &viewing_geo) const;

        void not_implemented_error(const std::string& type) const;

        void construct_viewing_geometry(const sasktran2::Geometry1D& geometry,
                                        const SKTRAN_LineOfSightArray_V21& linesofsight,
                                        nxGeodetic& geodetic,
                                        std::unique_ptr<sasktran2::viewinggeometry::ViewingGeometryContainer>& viewing_geo) const;

    public:
        GeometryConstructor() = default;
        ~GeometryConstructor() = default;

        void construct_from_config(const sasktran2::Config& config,
                                   const Eigen::VectorXd& altitude_grid,
                                   nxGeodetic& geodetic,
                                   const SKTRAN_LineOfSightArray_V21& linesofsight,
                                   std::unique_ptr<sasktran2::viewinggeometry::ViewingGeometryContainer>& viewing_geo,
                                   std::unique_ptr<sasktran2::Geometry1D>& geometry,
                                   const nxVector& sun_unit
                                   );


        const GEODETIC_INSTANT& reference_point() const { return m_refpt; }

    };
}