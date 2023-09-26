#pragma once

#include <sasktran2/config.h>
#include <sasktran2/atmosphere/atmosphere.h>
#include <sasktran2/viewinggeometry.h>
#include <sasktran2/raytracing.h>
#include <sasktran2/solartransmission.h>
#include <sasktran2/source_integrator.h>
#include <sasktran2/output.h>

/** Internal essentially void class that the main Sasktran2 inherits from.  Necessary to remove the NSTOKES template
 *  to interface with SWIG.
 *
 */
class Sasktran2Interface {
public:
    virtual ~Sasktran2Interface() {}
};

/** The main object for SASKTRAN2 that performs the radiative transfer calculation.  The model is initialized with
 *  a Config object that contains user settings, a Geometry object to specify grids and coordinate information,
 *  and a ViewingGeometryContainer object that contains information on the viewing geometry.
 *
 *  After construction, the model may be called using an Atmosphere object.
 *
 * @tparam NSTOKES Number of stokes parameters, either 1 or 3
 */
template <int NSTOKES>
class Sasktran2 : public Sasktran2Interface {
private:
    const sasktran2::Config& m_config; /**< Internal reference to the config */
    const sasktran2::viewinggeometry::ViewingGeometryContainer& m_viewing_geometry; /**< Internal reference to the viewing geometry */
    const sasktran2::Geometry1D* m_geometry; /**< Internal reference to the model geometry */

    std::unique_ptr<const sasktran2::raytracing::RayTracerBase> m_raytracer; /**< Ray tracer that is internally constructed */

    std::vector<sasktran2::raytracing::TracedRay> m_traced_rays; /**< Traced observer rays */
    std::unique_ptr<sasktran2::SourceIntegrator<NSTOKES>> m_source_integrator; /** integrator for the source terms */

    std::vector<Eigen::SparseMatrix<double, Eigen::RowMajor>> m_traced_ray_od_matrix; /**< Vector of matrices A such that A * atmosphere_extinction = OD for each layer in that ray */

    // We split source terms by characterization
    std::vector<std::unique_ptr<SourceTermInterface<NSTOKES>>> m_source_terms; /**< All of the source terms */
    std::vector<SourceTermInterface<NSTOKES>*> m_los_source_terms; /**< Source terms that are integrated over the line of sight */
    std::vector<SourceTermInterface<NSTOKES>*> m_thermal_source; /**< Thermal source terms (planck, photochemical) */

    /** Internal method to calculate all terms inside the engine that are only geometry dependent
     */
    void calculate_geometry();

    /** Internal method to construct the internal ray tracer
     */
    void construct_raytracer();

    /** Internal method to construct all of the source terms from the config
     */
    void construct_source_terms();

    /** Internal method to construct the source term integrator from the config
     *
     */
    void construct_integrator();

    /** Internal method to validate that the input atmosphere is in the correct format
     *  and contains all of the necessary information.
     */
    void validate_input_atmosphere(const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) const;

public:
    /** Constructs the model
     *
     * @param config user configuration options
     * @param geometry geometry and coordinate information
     * @param viewing_rays information on the viewing geometry
     */
    Sasktran2(const sasktran2::Config& config,
              const sasktran2::Geometry1D* geometry,
              const sasktran2::viewinggeometry::ViewingGeometryContainer& viewing_rays) :
              m_config(config),
              m_viewing_geometry(viewing_rays),
              m_geometry(geometry)
              {
                  // First create the ray tracer
                  construct_raytracer();

                  // Then the integrator
                  construct_integrator();

                  // Then all of our source terms
                  construct_source_terms();

                  // And finally create all of the geometry only information
                  calculate_geometry();
              };

    virtual ~Sasktran2() {}

    /** Calculates the radiance for a given Atmosphere object.  May be called any number of times with different
     *  Atmosphere objects.  Results are stored in output.
     *
     * @param atmosphere
     * @param output
     */
    void calculate_radiance(const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere, sasktran2::Output<NSTOKES>& output) const;
};