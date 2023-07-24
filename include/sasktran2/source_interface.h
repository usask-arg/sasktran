#pragma once

#include <sasktran2/internal_common.h>
#include <sasktran2/raytracing.h>
#include <sasktran2/dual.h>
#include <sasktran2/config.h>
#include <sasktran2/atmosphere/atmosphere.h>

/** Base interface class that provides source term functionality to the Engine
 *
 */
template<int NSTOKES>
class SourceTermInterface {
protected:


public:
    virtual ~SourceTermInterface() {};

    virtual void initialize_config(const sasktran2::Config& config) {};

    /** Initializes any geometry information that is required for calculating the source term.  This method is called
     *  after the line of sight rays ar traced.
     *
     * @param los_rays The traced line of sight rays
     */
    virtual void initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay>& los_rays) {};

    /**
     *
     */
    virtual void initialize_atmosphere(const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) {};

    /** Triggers an internal calculation of the source term.  This method is called at the beginning of each
     *  'wavelength' calculation.
     *
     * @param wavelidx Index of the wavelength being calculated
     */
    virtual void calculate(int wavelidx, int threadidx) {};

    // TODO: Is Dual proper here? what about when the source term derivative is sparse? Maybe it isn't that important...
    // Should we be templated over NSTOKES in the interface?

    /** Calculates the integrated source term for a given layer.
     *
     * @param losidx Raw index pointing to the ray that was previously passed in initialize_geometry
     * @param layeridx Raw index pointing to the layer that was previosuly passed in initialize_geometry
     * @param layer The layer that we are integrating over
     * @param source The returned source term
     */
    virtual void integrated_source(int wavelidx,  int losidx, int layeridx, int threadidx, const sasktran2::raytracing::SphericalLayer& layer,
                                   const sasktran2::SparseODDualView& shell_od,
                                   sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source) const = 0;

    /** Calculates the source term at the end of the ray.  Common examples of this are ground scattering, ground emission,
     *  or the solar radiance if looking directly at the sun.
     *
     * @param wavelidx Raw index for the wavelength we are calculating
     * @param losidx Raw index pointing to the ray that was previously passed in initialize_geometry
     * @param source The returned source term
     */
    virtual void end_of_ray_source(int wavelidx, int losidx, int threadidx, sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source) const = 0;
};

