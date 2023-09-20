#pragma once
#include "sktran_disco/sktran_do.h"

/* Standalone functions for postprocessing of discrete ordinates solutions
 *
 */

namespace sasktran_disco::postprocessing {
    // Layer multipliers

    /** Calculates the homogeneous solution multiplier, h_plus, sampled at a relative location
     *  layer_fraction within the layer.
     *
     * @param thickness Layer optical depth
     * @param eigval Eigenvalue corresponding to this homogenous solution
     * @param layer_relative_location Relative location in optical depth coordinates. 0 for top of the layer, 1 for bottom of the layer.
     */
     template<int CNSTR=-1>
    void h_plus_sampled(const sasktran_disco::LayerDual<double>& thickness,
                        const sasktran_disco::VectorLayerDual<double, CNSTR>& eigval,
                        sasktran_disco::SolutionIndex solution_index,
                        double layer_relative_location,
                        sasktran_disco::LayerDual<double>& h_plus
                        );

    /** Calculates the homogeneous solution multiplier, h_minus, sampled at a relative location
     *  layer_fraction within the layer.
     *
     * @param thickness Layer optical depth
     * @param eigval Eigenvalue corresponding to this homogenous solution
     * @param layer_relative_location Relative location in optical depth coordinates. 0 for top of the layer, 1 for bottom of the layer.
     */
    template<int CNSTR=-1>
    void h_minus_sampled(const sasktran_disco::LayerDual<double>& thickness,
                         const sasktran_disco::VectorLayerDual<double, CNSTR>& eigval,
                         sasktran_disco::SolutionIndex solution_index,
                         double layer_relative_location,
                         sasktran_disco::LayerDual<double>& h_minus
    );

    /** Calculates the green's function solution multipier, d_minus, sampled at a relative location layer_fraction
     *  within the layer
     *
     * @param thickness
     * @param eigval
     * @param solution_index
     * @param layer_relative_location
     * @param transmission
     * @param average_secant
     * @param layerderivstart
     * @param d_minus
     */
    template<int CNSTR=-1>
    void d_minus_sampled(const sasktran_disco::LayerDual<double>& thickness,
                         const sasktran_disco::VectorLayerDual<double, CNSTR>& eigval,
                         sasktran_disco::SolutionIndex solution_index,
                         double layer_relative_location,
                         const sasktran_disco::Dual<double>& transmission,
                         const sasktran_disco::Dual<double>& average_secant,
                         int layerderivstart,
                         sasktran_disco::Dual<double>& d_minus
                         );

    /** Calculates the green's function solution multipier, d_plus, sampled at a relative location layer_fraction
     *  within the layer
     *
     * @param thickness
     * @param eigval
     * @param solution_index
     * @param layer_relative_location
     * @param transmission
     * @param average_secant
     * @param layerderivstart
     * @param d_plus
     */
    template<int CNSTR=-1>
    void d_plus_sampled(const sasktran_disco::LayerDual<double>& thickness,
                        const sasktran_disco::VectorLayerDual<double, CNSTR>& eigval,
                        sasktran_disco::SolutionIndex solution_index,
                        double layer_relative_location,
                        const sasktran_disco::Dual<double>& transmission,
                        const sasktran_disco::Dual<double>& average_secant,
                        int layerderivstart,
                        sasktran_disco::Dual<double>& d_plus
    );
}