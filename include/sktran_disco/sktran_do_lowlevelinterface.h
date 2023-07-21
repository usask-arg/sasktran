#pragma once
//#include "sktran_disco/sktran_do.h"
#include <vector>

// Low level interface to provide a C api to the model

namespace sasktran_disco_lowlevel {
    struct Atmosphere {
        double* od;                                 // [lyr, wavel]
        double* ssa;                                // [lyr, wavel]
        double* f;                                  // [lyr, wavel]
        double* a1;                                 // [nstr, lyr, wavel]
        double* a2;                                 // [nstr, lyr, wavel]
        double* a3;                                 // [nstr, lyr, wavel]
        double* a4;                                 // [nstr, lyr, wavel]
        double* b1;                                 // [nstr, lyr, wavel]
        double* b2;                                 // [nstr, lyr, wavel]

        double* ss_phase;                           // [stokes, nlos, lyr, wavel]

        // BRDF coefficients?
        double* albedo;                             // [wavel]

        double* layerboundaryaltitude;              // [lyr]

        double earthradius;
    };

    struct WeightingFunctions {
        double* d_od;                                 // [deriv, wavel]
        double* d_ssa;                                // [deriv, wavel]
        double* d_f;                                  // [deriv, wavel]
        double* d_a1;                                 // [deriv, nstr, wavel]
        double* d_a2;                                 // [deriv, nstr, wavel]
        double* d_a3;                                 // [deriv, nstr, wavel]
        double* d_a4;                                 // [deriv, nstr, wavel]
        double* d_b1;                                 // [deriv, nstr, wavel]
        double* d_b2;                                 // [deriv, nstr, wavel]

        double* d_ss_phase;                           // [deriv, stokes, nlos, wavel]

        double* d_albedo;                             // [deriv, wavel]

        int* d_layerindex;                              // [deriv], assuming the same for every wavel

        int numderiv;
    };

    struct Config {
        int nstr;
        int nwavel;
        int nlyr;
        int nstokes;
        int nthreads;

        bool useexactsinglescatter;
        bool usepseudospherical;

        int numazimuthexpansion;
    };

    struct ViewingGeometry {
        double* cos_vza;
        double cos_sza;
        double* saa;

        double* viewingaltitude;

        int nlos;
    };

    struct Output {
        double* radiance;           // [stokes, los, wavel] (preallocated)
        double* d_radiance;         // [deriv, stokes, los, wavel] (preallocated)
        // int debug;
    };

    // Exposed, untemplated main interface function, basically just validates the input and then dispatches to
    // a templated internal calculate
    void calculate(const Atmosphere* atmosphere, const Config* config,
                   const WeightingFunctions* weightingfunctions,
                   const ViewingGeometry* geometry, const Output* output);


    // Internal templated function
    template<int NSTOKES, int CNSTR=-1>
    void internal_calculate(const Atmosphere* atmosphere, const Config* config,
                            const WeightingFunctions* weightingfunctions,
                            const ViewingGeometry* geometry, const Output* output);

    // CPP Helper classe to construct the inputs to the C api
    class CPPApi {
    private:
        std::vector<double> m_od;                                 // [lyr, wavel]
        std::vector<double> m_ssa;                                // [lyr, wavel]
        std::vector<double> m_f;                                  // [lyr, wavel]
        std::vector<double> m_a1;                                 // [nstr, lyr, wavel]
        std::vector<double> m_a2;                                 // [nstr, lyr, wavel]
        std::vector<double> m_a3;                                 // [nstr, lyr, wavel]
        std::vector<double> m_a4;                                 // [nstr, lyr, wavel]
        std::vector<double> m_b1;                                 // [nstr, lyr, wavel]
        std::vector<double> m_b2;                                 // [nstr, lyr, wavel]

        // BRDF coefficients?
        std::vector<double> m_albedo;                             // [wavel]

        std::vector<double> m_layerboundaryaltitude;              // [lyr]

        // Derivative options
        std::vector<double> m_d_ssa;                                // [deriv, wavel]
        std::vector<double> m_d_od;                                 // [deriv, wavel]
        std::vector<double> m_d_f;                                  // [deriv, wavel]
        std::vector<double> m_d_a1;                                 // [deriv, nstr, wavel]
        std::vector<double> m_d_a2;                                 // [deriv, nstr, wavel]
        std::vector<double> m_d_a3;                                 // [deriv, nstr, wavel]
        std::vector<double> m_d_a4;                                 // [deriv, nstr, wavel]
        std::vector<double> m_d_b1;                                 // [deriv, nstr, wavel]
        std::vector<double> m_d_b2;                                 // [deriv, nstr, wavel]
        std::vector<double> m_d_ss_phase;                           // [nstokes, deriv, nlos, nwavel]

        std::vector<double> m_d_albedo;                             // [deriv, wavel]

        std::vector<int> m_d_layerindex;                              // [deriv], assuming the same for every wavel

        double m_earthradius;

        int m_nstr;
        int m_nlyr;
        int m_nwavel;
        int m_nstokes;
        int m_nthreads;
        int m_nderiv;

        std::vector<double> m_cos_vza;                            // [los]
        double m_cos_sza;
        std::vector<double> m_saa;                                // [los] (relative azimuth angle in radians)
        std::vector<double> m_viewingaltitude;                    // [los]

        int m_nlos;

        std::vector<double> m_radiance;                           // [stokes, los, wavel]
        std::vector<double> m_d_radiance;                         // [deriv, stokes, los, wavel]

        int linear_index(int lyr, int wav) {
            return lyr + wav * m_nlyr;
        }
        int linear_index(int str, int lyr, int wav) {
            return str + lyr * m_nstr + wav * m_nstr * m_nlyr;
        }

        int d_linear_index(int deriv, int wav) {
            return wav * m_nderiv + deriv;
        }

    public:
        CPPApi(int nstr, int nlyr, int nwavel, int nstokes, int nlos, int nderiv) {
            this->m_nstr = nstr;
            this->m_nlyr = nlyr;
            this->m_nwavel = nwavel;
            this->m_nstokes = nstokes;
            this->m_nlos = nlos;
            this->m_nderiv = nderiv;

            m_od.resize(nlyr * nwavel, 0);
            m_ssa.resize(nlyr * nwavel, 0);
            m_f.resize(nlyr * nwavel, 0);
            m_a1.resize(nlyr * nwavel * nstr, 0);
            m_a2.resize(nlyr * nwavel * nstr, 0);
            m_a3.resize(nlyr * nwavel * nstr, 0);
            m_a4.resize(nlyr * nwavel * nstr, 0);
            m_b1.resize(nlyr * nwavel * nstr, 0);
            m_b2.resize(nlyr * nwavel * nstr, 0);
            m_albedo.resize(nwavel, 0);
            m_layerboundaryaltitude.resize(nlyr, 0);

            m_d_ssa.resize(nderiv * nwavel, 0);
            m_d_od.resize(nderiv * nwavel, 0);
            m_d_f.resize(nderiv * nwavel, 0);
            m_d_a1.resize(nderiv * nstr * nwavel, 0);
            m_d_a2.resize(nderiv * nstr * nwavel, 0);
            m_d_a3.resize(nderiv * nstr * nwavel, 0);
            m_d_a4.resize(nderiv * nstr * nwavel, 0);
            m_d_b1.resize(nderiv * nstr * nwavel, 0);
            m_d_b2.resize(nderiv * nstr * nwavel, 0);
            m_d_ss_phase.resize(nderiv * nstokes * nlos * nwavel, 0);

            m_d_layerindex.resize(nderiv);

            m_d_albedo.resize(nderiv * nwavel, 0);

            m_radiance.resize(nstokes * nlos * nwavel);
            m_d_radiance.resize(nderiv * nstokes * nlos * nwavel);

            m_cos_vza.resize(nlos);
            m_saa.resize(nlos);
            m_viewingaltitude.resize(nlos, -1);

            // Default configuration options
            m_nthreads = -1;
        }

        void initialize_c_api(Atmosphere* atmosphere, Config* config, ViewingGeometry* geometry, WeightingFunctions* weightingfunctions, Output* output) {
            // Atmosphere construction
            atmosphere->a1 = &m_a1[0];
            atmosphere->a2 = &m_a2[0];
            atmosphere->a3 = &m_a3[0];
            atmosphere->a4 = &m_a4[0];
            atmosphere->b1 = &m_b1[0];
            atmosphere->b2 = &m_b2[0];
            atmosphere->albedo = &m_albedo[0];
            atmosphere->od = &m_od[0];
            atmosphere->ssa = &m_ssa[0];
            atmosphere->f = &m_f[0];
            atmosphere->layerboundaryaltitude = &m_layerboundaryaltitude[0];
            atmosphere->earthradius = m_earthradius;

            atmosphere->ss_phase = nullptr;

            // Config construction
            config->nlyr = m_nlyr;
            config->nstokes = m_nstokes;
            config->nstr = m_nstr;
            config->nwavel = m_nwavel;
            config->nthreads = m_nthreads;

            config->useexactsinglescatter = false;
            config->numazimuthexpansion = 0;
            config->usepseudospherical = true;

            // Viewing geometry construction
            geometry->cos_sza = m_cos_sza;
            geometry->cos_vza = &m_cos_vza[0];
            geometry->saa = &m_saa[0];
            geometry->viewingaltitude = &m_viewingaltitude[0];
            geometry->nlos = m_nlos;

            output->radiance = &m_radiance[0];
            if (m_nderiv > 0) {
                output->d_radiance = &m_d_radiance[0];
                // Weighting function input construction
                weightingfunctions->d_a1 = &m_d_a1[0];
                weightingfunctions->d_a2 = &m_d_a2[0];
                weightingfunctions->d_a3 = &m_d_a3[0];
                weightingfunctions->d_a4 = &m_d_a4[0];
                weightingfunctions->d_b1 = &m_d_b1[0];
                weightingfunctions->d_b2 = &m_d_b2[0];
                weightingfunctions->d_ss_phase = &m_d_ss_phase[0];

                weightingfunctions->d_albedo = &m_d_albedo[0];
                weightingfunctions->d_od = &m_d_od[0];
                weightingfunctions->d_ssa = &m_d_ssa[0];
                weightingfunctions->d_f = &m_d_f[0];
                weightingfunctions->d_layerindex = &m_d_layerindex[0];
            }
            else {
                output->d_radiance = nullptr;
            }
            weightingfunctions->numderiv = m_nderiv;
        }

        double& od(int lyr, int wav) {
            return m_od[linear_index(lyr, wav)];
        }

        double& d_od(int deriv, int wav) {
            return m_d_od[d_linear_index(deriv, wav)];
        }

        int& d_layerindex(int deriv) {
            return m_d_layerindex[deriv];
        }

        double& ssa(int lyr, int wav) {
            return m_ssa[linear_index(lyr, wav)];
        }

        double& d_ssa(int deriv, int wav) {
            return m_d_ssa[d_linear_index(deriv, wav)];
        }

        double& f(int lyr, int wav) {
            return m_f[linear_index(lyr, wav)];
        }

        double& albedo(int wav) {
            return m_albedo[wav];
        }

        double& a1(int str, int lyr, int wav) {
            return m_a1[linear_index(str, lyr, wav)];
        }

        double& a2(int str, int lyr, int wav) {
            return m_a2[linear_index(str, lyr, wav)];
        }

        double& a3(int str, int lyr, int wav) {
            return m_a3[linear_index(str, lyr, wav)];
        }

        double& a4(int str, int lyr, int wav) {
            return m_a4[linear_index(str, lyr, wav)];
        }

        double& b1(int str, int lyr, int wav) {
            return m_b1[linear_index(str, lyr, wav)];
        }

        double& b2(int str, int lyr, int wav) {
            return m_b2[linear_index(str, lyr, wav)];
        }

        double& layerboundaryaltitude(int lyr) {
            return m_layerboundaryaltitude[lyr];
        }

        double& cos_sza() {
            return m_cos_sza;
        }

        double& cos_vza(int los) {
            return m_cos_vza[los];
        }

        double& saa(int los) {
            return m_saa[los];
        }

        double& earth_radius() {
            return m_earthradius;
        }


    };

}

