#pragma once
#include "modules/sktran_do_deprecated/include/sktran_do.h"
#include "omp.h"

namespace sktran_do_detail {
    // Temporaries needed during PostProcessing, create one for each layer even though some of it could be
    // shared between layers in theory
    template <int NSTOKES, int CNSTR=-1>
    struct PostProcessingCache {
        VectorLayerDual<double> dual_lpsum_plus, dual_lpsum_minus;
        VectorLayerDual<double> Y_plus, Y_minus;
        Radiance<NSTOKES> V, Q, J;

        sktran_do_detail::InhomogeneousSourceHolder<NSTOKES> Qtemp, temp;

        LayerDual<double> hp, hm;

        Dual<double> Dm, Dp, Eform;

        void resize(uint NSTR, LayerIndex p, uint numlayerderiv, uint layerstart, uint numtotalderiv) {
            // Vector Layer Duals
            dual_lpsum_plus.resize(NSTR/2 * NSTOKES * NSTOKES, numlayerderiv, p, layerstart);
            dual_lpsum_minus.resize(NSTR/2 * NSTOKES * NSTOKES, numlayerderiv, p, layerstart);
            Y_plus.resize(NSTR/2 * NSTOKES * NSTOKES, numlayerderiv, p, layerstart);
            Y_minus.resize(NSTR/2 * NSTOKES * NSTOKES, numlayerderiv, p, layerstart);

            // Radiance containers
            V.resize(numtotalderiv, false);
            Q.resize(numtotalderiv, false);
            J.resize(numtotalderiv, false);

            // Source holders
            Qtemp.resize(NSTR);
            temp.resize(NSTR);

            // Layer Duals
            hp.resize(numlayerderiv);
            hm.resize(numlayerderiv);

            // Full Duals
            Dm.resize(numtotalderiv);
            Dp.resize(numtotalderiv);
            Eform.resize(numtotalderiv);
        }
    };

    // Temporaries that are used inside every OpticalLayer
    template <int NSTOKES, int CNSTR=-1>
    struct LayerCache {
        LayerDual<double>						    dual_thickness;
        Dual<double>								average_secant;
        Dual<double>								dual_bt_floor;
        Dual<double>								dual_bt_ceiling;

        TripleProductDerivativeHolder<NSTOKES> triple_product_holder;
        LPTripleProduct<NSTOKES> triple_product;

        LayerCache(uint NSTR) : triple_product(NSTR), triple_product_holder(NSTR) {
        }
    };

    // Temporaries needed for the RTE solution
    template <int NSTOKES, int CNSTR>
    struct RTEMemoryCache {
        using HomogType = typename std::conditional<NSTOKES != 5, double, std::complex<double>>::type;

        using Matrix = typename std::conditional<CNSTR != -1, Eigen::Matrix<double, CNSTR/2 * NSTOKES, CNSTR/2 * NSTOKES>, Eigen::MatrixXd>::type;
        using MatrixH = typename std::conditional<CNSTR != -1, Eigen::Matrix<HomogType, CNSTR/2 * NSTOKES, CNSTR/2 * NSTOKES>, Eigen::MatrixXd>::type;
        using MatrixView = typename Eigen::Map<Matrix>;
        using MatrixViewH = typename Eigen::Map<MatrixH>;

        using Vector = typename std::conditional<CNSTR != -1, Eigen::Vector<double, CNSTR/2 * NSTOKES>, Eigen::VectorXd>::type;
        using VectorH = typename std::conditional<CNSTR != -1, Eigen::Vector<HomogType, CNSTR/2 * NSTOKES>, Eigen::VectorXd>::type;

        using VectorViewH = typename Eigen::Map<VectorH>;


        using MatrixHLHS = typename std::conditional<CNSTR != -1, Eigen::Matrix<double, CNSTR/2 * NSTOKES + 1, CNSTR/2 * NSTOKES + 1>, Eigen::MatrixXd>::type;
        using MatrixHRHS = typename std::conditional<CNSTR != -1, Eigen::Matrix<double, CNSTR/2 * NSTOKES + 1, -1>, Eigen::MatrixXd>::type;

        // These four aren't used anymore since we should just kill the non-greens function solution
        Eigen::MatrixXd particular_rhs;
        Matrix particular_A;
        Matrix particular_b;
        std::vector<Matrix> particular_d_A;

        Matrix h_eigmtx_destroy;
        Matrix h_MX_plus;
        Matrix h_MX_minus;
        Vector h_eigvalsq;
        Vector h_reigval_imag;
        Matrix h_identity;
        MatrixHLHS h_lhs;
        std::vector<MatrixHRHS> h_rhs;

        Eigen::PartialPivLU<MatrixHLHS> h_partiallu;
        Eigen::FullPivLU<MatrixHLHS> h_fullpivlu;


        std::vector<MatrixHRHS> h_d_X_d_k;

        std::vector<VectorLayerDual<double>> p_Qplus, p_Qminus;
        Dual<double> p_Cplus, p_Cminus;

        sktran_do_detail::TripleProductDerivativeHolder<NSTOKES> h_l_upwelling;
        sktran_do_detail::TripleProductDerivativeHolder<NSTOKES> h_l_downwelling;

        sktran_do_detail::InhomogeneousSourceHolder<NSTOKES> p_d_temp;
        sktran_do_detail::InhomogeneousSourceHolder<NSTOKES> p_d_temp2;

        VectorDim1<BVPMatrixDenseBlock<NSTOKES>> d_mat;
        VectorDim1<Eigen::VectorXd> d_b;
    };

    // Data that each thread will need, reused across wavelengths.
    template <int NSTOKES, int CNSTR=-1>
    class ThreadData {
    private:
        mutable VectorDim2<LayerSolution<NSTOKES, CNSTR>> m_rte_solution;
        mutable VectorDim2<LegendreSumMatrixStorage<NSTOKES>> m_legendre_sum_storage;
        mutable VectorDim1<PostProcessingCache<NSTOKES>> m_postprocessing_cache;
        mutable VectorDim1<LayerCache<NSTOKES>> m_layer_cache;
        mutable InputDerivatives<NSTOKES> m_input_derivatives;
        mutable RTEMemoryCache<NSTOKES, CNSTR> m_rte_cache;

    public:
        ThreadData(uint NLYR, uint NSTR) {
            m_rte_solution.resize(NLYR);
            for(auto& soln : m_rte_solution) {
                soln.resize(NSTR);
            }
            m_legendre_sum_storage.resize(NLYR);
            for(auto& leg : m_legendre_sum_storage) {
                leg.resize(NSTR);
            }

            m_postprocessing_cache.resize(NLYR);
            m_layer_cache.resize(NLYR, NSTR);
        }

        std::vector<LayerSolution<NSTOKES, CNSTR>>& rte_solution(uint layerindex) const {
            return m_rte_solution[layerindex];
        }

        std::vector<LegendreSumMatrixStorage<NSTOKES>>& legendre_sum_storage(uint layerindex) const {
            return m_legendre_sum_storage[layerindex];
        }

        PostProcessingCache<NSTOKES>& postprocessing_cache(uint layerindex) const {
            return m_postprocessing_cache[layerindex];
        }

        LayerCache<NSTOKES>& layer_cache(uint layerindex) const {
            return m_layer_cache[layerindex];
        }

        InputDerivatives<NSTOKES>& input_derivatives() const {
            return m_input_derivatives;
        }

        RTEMemoryCache<NSTOKES, CNSTR>& rte_cache() const {
            return m_rte_cache;
        }

    };

    // Pool of memory for a single Engine instance.  One ThreadData object is instantiated for each thread that
    // is then intended to be reused across wavelengths.  Currently the engine does not know exactly how many threads
    // will be used at the time of calling init, and so we instantiate the maximum number of available threads, but
    // this should be fixed
    template <int NSTOKES, int CNSTR=-1>
    class MemoryPool {
    private:
        // Use a map because we emplace a new object and we don't want reallocs triggering a move
        mutable std::map<int, ThreadData<NSTOKES, CNSTR>> m_threaddata;
        uint m_nlyr;
        uint m_nstr;

    public:
        MemoryPool() {}

        MemoryPool(uint NLYR, uint NSTR) {
            init(NLYR, NSTR);
        }

        void init(uint NLYR, uint NSTR) {
            int num_threads = omp_get_max_threads();

            m_nlyr = NLYR;
            m_nstr = NSTR;

            for(int k = 0; k < num_threads; ++k) {
                m_threaddata.emplace(k, ThreadData<NSTOKES, CNSTR>(NLYR, NSTR));
            }
        }

        ThreadData<NSTOKES, CNSTR>& thread_data() const {
            int thread_num = omp_get_thread_num();

            auto it = m_threaddata.find(thread_num);
            if( it == m_threaddata.end()) {
                m_threaddata.emplace(thread_num, ThreadData<NSTOKES, CNSTR>(m_nlyr, m_nstr));
            }

            return m_threaddata.at(thread_num);
        }
    };

}