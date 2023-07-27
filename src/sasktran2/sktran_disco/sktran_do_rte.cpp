#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_rte.h"
#include <thread>


template <int NSTOKES, int CNSTR>
sasktran_disco::RTESolver<NSTOKES, CNSTR>::RTESolver(const PersistentConfiguration<NSTOKES, CNSTR>& config, OpticalLayerArray<NSTOKES, CNSTR>& layers)
	: RTESProperties<NSTOKES>(config),
	m_layers(layers),
    m_cache(config.pool().thread_data().rte_cache())
{
	// Configure azimuth expansion
	registerAzimuthDependency(*this->M_LP_CSZ);
	registerAzimuthDependency(layers);
	// Initialize tracker for which orders have been solved
	m_is_solved.resize(this->M_NSTR, false);

    m_use_greens_function = true;

	configureCache();
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::RTESolver<NSTOKES, CNSTR>::configureCache() {
	const uint N = this->M_NSTR / 2;
	m_cache.h_eigmtx_destroy.resize(N * NSTOKES, N * NSTOKES);
	m_cache.h_MX_minus.resize(N * NSTOKES, N * NSTOKES);
	m_cache.h_MX_plus.resize(N * NSTOKES, N * NSTOKES);
	m_cache.h_eigvalsq.resize(N * NSTOKES);
	m_cache.h_reigval_imag.resize(N * NSTOKES);
	m_cache.h_identity.resize(N * NSTOKES, N * NSTOKES);
	m_cache.h_identity.setIdentity();
	m_cache.h_lhs.resize(N * NSTOKES + 1, N * NSTOKES + 1);

	m_cache.h_rhs.resize(m_layers.numLayers());
	m_cache.h_d_X_d_k.resize(m_layers.numLayers());

    m_cache.p_Qplus.resize(m_layers.numLayers());
    m_cache.p_Qminus.resize(m_layers.numLayers());

	for (uint k = 0; k < m_layers.numLayers(); ++k) {
		const uint numderiv = (uint)m_layers.inputDerivatives().numDerivativeLayer(k);
        const uint layerstart = (uint)m_layers.inputDerivatives().layerStartIndex(k);
		m_cache.h_rhs[k].resize(N * NSTOKES + 1, numderiv);
		m_cache.h_d_X_d_k[k].resize(N * NSTOKES + 1, numderiv);

        m_cache.p_Qplus[k].resize(this->M_NSTR/2 * NSTOKES, numderiv, k, layerstart);
        m_cache.p_Qminus[k].resize(this->M_NSTR/2 * NSTOKES, numderiv, k, layerstart);

        m_cache.p_norm.emplace_back(LayerDual<double>(numderiv, k, layerstart));
    }

    m_cache.h_l_downwelling.resize(this->M_NSTR);
    m_cache.h_l_upwelling.resize(this->M_NSTR);

    m_cache.p_d_temp.resize(this->M_NSTR);
    m_cache.p_d_temp2.resize(this->M_NSTR);

    uint numDeriv = (uint)m_layers.inputDerivatives().numDerivative();

    m_cache.p_Cplus.resize(numDeriv, false);
    m_cache.p_Cminus.resize(numDeriv, false);

    if(m_cache.d_mat.size() == 0) {
        m_cache.d_mat.reserve(numDeriv);
        m_cache.d_b.reserve(numDeriv);
        for (uint i = 0; i < numDeriv; ++i) {
            m_cache.d_mat.emplace_back(BVPMatrixDenseBlock<NSTOKES>(m_layers.inputDerivatives().layerDerivatives()[i].layer_index, this->M_NSTR, this->M_NLYR));
            m_cache.d_b.emplace_back(Eigen::VectorXd(this->M_NSTR * NSTOKES * this->M_NLYR));
        }
    }

    m_cache.bvp_b.resize(this->M_NSTR * NSTOKES * this->M_NLYR);
    m_cache.bvp_mat = std::make_unique<la::BVPMatrix<NSTOKES>>(this->M_NSTR, this->M_NLYR);

    m_cache.bvp_d_rhs.resize(m_cache.bvp_mat->N(), numDeriv);
    m_cache.bvp_temp.resize(m_cache.bvp_b.size());
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::RTESolver<NSTOKES, CNSTR>::solve(AEOrder m)
{
	// Environment setup
	using namespace sasktran_disco;
	// Check if its already solved. If so then return true
	if(m_is_solved[m]) return;

	// Otherwise do the calculation
	configureAEOrder(m);
	for(int p = 0; p < static_cast<int>(this->M_NLYR); ++p) {
		auto& layer = m_layers[p];
		layer.solution(m).configure(this->M_NSTR, p, m_layers.inputDerivatives());
		// Calculate homogeneous and particular solutions
		solveHomogeneous(m, layer);
        solveParticularGreen(m, layer);
	}
	// Calculate coefficients that satisfy boundary conditions
	solveBVP(m);
	m_is_solved[m] = true;

	postProcessAEOrder(m);
}

template<int NSTOKES, int CNSTR>
void sasktran_disco::RTESolver<NSTOKES, CNSTR>::assignHomogenousSplusMinus(AEOrder m, OpticalLayer<NSTOKES, CNSTR>& layer
) {
    // Constructs the Splus and Sminus matrices in the homogeneous solution

    const uint N = this->M_NSTR / 2;
    uint numDeriv = (uint)m_layers.inputDerivatives().numDerivativeLayer(layer.index());
    const auto& derivIter = m_layers.inputDerivatives().layerDerivatives().begin() + m_layers.inputDerivatives().layerStartIndex(layer.index());

    double zeta, eta, d_zeta, d_eta;

    typename std::conditional<NSTOKES==1, double, Eigen::Matrix<double, NSTOKES, NSTOKES>>::type d_l_upwelling, d_l_downwelling;

    sasktran_disco::TripleProductDerivativeHolder<NSTOKES>& l_upwelling = m_cache.h_l_upwelling;
    sasktran_disco::TripleProductDerivativeHolder<NSTOKES>& l_downwelling = m_cache.h_l_downwelling;

    // i,j are indexes in the block matrix that is composed of 4x4 stokes blocks
    for(uint i = 0; i < N; ++i) {
        for(uint j = 0; j < N; ++j) {
            // Cache upwelling and downwelling legendre sums
            layer.inplace_scatPhaseFAndDerivative(m, i, j, l_upwelling);
            layer.inplace_scatPhaseFAndDerivative(m, i, j + N, l_downwelling);

            // TODO: We can probably refactor this to emplace the values directly in S_plus/S_minus and derivative
            // the values directly in as well
            for(uint s1 = 0; s1 < NSTOKES; ++s1) {
                for(uint s2 = 0; s2 < NSTOKES; ++s2) {
                    // Calculate indicies in the full matrix
                    uint ii, jj;
                    if (j > i) {
                        // Have to transpose the matrix
                        ii = NSTOKES * i + s2;
                        jj = NSTOKES * j + s1;
                    }
                    else {
                        ii = NSTOKES * i + s1;
                        jj = NSTOKES * j + s2;
                    }

                    if constexpr(NSTOKES == 1) {
                        zeta = ((*this->M_WT)[j] * l_upwelling.value - kronDelta(ii, jj)) / (*this->M_MU)[i];
                        eta = (*this->M_WT)[j] * l_downwelling.value / (*this->M_MU)[i];
                    } else {
                        zeta = ((*this->M_WT)[j] * l_upwelling.value(s1, s2) - kronDelta(ii, jj)) / (*this->M_MU)[i];
                        eta = (*this->M_WT)[j] * l_downwelling.value(s1, s2) / (*this->M_MU)[i];
                    }

                    layer.solution(m).cache.s_plus()(ii, jj) = -1.0 * (zeta + eta);
                    layer.solution(m).cache.s_minus()(ii, jj) = -1.0 * (zeta - eta);
                }
            }

            // Because the derivative propagation is so simple we just do it manually here
            for (uint k = 0; k < numDeriv; ++k) {
                l_upwelling.reduce(*(derivIter + k), d_l_upwelling);
                l_downwelling.reduce(*(derivIter + k), d_l_downwelling);

                for(uint s1 = 0; s1 < NSTOKES; ++s1) {
                    for (uint s2 = 0; s2 < NSTOKES; ++s2) {
                        int ii,jj;
                        if (j > i) {
                            // Have to transpose the matrix
                            ii = NSTOKES * i + s2;
                            jj = NSTOKES * j + s1;
                        }
                        else {
                            ii = NSTOKES * i + s1;
                            jj = NSTOKES * j + s2;
                        }

                        if constexpr(NSTOKES == 1) {
                            d_zeta = ((*this->M_WT)[j] * d_l_upwelling) / (*this->M_MU)[i];
                            d_eta = ((*this->M_WT)[j] * d_l_downwelling) / (*this->M_MU)[i];
                        } else {
                            d_zeta = ((*this->M_WT)[j] * d_l_upwelling(s1, s2)) / (*this->M_MU)[i];
                            d_eta = (*this->M_WT)[j] * d_l_downwelling(s1, s2) / (*this->M_MU)[i];
                        }

                        layer.solution(m).d_cache[k].s_plus()(ii, jj) = -1.0 * (d_zeta + d_eta);
                        layer.solution(m).d_cache[k].s_minus()(ii, jj) = -1.0 * (d_zeta - d_eta);
                    }
                }
            }
        }
    }
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::RTESolver<NSTOKES, CNSTR>::linearizeHomogeneous(AEOrder m,
                                                                     OpticalLayer<NSTOKES, CNSTR>& layer)
{
    // Linearizes the homogeneous solution

    auto& solution = layer.solution(m);	
    const uint N = this->M_NSTR / 2;
    uint numDeriv = (uint)m_layers.inputDerivatives().numDerivativeLayer(layer.index());

    if(numDeriv == 0) {
        return;
    }

    Matrix& S_plus = solution.cache.s_plus();
    Matrix& eigmtx = solution.cache.eigmtx();

    Vector& eigvalsq = m_cache.h_eigvalsq;
    VectorViewH eigval(solution.value.eigval().data(), N*NSTOKES, 1);
    Matrix& MX_plus = m_cache.h_MX_plus;

    // Now calculate the linearization of homog_plus and homog_minus
    // For each solution we have a N+1 system of equations for d_X and d_k
    Matrix& identity = m_cache.h_identity;
    for (SolutionIndex j = 0; j < N * NSTOKES; ++j)
    {
        // The LHS of the system is the same for every linearization
        MatrixHLHS& lhs = m_cache.h_lhs;
        // For each row add in an equation
        for (uint row = 0; row < N * NSTOKES; ++row) {

            lhs(row, Eigen::seq(0, Eigen::last - 1)).noalias() = eigmtx(row, Eigen::all);
            lhs(row, row) -= eigvalsq(j);

            lhs(row, N*NSTOKES) = -2 * eigval(j) * MX_plus(row, j);
        }
        // The final equation is that d_X is orthogonal to X,
        lhs(N*NSTOKES, Eigen::seq(0, Eigen::last - 1)).noalias() = MX_plus(Eigen::all, j);
        // Have to make sure to set the last element to 0
        lhs(N*NSTOKES, N*NSTOKES) = 0;

        // for each linearization we have a different RHS
        MatrixHRHS& rhs = m_cache.h_rhs[layer.index()];
        for (uint i = 0; i < numDeriv; ++i) {
            rhs(Eigen::seq(0, Eigen::last - 1), i).noalias() = -1.0 * (layer.solution(m).d_cache[i].eigmtx()) * MX_plus(Eigen::all, j); // Regular system
            rhs(N*NSTOKES, i) = 0.0; // Orthogonality
        }
        MatrixHRHS& d_X_d_k = m_cache.h_d_X_d_k[layer.index()];
        if (NSTOKES == 1) {
            m_cache.h_partiallu.compute(lhs);
            d_X_d_k.noalias() = m_cache.h_partiallu.solve(rhs);
        }
        else {
            // We can get a degeneracy where Q or U is identical to 0 that requires the full pivot
            m_cache.h_fullpivlu.compute(lhs);
            d_X_d_k.noalias() = m_cache.h_fullpivlu.solve(rhs);
        }

        // Store the derivative result for each linearization
        for (uint i = 0; i < numDeriv; ++i) {
            // Store the derivative for the eigen value
            double d_k = d_X_d_k(Eigen::last, i);
            solution.value.dual_eigval().deriv(i, j) = d_k;

            // Begin constructing the derivatives for W_plus and W_minus
            // These get stored in one row of the resulting derivative
            auto d_W_plus = solution.value.dual_homog_plus().deriv.row(i);
            auto d_W_minus = solution.value.dual_homog_minus().deriv.row(i);

            d_W_plus(Eigen::seq(j * N * NSTOKES, (j + 1) * N * NSTOKES - 1)).noalias() = 0.5 * ((-d_k / eigvalsq(j)) * S_plus + 1 / eigval(j) * (layer.solution(m).d_cache[i].s_plus())).lazyProduct(MX_plus(Eigen::all, j));
            d_W_plus(Eigen::seq(j * N * NSTOKES, (j + 1) * N * NSTOKES - 1)).noalias() += 0.5 * (identity + 1 / eigval(j) * S_plus).lazyProduct(d_X_d_k(Eigen::seq(0, Eigen::last - 1), i));

            d_W_minus(Eigen::seq(j * N * NSTOKES, (j + 1) * N * NSTOKES - 1)).noalias() = 0.5 * ((d_k / eigvalsq(j)) * S_plus - 1 / eigval(j) * (layer.solution(m).d_cache[i].s_plus())).lazyProduct(MX_plus(Eigen::all, j));
            d_W_minus(Eigen::seq(j * N * NSTOKES, (j + 1) * N * NSTOKES - 1)).noalias() += 0.5 * (identity - 1 / eigval(j) * S_plus).lazyProduct(d_X_d_k(Eigen::seq(0, Eigen::last - 1), i));
        }
    }
}

#ifdef SASKTRAN_DISCO_FULL_COMPILE
template <>
void sasktran_disco::RTESolver<1, 2>::linearizeHomogeneous(AEOrder m, OpticalLayer<1, 2>& layer)
{
    // Linearizes the homogeneous solution for the 2 stream special case

    auto& solution = layer.solution(m);
    uint numDeriv = (uint)m_layers.inputDerivatives().numDerivativeLayer(layer.index());

    if(numDeriv == 0) {
        return;
    }

    Matrix& S_plus = solution.cache.s_plus();
    Matrix& eigmtx = solution.cache.eigmtx();

    Vector& eigvalsq = m_cache.h_eigvalsq;
    VectorViewH eigval(solution.value.eigval().data(), 1, 1);
    Matrix& MX_plus = m_cache.h_MX_plus;

    // Now calculate the linearization of homog_plus and homog_minus
    // For each solution we have a N+1 system of equations for d_X and d_k
    Matrix& identity = m_cache.h_identity;
    // The LHS of the system is the same for every linearization
    MatrixHLHS& lhs = m_cache.h_lhs;
    lhs(0, 0) = eigmtx(0, 0) - eigvalsq(0);
    lhs(0, 1) = -2 * eigval(0) * MX_plus(0, 0);

    // The final equation is that d_X is orthogonal to X,
    lhs(1, 0) = MX_plus(0, 0);
    // Have to make sure to set the last element to 0
    lhs(1, 1) = 0;

    // for each linearization we have a different RHS
    MatrixHRHS& rhs = m_cache.h_rhs[layer.index()];
    for (uint i = 0; i < numDeriv; ++i) {
        rhs(0, i) = -1.0 * (layer.solution(m).d_cache[i].eigmtx()(0, 0)) * MX_plus(0, 0); // Regular system
        rhs(1, i) = 0.0; // Orthogonality
    }
    MatrixHRHS& d_X_d_k = m_cache.h_d_X_d_k[layer.index()];
    //m_cache.h_partiallu.compute(lhs);
    //d_X_d_k.noalias() = m_cache.h_partiallu.solve(rhs);

    double factor = (lhs(0, 1) - lhs(0, 0) * lhs(1, 1)) / lhs(1, 0);
    double factor2 = -lhs(1, 1) / lhs(1, 0);
    for(uint i = 0; i < numDeriv; ++i) {
        d_X_d_k(0, i) = 0;
        d_X_d_k(1, i) = rhs(0, i) / lhs(0, 1);
    }


    // Store the derivative result for each linearization
    for (uint i = 0; i < numDeriv; ++i) {
        // Store the derivative for the eigen value
        double d_k = d_X_d_k(Eigen::last, i);
        solution.value.dual_eigval().deriv(i, 0) = d_k;

        // Begin constructing the derivatives for W_plus and W_minus
        // These get stored in one row of the resulting derivative
        auto d_W_plus = solution.value.dual_homog_plus().deriv.row(i);
        auto d_W_minus = solution.value.dual_homog_minus().deriv.row(i);

        d_W_plus(0) = 0.5 * ((-d_k / eigvalsq(0)) * S_plus(0, 0) + 1 / eigval(0) * (layer.solution(m).d_cache[i].s_plus()(0, 0))) * (MX_plus(0, 0));
        d_W_plus(0) += 0.5 * (1 + 1 / eigval(0) * S_plus(0, 0)) * (d_X_d_k(0, i));

        d_W_minus(0) = 0.5 * ((d_k / eigvalsq(0)) * S_plus(0 , 0) - 1 / eigval(0) * (layer.solution(m).d_cache[i].s_plus())(0, 0)) * (MX_plus(0, 0));
        d_W_minus(0) += 0.5 * (1 - 1 / eigval(0) * S_plus(0, 0)) * (d_X_d_k(0, i));
    }
}
#endif

template <int NSTOKES, int CNSTR>
void sasktran_disco::RTESolver<NSTOKES, CNSTR>::solveHomogeneous(AEOrder m, OpticalLayer<NSTOKES, CNSTR>& layer)
{
	// Setup up calculation environment
	const uint N = this->M_NSTR / 2;
	auto& solution = layer.solution(m);
	
	uint numDeriv = (uint)m_layers.inputDerivatives().numDerivativeLayer(layer.index());
	const auto& derivIter = m_layers.inputDerivatives().layerDerivatives().begin() + m_layers.inputDerivatives().layerStartIndex(layer.index());

	// Sum and difference matricies, S_plus and S_minus, See eq 14
    Matrix& S_plus = solution.cache.s_plus();
    Matrix& S_minus = solution.cache.s_minus();

    assignHomogenousSplusMinus(m, layer);

	// Eigenvalue problem (16)
    Matrix& eigmtx = solution.cache.eigmtx();
	Matrix& eigmtx_destroyable = m_cache.h_eigmtx_destroy; // Destructible copy of eigmtx

	// Do eigenmatrix = S_minus*S_plus
    eigmtx.noalias() = S_minus * S_plus;

	// Calculate the derivative of the eigenmatrix
	// d_eigenmatrix = d_S_minus*S_plus + S_minus*d_S_plus
	for (uint i = 0; i < numDeriv; ++i) {
		// First multiply d_S_minus*S_plus and store it in eigmtx
        solution.d_cache[i].eigmtx().noalias() = solution.d_cache[i].s_minus() * S_plus;

		// Then do the second multiplication and add it to the previous result
        solution.d_cache[i].eigmtx() += S_minus * solution.d_cache[i].s_plus();
	}
	
	// Copy matrix
	eigmtx_destroyable = eigmtx;

	// Sum and difference eigenvectors, these are MX_plus and MX_minus in eq (16)
	Matrix& MX_plus = m_cache.h_MX_plus;
	Matrix& MX_minus = m_cache.h_MX_minus;

	// Squared eigenvalues
	Vector& eigvalsq = m_cache.h_eigvalsq; // eigvalsq(j) = eigval(j) * eigval*(j)
	Vector& reigval_imag = m_cache.h_reigval_imag; // Imaginary components of eigenvalues, zero in the case of NSTOKES=1 or pure Rayleigh scattering

    if(SASKTRAN_DISCO_USE_EIGEN_EIGENSOLVER || CNSTR != -1) {
        Eigen::EigenSolver<Matrix> es(eigmtx_destroyable);

        auto eigeninfo = es.info();
        if(eigeninfo != Eigen::Success) {
            throw InternalRuntimeError("Error computing the homogeneous solution");
        }

        eigvalsq = es.eigenvalues().real();
        MX_plus = es.eigenvectors().real();
        reigval_imag.setZero();
    }
    else {
        // Use lapack computation directly
        // Compute eigensolution
        lapack_int errorcode;

        #ifdef SKTRAN_USE_ACCELERATE
            Eigen::VectorXd work(N*NSTOKES*4);
            int worksize = N*NSTOKES*4;
            char jobvl = 'N';
            char jobvr = 'V';
            int n = N * NSTOKES;
            int one = 1;

            dgeev_(&jobvl, &jobvr, &n, eigmtx_destroyable.data(),
                   &n, eigvalsq.data(), reigval_imag.data(),
                   nullptr,
                   &one, MX_plus.data(), &n,
                   work.data(),
                   &worksize, &errorcode);

        #else
            errorcode = LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'V', N * NSTOKES, eigmtx_destroyable.data(), N * NSTOKES, eigvalsq.data(), reigval_imag.data(), nullptr, 1, MX_plus.data(), N * NSTOKES);
        #endif
        // Check for errors
        if(errorcode != 0) {
            if(errorcode < 0) {
                throw InternalRuntimeError("An argument to LAPACKE_dgeev had an illegal argument in sasktran_disco::RTESolver::SolveHomogeneous");
            } else {
                throw InternalRuntimeError("LAPACKE_dgeev failed to compute all solutions");
            }
        }
    }

	// MX_minus = S_plus * MX_plus * inv{diag{eigenvalues}} from eq (14). Do S_plus * MX_plus here
    MX_minus = S_plus * MX_plus;

	VectorViewH eigval(solution.value.eigval().data(), N*NSTOKES, 1);
	MatrixViewH homog_plus(solution.value.homog_plus().data(), N * NSTOKES, N * NSTOKES);
	MatrixViewH homog_minus(solution.value.homog_minus().data(), N * NSTOKES, N * NSTOKES);

	// Calculate the homogeneous solutions W_plus and W_minus and store them in the layer solution
	for (SolutionIndex j = 0; j < N * NSTOKES; ++j)
	{
		if (eigvalsq(j) <= 0 && NSTOKES == 1)
			throw InternalRuntimeError("An homogeneous solution was found to be imaginary. An insufficient number of streams is likely.");
		eigval(j) = sqrt(abs(eigvalsq(j)));
		for (uint i = 0; i < N * NSTOKES; ++i) {
            homog_plus(i, j) = 0.5 * (MX_plus(i, j) + MX_minus(i, j) / eigval(j));
            homog_minus(i, j) = 0.5 * (MX_plus(i, j) - MX_minus(i, j) / eigval(j));
		}
	}

    linearizeHomogeneous(m, layer);
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::RTESolver<NSTOKES, CNSTR>::assignParticularQ(AEOrder m, const OpticalLayer<NSTOKES, CNSTR> &layer, VectorLayerDual<double>&Qplus,
                                                                  VectorLayerDual<double>&Qminus) {
    InhomogeneousSourceHolder<NSTOKES>& d_temp = m_cache.p_d_temp;
    InhomogeneousSourceHolder<NSTOKES>& d_temp2 = m_cache.p_d_temp2;

    typename std::conditional<NSTOKES==1, double, Eigen::Vector<double, NSTOKES>>::type dholder1, dholder2;

    const auto &derivIter = m_layers.inputDerivatives().layerDerivatives().begin() +
                            m_layers.inputDerivatives().layerStartIndex(layer.index());
    uint numLayerDeriv = (uint) m_layers.inputDerivatives().numDerivativeLayer(layer.index());

    for (uint i = 0; i < this->M_NSTR / 2; ++i) {
        layer.singleScatST(m, (*this->M_LP_MU)[m][i], d_temp, d_temp2);

        for(int s = 0; s < NSTOKES; ++s) {
            int ii = i*NSTOKES + s;

            if constexpr(NSTOKES==1) {
                Qplus.value(ii) = d_temp2.value * (*this->M_WT)[i];
                Qminus.value(ii) = d_temp.value * (*this->M_WT)[i];
            } else {
                Qplus.value(ii) = d_temp2.value[s] * (*this->M_WT)[i];
                Qminus.value(ii) = d_temp.value[s] * (*this->M_WT)[i];
            }

        }

        for(uint k = 0; k < numLayerDeriv; ++k) {
            d_temp2.reduce(*(derivIter + k), dholder2);
            d_temp.reduce(*(derivIter + k), dholder1);

            for(int s = 0; s < NSTOKES; ++s) {
                int ii = i*NSTOKES + s;

                if constexpr(NSTOKES==1) {
                    Qplus.deriv(k, ii) = dholder2 * (*this->M_WT)[i];
                    Qminus.deriv(k, ii) = dholder1 * (*this->M_WT)[i];
                } else {
                    Qplus.deriv(k, ii) = dholder2[s] * (*this->M_WT)[i];
                    Qminus.deriv(k, ii) = dholder1[s] * (*this->M_WT)[i];
                }
            }
        }
    }
}

#ifdef SASKTRAN_DISCO_FULL_COMPILE
template <>
void sasktran_disco::RTESolver<1, 2>::solveParticularGreen(AEOrder m, OpticalLayer<1, 2>& layer) {
    // Setup up calculation enviroment
    LayerIndex p = layer.index();
    auto &solution = layer.solution(m);

    solution.value.set_use_green_function(true);

    uint numLayerDeriv = (uint) m_layers.inputDerivatives().numDerivativeLayer(layer.index());
    uint numTotalDeriv = (uint) m_layers.inputDerivatives().numDerivative();
    uint layerStart = (uint) m_layers.inputDerivatives().layerStartIndex(p);

    const auto &average_secant = layer.dual_average_secant();
    const auto &thickness = layer.dual_thickness();
    const auto &transmission = layer.ceiling_beam_transmittanc();

    // Construct the Q vectors, we immediately scale by the
    // quadrature weights
    auto& Qplus = m_cache.p_Qplus[p];
    auto& Qminus = m_cache.p_Qminus[p];

    assignParticularQ(m, layer, Qplus, Qminus);

    const auto& homog_plus = solution.value.dual_homog_plus();
    const auto& homog_minus = solution.value.dual_homog_minus();

    auto& Aplus = solution.value.dual_green_A_plus();
    auto& Aminus = solution.value.dual_green_A_minus();

    auto& Gplus_bottom = solution.value.dual_Gplus_bottom();
    auto& Gplus_top = solution.value.dual_Gplus_top();

    auto& Gminus_bottom = solution.value.dual_Gminus_bottom();
    auto& Gminus_top = solution.value.dual_Gminus_top();

    // Normalization constant
    auto& norm = m_cache.p_norm[p];

    // Multipliers
    auto& Cplus = m_cache.p_Cplus;
    auto& Cminus = m_cache.p_Cminus;

    // For each homogeneous solution, add it's contribution in

    const auto& eigval = solution.value.dual_eigval();

    // Normalization constant calculation
    norm.value = (*this->M_WT)[0] * (*this->M_MU)[0] * (homog_plus.value(0) * homog_plus.value(0) - homog_minus.value(0) * homog_minus.value(0));

    Aplus.value(0) = Qplus.value(0) * homog_plus.value(0) + Qminus.value(0) * homog_minus.value(0);
    Aminus.value(0) = Qminus.value(0) * homog_plus.value(0)+ Qplus.value(0) * homog_minus.value(0);

    for(uint k = 0; k < numLayerDeriv; ++k) {
        norm.deriv(k) = (*this->M_WT)[0] * (*this->M_MU)[0] * (2.0 * homog_plus.deriv(k, 0) * homog_plus.value(0) - 2.0 * homog_minus.deriv(k, 0) * homog_minus.value(0));
        Aplus.deriv(k, 0) = Qplus.value(0) * homog_plus.deriv(k, 0) + Qplus.deriv(k, 0) * homog_plus.value(0);
        Aplus.deriv(k, 0) += Qminus.value(0) * homog_minus.deriv(k, 0) + Qminus.deriv(k, 0) * homog_minus.value(0);

        Aminus.deriv(k, 0) = Qminus.value(0) * homog_plus.deriv(k, 0) + Qminus.deriv(k, 0) * homog_plus.value(0);
        Aminus.deriv(k, 0) += Qplus.value(0) * homog_minus.deriv(k, 0) + Qplus.deriv(k, 0) * homog_minus.value(0);
    }

    Aplus.value(0) /= norm.value;
    Aminus.value(0) /= norm.value;

    for(uint k = 0; k < numLayerDeriv; ++k) {
        Aplus.deriv(k, 0) = Aplus.deriv(k, 0) / norm.value - norm.deriv(k) * Aplus.value(0) / norm.value;
        Aminus.deriv(k, 0) = Aminus.deriv(k, 0) / norm.value - norm.deriv(k) * Aminus.value(0) / norm.value;
    }

    Cplus.value = 0.0;
    Cminus.value = 0.0;

    double exp_thickness_eigval = exp(-thickness.value * eigval.value(0));
    double exp_thickness_secant = exp(-thickness.value*average_secant.value);

    // If average secant is close to eigval then we evaluate Cplus or Cminus with a taylor series expansion instead
    if(abs(average_secant.value - eigval.value(0)) > SKTRAN_DO_GREENS_EPS) {
        Cplus.value = transmission.value * (exp_thickness_eigval - exp_thickness_secant) /
                      (average_secant.value - eigval.value(0));

        Cplus.deriv.noalias() = (exp_thickness_eigval - exp_thickness_secant) /
                                 (average_secant.value - eigval.value(0)) * transmission.deriv;
        Cplus.deriv.noalias() += average_secant.deriv * (transmission.value * thickness.value * exp(-1.0 * thickness.value * average_secant.value) - Cplus.value) / (average_secant.value - eigval.value(0));

        for(uint k = 0; k < numLayerDeriv; ++k) {
            Cplus.deriv(k + layerStart) += eigval.deriv(k, 0) / (average_secant.value - eigval.value(0)) *
                                           (Cplus.value - transmission.value * thickness.value *
                                                                  exp_thickness_eigval);
            Cplus.deriv(k + layerStart) -=
                    thickness.deriv(k) * transmission.value / (average_secant.value - eigval.value(0)) *
                    (eigval.value(0) * exp_thickness_eigval -
                     average_secant.value * exp_thickness_secant);
        }
    } else {
        // Second order taylor expansion of Cplus
        Cplus.value = transmission.value * exp(-1.0 * thickness.value * eigval.value(0)) * thickness.value * (1 - thickness.value / 2 * (average_secant.value - eigval.value(0)));

        Cplus.deriv.noalias() = exp(-1.0 * thickness.value * eigval.value(0)) * thickness.value * (1 - thickness.value / 2 * (average_secant.value - eigval.value(0))) * transmission.deriv;;
        Cplus.deriv.noalias() += -1.0 * average_secant.deriv * thickness.value / 2.0 * thickness.value * transmission.value * exp(-1.0 * thickness.value * eigval.value(0));

        for(uint k = 0; k < numLayerDeriv; ++k) {
            Cplus.deriv(k + layerStart) += eigval.deriv(k, 0) * Cplus.value * -1.0 * thickness.value;
            Cplus.deriv(k + layerStart) += eigval.deriv(k, 0) * thickness.value / 2.0 * thickness.value * transmission.value * exp(-1.0 * thickness.value * eigval.value(0));

            Cplus.deriv(k + layerStart) += thickness.deriv(k) * transmission.value * exp(-1.0 * thickness.value * eigval.value(0)) * (1 - thickness.value * (average_secant.value - eigval.value(0)));
            Cplus.deriv(k + layerStart) += -1.0 * eigval.value(0) * Cplus.value * thickness.deriv(k);
        }
    }

    if(abs(average_secant.value + eigval.value(0)) > SKTRAN_DO_GREENS_EPS) {

        Cminus.value = transmission.value * (1 - exp_thickness_secant *
                                                         exp_thickness_eigval) /
                       (average_secant.value + eigval.value(0));

        Cminus.deriv.noalias() =  (1 - exp_thickness_secant *
                                               exp_thickness_eigval) /
                                   (average_secant.value + eigval.value(0)) * transmission.deriv;
        Cminus.deriv.noalias() += average_secant.deriv / (average_secant.value + eigval.value(0)) *
                                  (transmission.value * thickness.value *
                                   exp_thickness_eigval*exp_thickness_secant - Cminus.value);

        for (uint k = 0; k < numLayerDeriv; ++k) {
            Cminus.deriv(k + layerStart) += eigval.deriv(k, 0) / (average_secant.value + eigval.value(0)) *
                                            (transmission.value * thickness.value *
                                                     exp_thickness_eigval*exp_thickness_secant -
                                             Cminus.value);
            Cminus.deriv(k + layerStart) += thickness.deriv(k) * transmission.value *
                    exp_thickness_eigval*exp_thickness_secant;
        }
    } else {
        // Second order taylor expansion of CMinus
        Cminus.value = transmission.value * thickness.value * (1 - thickness.value / 2 * (average_secant.value + eigval.value(0)));

        Cminus.deriv = transmission.deriv * thickness.value * (1 - thickness.value / 2 * (average_secant.value + eigval.value(0)));
        Cminus.deriv += average_secant.deriv * -1.0 * thickness.value / 2 * thickness.value * transmission.value;

        for(uint k = 0; k < numLayerDeriv; ++k) {
            Cminus.deriv(k + layerStart) += eigval.deriv(k, 0) * transmission.value * thickness.value * thickness.value / 2 * -1.0;
            Cminus.deriv(k + layerStart) += thickness.deriv(k) * transmission.value * (1 - thickness.value * (average_secant.value + eigval.value(0)));
        }
    }

    Gplus_top.value(0) = Aminus.value(0) * Cminus.value * homog_minus.value(0);
    Gminus_top.value(0) = Aminus.value(0) * Cminus.value * homog_plus.value(0);

    Gplus_bottom.value(0) = Aplus.value(0) * Cplus.value * homog_plus.value(0);
    Gminus_bottom.value(0) = Aplus.value(0) * Cplus.value * homog_minus.value(0);

    Gplus_top.deriv(Eigen::all, 0).noalias() = Aminus.value(0) * homog_minus.value(0) * Cminus.deriv;
    Gminus_top.deriv(Eigen::all, 0).noalias() = Aminus.value(0) * homog_plus.value(0) * Cminus.deriv;
    Gplus_bottom.deriv(Eigen::all, 0).noalias() = Aplus.value(0) * homog_plus.value(0) * Cplus.deriv;
    Gminus_bottom.deriv(Eigen::all, 0).noalias() = Aplus.value(0) * homog_minus.value(0) * Cplus.deriv;

    for(uint k = 0; k < numLayerDeriv; ++k) {
        Gplus_top.deriv(layerStart + k, 0) += Aminus.deriv(k, 0) * Cminus.value * homog_minus.value(0);
        Gminus_top.deriv(layerStart + k, 0) += Aminus.deriv(k, 0) * Cminus.value * homog_plus.value(0);
        Gplus_bottom.deriv(layerStart + k, 0) += Aplus.deriv(k, 0) * Cplus.value * homog_plus.value(0);
        Gminus_bottom.deriv(layerStart + k, 0) += Aplus.deriv(k, 0) * Cplus.value * homog_minus.value(0);

        Gplus_top.deriv(layerStart + k, 0) += Aminus.value(0) * Cminus.value * homog_minus.deriv(k, 0);
        Gminus_top.deriv(layerStart + k, 0) += Aminus.value(0) * Cminus.value * homog_plus.deriv(k, 0);
        Gplus_bottom.deriv(layerStart + k, 0) += Aplus.value(0) * Cplus.value * homog_plus.deriv(k, 0);
        Gminus_bottom.deriv(layerStart + k, 0) += Aplus.value(0) * Cplus.value * homog_minus.deriv(k, 0);
    }
}
#endif

template <int NSTOKES, int CNSTR>
void sasktran_disco::RTESolver<NSTOKES, CNSTR>::solveParticularGreen(AEOrder m, OpticalLayer<NSTOKES, CNSTR>& layer) {
    // Setup up calculation enviroment
    LayerIndex p = layer.index();
    auto &solution = layer.solution(m);

    solution.value.set_use_green_function(true);

    const uint N = this->M_NSTR / 2;

    uint numLayerDeriv = (uint) m_layers.inputDerivatives().numDerivativeLayer(layer.index());
    uint numTotalDeriv = (uint) m_layers.inputDerivatives().numDerivative();
    uint layerStart = (uint) m_layers.inputDerivatives().layerStartIndex(p);

    const auto &average_secant = layer.dual_average_secant();
    const auto &thickness = layer.dual_thickness();
    const auto &transmission = layer.ceiling_beam_transmittanc();

    // Construct the Q vectors, we immediately scale by the
    // quadrature weights
    auto& Qplus = m_cache.p_Qplus[p];
    auto& Qminus = m_cache.p_Qminus[p];

    assignParticularQ(m, layer, Qplus, Qminus);

    const auto& homog_plus = solution.value.dual_homog_plus();
    const auto& homog_minus = solution.value.dual_homog_minus();

    auto& Aplus = solution.value.dual_green_A_plus();
    auto& Aminus = solution.value.dual_green_A_minus();

    Aplus.value.setZero();
    Aminus.value.setZero();

    Aplus.deriv.setZero();
    Aminus.deriv.setZero();

    auto& Gplus_bottom = solution.value.dual_Gplus_bottom();
    auto& Gplus_top = solution.value.dual_Gplus_top();

    auto& Gminus_bottom = solution.value.dual_Gminus_bottom();
    auto& Gminus_top = solution.value.dual_Gminus_top();

    Gplus_bottom.value.setZero();
    Gplus_top.value.setZero();
    Gminus_bottom.value.setZero();
    Gminus_top.value.setZero();

    Gplus_bottom.deriv.setZero();
    Gplus_top.deriv.setZero();
    Gminus_bottom.deriv.setZero();
    Gminus_top.deriv.setZero();

    // Normalization constant
    LayerDual<double> norm(numLayerDeriv, p, layerStart);

    // Multipliers
    auto& Cplus = m_cache.p_Cplus;
    auto& Cminus = m_cache.p_Cminus;

    // For each homogeneous solution, add it's contribution in
    for(SolutionIndex i = 0; i < N * NSTOKES; ++i) {
        uint h_start = i * N * NSTOKES;
        const auto& eigval = solution.value.dual_eigval();

        // Normalization constant calculation
        norm.value = 0.0;
        norm.deriv.setZero();

        for(uint j = 0; j < N * NSTOKES; ++j) {
            //double negation = sasktran_disco::stokes_negation_factor<NSTOKES>(j);
            double negation = 1.0;
            norm.value += (*this->M_WT)[j / NSTOKES] * (*this->M_MU)[j / NSTOKES] * (homog_plus.value(h_start + j) * homog_plus.value(h_start + j) - homog_minus.value(h_start + j) * homog_minus.value(h_start + j));

            Aplus.value(i) += Qplus.value(j) * homog_plus.value(h_start + j) + Qminus.value(j) * homog_minus.value(h_start + j) * negation;
            Aminus.value(i) += Qminus.value(j) * homog_plus.value(h_start + j) * negation + Qplus.value(j) * homog_minus.value(h_start + j);

            for(uint k = 0; k < numLayerDeriv; ++k) {
                norm.deriv(k) += (*this->M_WT)[j / NSTOKES] * (*this->M_MU)[j / NSTOKES] * (2.0 * homog_plus.deriv(k, h_start + j) * homog_plus.value(h_start + j) - 2.0 * homog_minus.deriv(k, h_start + j) * homog_minus.value(h_start + j));
                Aplus.deriv(k, i) += Qplus.value(j) * homog_plus.deriv(k, h_start + j) + Qplus.deriv(k, j) * homog_plus.value(h_start + j);
                Aplus.deriv(k, i) += Qminus.value(j) * homog_minus.deriv(k, h_start + j) * negation + Qminus.deriv(k, j) * homog_minus.value(h_start + j) * negation;

                Aminus.deriv(k, i) += Qminus.value(j) * homog_plus.deriv(k, h_start + j) * negation + Qminus.deriv(k, j) * homog_plus.value(h_start + j) * negation;
                Aminus.deriv(k, i) += Qplus.value(j) * homog_minus.deriv(k, h_start + j) + Qplus.deriv(k, j) * homog_minus.value(h_start + j);
            }

        }
        Aplus.value(i) /= norm.value;
        Aminus.value(i) /= norm.value;

        for(uint k = 0; k < numLayerDeriv; ++k) {
            Aplus.deriv(k, i) = Aplus.deriv(k, i) / norm.value - norm.deriv(k) * Aplus.value(i) / norm.value;
            Aminus.deriv(k, i) = Aminus.deriv(k, i) / norm.value - norm.deriv(k) * Aminus.value(i) / norm.value;
        }

        Cplus.value = 0.0;
        Cminus.value = 0.0;

        Cplus.deriv.setZero();
        Cminus.deriv.setZero();

        // If average secant is close to eigval then we evaluate Cplus or Cminus with a taylor series expansion instead
        if(abs(average_secant.value - eigval.value(i)) > SKTRAN_DO_GREENS_EPS) {
            Cplus.value = transmission.value * (exp(-thickness.value * eigval.value(i)) - exp(-thickness.value*average_secant.value)) /
                          (average_secant.value - eigval.value(i));

			Cplus.deriv.noalias() += (exp(-thickness.value * eigval.value(i)) - exp(-thickness.value * average_secant.value)) /
				(average_secant.value - eigval.value(i)) * transmission.deriv;
            Cplus.deriv.noalias() += average_secant.deriv * (transmission.value * thickness.value * exp(-1.0 * thickness.value * average_secant.value) - Cplus.value) / (average_secant.value - eigval.value(i));

            for(uint k = 0; k < numLayerDeriv; ++k) {
                Cplus.deriv(k + layerStart) += eigval.deriv(k, i) / (average_secant.value - eigval.value(i)) *
                                               (Cplus.value - transmission.value * thickness.value *
                                                              exp(-1.0 * thickness.value * eigval.value(i)));
                Cplus.deriv(k + layerStart) -=
                        thickness.deriv(k) * transmission.value / (average_secant.value - eigval.value(i)) *
                        (eigval.value(i) * exp(-1.0 * thickness.value * eigval.value(i)) -
                         average_secant.value * exp(-1.0 * thickness.value * average_secant.value));
            }
        } else {
            // Second order taylor expansion of Cplus
            Cplus.value = transmission.value * exp(-1.0 * thickness.value * eigval.value(i)) * thickness.value * (1 - thickness.value / 2 * (average_secant.value - eigval.value(i)));

			Cplus.deriv.noalias() += exp(-1.0 * thickness.value * eigval.value(i)) * thickness.value * (1 - thickness.value / 2 * (average_secant.value - eigval.value(i))) * transmission.deriv;;
            Cplus.deriv.noalias() += -1.0 * average_secant.deriv * thickness.value / 2.0 * thickness.value * transmission.value * exp(-1.0 * thickness.value * eigval.value(i));

            for(uint k = 0; k < numLayerDeriv; ++k) {
                Cplus.deriv(k + layerStart) += eigval.deriv(k, i) * Cplus.value * -1.0 * thickness.value;
                Cplus.deriv(k + layerStart) += eigval.deriv(k, i) * thickness.value / 2.0 * thickness.value * transmission.value * exp(-1.0 * thickness.value * eigval.value(i));

                Cplus.deriv(k + layerStart) += thickness.deriv(k) * transmission.value * exp(-1.0 * thickness.value * eigval.value(i)) * (1 - thickness.value * (average_secant.value - eigval.value(i)));
				Cplus.deriv(k + layerStart) += -1.0 * eigval.value(i) * Cplus.value * thickness.deriv(k);
            }
        }

        if(abs(average_secant.value + eigval.value(i)) > SKTRAN_DO_GREENS_EPS) {

            Cminus.value = transmission.value * (1 - exp(-thickness.value * average_secant.value) *
                                                     exp(-thickness.value * eigval.value(i))) /
                           (average_secant.value + eigval.value(i));

            Cminus.deriv.noalias() +=  (1 - exp(-thickness.value * average_secant.value) *
				exp(-thickness.value * eigval.value(i))) /
				(average_secant.value + eigval.value(i)) * transmission.deriv;
            Cminus.deriv.noalias() += average_secant.deriv / (average_secant.value + eigval.value(i)) *
                                      (transmission.value * thickness.value *
                                       exp(-thickness.value * (average_secant.value + eigval.value(i))) - Cminus.value);

            for (uint k = 0; k < numLayerDeriv; ++k) {
                Cminus.deriv(k + layerStart) += eigval.deriv(k, i) / (average_secant.value + eigval.value(i)) *
                                                (transmission.value * thickness.value *
                                                 exp(-thickness.value * (average_secant.value + eigval.value(i))) -
                                                 Cminus.value);
                Cminus.deriv(k + layerStart) += thickness.deriv(k) * transmission.value *
                                                exp(-1 * thickness.value * (average_secant.value + eigval.value(i)));
            }
        } else {
			// Second order taylor expansion of CMinus
            Cminus.value = transmission.value * thickness.value * (1 - thickness.value / 2 * (average_secant.value + eigval.value(i)));

            Cminus.deriv += transmission.deriv * thickness.value * (1 - thickness.value / 2 * (average_secant.value + eigval.value(i)));
            Cminus.deriv += average_secant.deriv * -1.0 * thickness.value / 2 * thickness.value * transmission.value;

            for(uint k = 0; k < numLayerDeriv; ++k) {
                Cminus.deriv(k + layerStart) += eigval.deriv(k, i) * transmission.value * thickness.value * thickness.value / 2 * -1.0;
                Cminus.deriv(k + layerStart) += thickness.deriv(k) * transmission.value * (1 - thickness.value * (average_secant.value + eigval.value(i)));
            }
        }

        for(uint j = 0; j < N * NSTOKES; ++j) {
            double negation = sasktran_disco::stokes_negation_factor<NSTOKES>(j);

            Gplus_top.value(j) += Aminus.value(i) * Cminus.value * homog_minus.value(h_start + j);
            Gminus_top.value(j) += negation * Aminus.value(i) * Cminus.value * homog_plus.value(h_start + j);

            Gplus_bottom.value(j) += Aplus.value(i) * Cplus.value * homog_plus.value(h_start + j);
            Gminus_bottom.value(j) += negation * Aplus.value(i) * Cplus.value * homog_minus.value(h_start + j);

            Gplus_top.deriv(Eigen::all, j).noalias() += Aminus.value(i) * homog_minus.value(h_start + j) * Cminus.deriv;
            Gminus_top.deriv(Eigen::all, j).noalias() += negation * Aminus.value(i) * homog_plus.value(h_start + j) * Cminus.deriv;
            Gplus_bottom.deriv(Eigen::all, j).noalias() += Aplus.value(i) * homog_plus.value(h_start + j) * Cplus.deriv;
            Gminus_bottom.deriv(Eigen::all, j).noalias() += negation * Aplus.value(i) * homog_minus.value(h_start + j) * Cplus.deriv;

            for(uint k = 0; k < numLayerDeriv; ++k) {
                Gplus_top.deriv(layerStart + k, j) += Aminus.deriv(k, i) * Cminus.value * homog_minus.value(h_start + j);
                Gminus_top.deriv(layerStart + k, j) += negation * Aminus.deriv(k, i) * Cminus.value * homog_plus.value(h_start + j);
                Gplus_bottom.deriv(layerStart + k, j) += Aplus.deriv(k, i) * Cplus.value * homog_plus.value(h_start + j);
                Gminus_bottom.deriv(layerStart + k, j) += negation * Aplus.deriv(k, i) * Cplus.value * homog_minus.value(h_start + j);

                Gplus_top.deriv(layerStart + k, j) += Aminus.value(i) * Cminus.value * homog_minus.deriv(k, h_start + j);
                Gminus_top.deriv(layerStart + k, j) += negation * Aminus.value(i) * Cminus.value * homog_plus.deriv(k, h_start + j);
                Gplus_bottom.deriv(layerStart + k, j) += Aplus.value(i) * Cplus.value * homog_plus.deriv(k, h_start + j);
                Gminus_bottom.deriv(layerStart + k, j) += negation * Aplus.value(i) * Cplus.value * homog_minus.deriv(k, h_start + j);
            }
        }
    }
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::RTESolver<NSTOKES, CNSTR>::solveBVP(AEOrder m)
{
	// Setup calculation environment
    // TODO: These should be moved to the cache
	la::BVPMatrix<NSTOKES>& mat = *m_cache.bvp_mat;  // Matrix that defines the system for the boundary value problem
	mat.setZero();
	
	Eigen::VectorXd& b = m_cache.bvp_b;  // Right hand side of the boundary value problem

	// Create derivative matricies that we need
	const auto& input_derivatives = m_layers.inputDerivatives();
	uint numDeriv = (uint)input_derivatives.numDerivative();

	// The derivatives of the BVP matrix are single dense blocks
    auto& d_mat = m_cache.d_mat;
    auto& d_b = m_cache.d_b;

    for (uint i = 0; i < numDeriv; ++i) {
        d_mat[i].setZero();
        d_b[i].setZero();
    }

	// Build the S.O.E for unperturbed system
	bvpTOACondition(m, 0, mat, d_mat);  // TOA Coefficients
	for (BoundaryIndex p = 1; p < this->M_NLYR; ++p) bvpContinuityCondition(m, p, mat, d_mat);  // Continuity coefficients
	bvpGroundCondition(m, this->M_NLYR, mat, d_mat); // Ground Coefficients
	uint loc = 0;
	bvpCouplingCondition_BC1(m, 0, loc, b, d_b); // TOA RHS
	for (BoundaryIndex p = 1; p < this->M_NLYR; ++p) bvpCouplingCondition_BC2(m, p, loc, b, d_b); // Continuity RHS
	bvpCouplingCondition_BC3(m, this->M_NLYR, loc, b, d_b); // Ground RHS

	lapack_int N = mat.N();
	lapack_int NCD = mat.NCD();
	lapack_int LDA = mat.LD();
	lapack_int errorcode;
    Eigen::VectorXi ipiv(this->M_NSTR * NSTOKES *  this->M_NLYR);

    if constexpr (CNSTR == 2 && NSTOKES==1 && SASKTRAN_DISCO_ENABLE_PENTADIAGONAL) {
        errorcode = la::dgbsv_pentadiagonal(N, 1, mat.data(), b.data(), N);
    }
    else {
    #ifdef SKTRAN_USE_ACCELERATE
        int one = 1;
        dgbsv_(&N, &NCD, &NCD, &one, mat.data(), &LDA, ipiv.data(), b.data(), &N, &errorcode);
    #else
        errorcode = LAPACKE_dgbsv(LAPACK_COL_MAJOR, N, NCD, NCD, 1, mat.data(), LDA, ipiv.data(), b.data(), N);
    #endif
    }



	// We we failed then return immediately
	if (errorcode != 0) {
		if (errorcode < 0)
			throw InternalRuntimeError("LAPACKE_dgbsv had an illegal argument in "
				"sasktran_disco::RTESolver::SolveBVP");
		else
			throw InternalRuntimeError("LAPACKE_dgbsv couldn't solve BVP since "
				"the coefficient matrix was singular.");
	}

	// Copy solutions to general solutions
	for (LayerIndex p = 0; p < this->M_NLYR; ++p) {
		auto& solution = m_layers[p].solution(m);
		for (SolutionIndex j = 0; j < this->M_NSTR*NSTOKES / 2; ++j) {
			solution.boundary.L_coeffs.value(j) = b[p * this->M_NSTR*NSTOKES + j];
			solution.boundary.M_coeffs.value(j) = b[p * this->M_NSTR*NSTOKES + j + this->M_NSTR / 2 * NSTOKES];
		}
	}

    if(numDeriv == 0) {
        // can immediately return
        return;
    }



	Eigen::VectorXd& temp = m_cache.bvp_temp;
	// Compute the linearizations of the boundary coefficients
	// Create an RHS matrix to store all of the RHS

	Eigen::MatrixXd& rhs = m_cache.bvp_d_rhs;
	for (uint i = 0; i < numDeriv; ++i) {
        /*
		d_mat[i].multiplyVector(b, temp);
		rhs(Eigen::all, i) = d_b[i] - temp;
         */
        d_mat[i].assign_rhs_d_bvp(rhs, i, d_b[i], b);
	}

    if constexpr (NSTOKES ==1 && CNSTR == 2 && SASKTRAN_DISCO_ENABLE_PENTADIAGONAL) {
        errorcode = la::dgbsv_pentadiagonal(N, numDeriv, mat.data(), rhs.data(), N);
    }
    else {
        // Solve using the same LHS as above
    #ifdef SKTRAN_USE_ACCELERATE
        char trans = 'N';
        int intnderiv = numDeriv;
        dgbtrs_(&trans, &N, &NCD, &NCD, &intnderiv, mat.data(), &LDA, ipiv.data(), rhs.data(), &N, &errorcode);
    #else
        errorcode = LAPACKE_dgbtrs(LAPACK_COL_MAJOR, 'N', N, NCD, NCD, numDeriv, mat.data(), LDA, ipiv.data(), rhs.data(), N);
    #endif
    }
	if (errorcode != 0) {
		if (errorcode < 0)
		{
			std::string error = "LAPACKE_dgbtrs had an illegal argument in sasktran_disco::RTESolver::SolveBVP when solving derivative, error code: " + std::to_string(errorcode);
			throw InternalRuntimeError(error.c_str());
		}
		else
			throw InternalRuntimeError("LAPACKE_dgetrs failed to compute particular solution derivative "
				"as the coefficient matrix was singular in "
				"sasktran_disco::RTESolver::SolveParticualr");
	}
    for (LayerIndex p = 0; p < this->M_NLYR; ++p) {
        auto& solution = m_layers[p].solution(m);
        int l_offset = p * this->M_NSTR * NSTOKES;
        int m_offset = l_offset + this->M_NSTR * NSTOKES / 2;

        for (uint i = 0; i < numDeriv; ++i) {
		    // Store the results in the solution derivative
			for (SolutionIndex k = 0; k < this->M_NSTR / 2 * NSTOKES; ++k) {
				solution.boundary.L_coeffs.deriv(i, k) = rhs(l_offset + k, i);
				solution.boundary.M_coeffs.deriv(i, k) = rhs(m_offset + k, i);
			}
		}
	}
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::RTESolver<NSTOKES, CNSTR>::bvpTOACondition(AEOrder m, BoundaryIndex p, BVPMatrix& A, VectorDim1<BVPMatrixDenseBlock<NSTOKES>>& d_A) const
{
	// Eq (29)

	// Setup enviroment
	assert(p == 0);
	typename BVPMatrix::Block mat = A.getBlock(p);
	const auto& solution = m_layers.top().solution(m).value;
	const OpticalLayer<NSTOKES, CNSTR>& layer = m_layers[p];

	const auto& input_deriv = m_layers.inputDerivatives().layerDerivatives();
	uint numDeriv = (uint)m_layers.inputDerivatives().numDerivativeLayer(layer.index());
	uint derivStartIndex = (uint)m_layers.inputDerivatives().layerStartIndex(layer.index());
	const auto& layer_deriv = m_layers.top().solution(m).value;

	const uint N = this->M_NSTR / 2;

	for (StreamIndex i = 0; i < N*NSTOKES; ++i) {
		for (SolutionIndex j = 0; j < N*NSTOKES; ++j) {
            mat(i, j) = solution.homog_plus(i, j);
            mat(i, j + N*NSTOKES) = solution.homog_minus(i, j) * layer.streamTransmittance(Location::INSIDE, m, j);

			// This condition does not contain any cross derivatives
			for (uint l = 0; l < numDeriv; ++l) {
				auto& d_mat = d_A[l + derivStartIndex].layer();

				d_mat(i, j) = solution.d_homog_plus(i, j, l);
				d_mat(i, j + N*NSTOKES) = solution.d_homog_minus(i, j, l) * layer.streamTransmittance(Location::INSIDE, m, j) +
					solution.homog_minus(i, j) * layer.d_streamTransmittance(Location::INSIDE, m, j, l, input_deriv[l + derivStartIndex]);
			}
		}
	}
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::RTESolver<NSTOKES, CNSTR>::bvpContinuityCondition(AEOrder m, BoundaryIndex p, BVPMatrix& A, VectorDim1<BVPMatrixDenseBlock<NSTOKES>>& d_A) const
{
	// Eq (31)
	// Setup enviroment
	assert(p != 0 && p != this->M_NLYR);
	const uint N = this->M_NSTR / 2;

	typename BVPMatrix::Block mat = A.getBlock(p);
	const OpticalLayer<NSTOKES, CNSTR>& layer = m_layers[p];
	const auto& input_deriv = m_layers.inputDerivatives().layerDerivatives();
	struct Layer
	{
		const RTEGeneralSolution<NSTOKES, CNSTR>& solution;
		const OpticalLayer<NSTOKES, CNSTR>& optical;
		uint numDeriv;
		uint layerDerivStart;
	};
	Layer upper = { m_layers[p - 1].solution(m).value, m_layers[p - 1], (uint)m_layers.inputDerivatives().numDerivativeLayer(p - 1), (uint)m_layers.inputDerivatives().layerStartIndex(p - 1) };
	Layer lower = { m_layers[p].solution(m).value, m_layers[p], (uint)m_layers.inputDerivatives().numDerivativeLayer(p), (uint)m_layers.inputDerivatives().layerStartIndex(p) };

	// Build the LHS continuity condition for coupling layers p-1 and p
	for (StreamIndex i = 0; i < N*NSTOKES; ++i) {
		for (SolutionIndex j = 0, j_top = 0, j_bot = this->M_NSTR*NSTOKES; j < N*NSTOKES; ++j, ++j_top, ++j_bot) {
			mat(i + N*NSTOKES, j_top) = upper.solution.homog_plus(i, j) * upper.optical.streamTransmittance(Location::INSIDE, m, j);
			mat(i + N*NSTOKES, j_bot) = -lower.solution.homog_plus(i, j);

			mat(i, j_top) = sasktran_disco::stokes_negation_factor<NSTOKES>(i) * upper.solution.homog_minus(i, j) * upper.optical.streamTransmittance(Location::INSIDE, m, j);
			mat(i, j_bot) = -sasktran_disco::stokes_negation_factor<NSTOKES>(i) * lower.solution.homog_minus(i, j);

			// We still have no cross derivatives, but we have to treat the derivatives for the upper and lower layers separately
			for (uint k = 0; k < upper.numDeriv; ++k) {
				auto& d_mat = d_A[upper.layerDerivStart + k].upper();

				d_mat(i + N*NSTOKES, j_top) = upper.solution.d_homog_plus(i, j, k) * upper.optical.streamTransmittance(Location::INSIDE, m, j) +
					upper.solution.homog_plus(i, j) * upper.optical.d_streamTransmittance(Location::INSIDE, m, j, k, input_deriv[upper.layerDerivStart + k]);

				d_mat(i, j_top) = sasktran_disco::stokes_negation_factor<NSTOKES>(i) * (upper.solution.d_homog_minus(i, j, k) * upper.optical.streamTransmittance(Location::INSIDE, m, j) +
					upper.solution.homog_minus(i, j) * upper.optical.d_streamTransmittance(Location::INSIDE, m, j, k, input_deriv[upper.layerDerivStart + k]));
			}
			for (uint k = 0; k < lower.numDeriv; ++k) {
				auto& d_mat = d_A[lower.layerDerivStart + k].layer();

				d_mat(i + N*NSTOKES, j_bot) = -lower.solution.d_homog_plus(i, j, k);
				d_mat(i, j_bot) = -sasktran_disco::stokes_negation_factor<NSTOKES>(i) * lower.solution.d_homog_minus(i, j, k);
			}

		}
		for (SolutionIndex j = 0, j_top = this->M_NSTR / 2 * NSTOKES, j_bot = this->M_NSTR*NSTOKES + this->M_NSTR / 2 * NSTOKES; j < N*NSTOKES; ++j, ++j_top, ++j_bot) {
			mat(i + N*NSTOKES, j_top) = upper.solution.homog_minus(i, j);
			mat(i + N*NSTOKES, j_bot) = -lower.solution.homog_minus(i, j) * lower.optical.streamTransmittance(Location::INSIDE, m, j);

			mat(i, j_top) = sasktran_disco::stokes_negation_factor<NSTOKES>(i) * upper.solution.homog_plus(i, j);
			mat(i, j_bot) = -sasktran_disco::stokes_negation_factor<NSTOKES>(i) * lower.solution.homog_plus(i, j) * lower.optical.streamTransmittance(Location::INSIDE, m, j);

			// We still have no cross derivatives, but we have to treat the derivatives for the upper and lower layers separately
			for (uint k = 0; k < lower.numDeriv; ++k) {
				auto& d_mat = d_A[lower.layerDerivStart + k].layer();

				d_mat(i + N*NSTOKES, j_bot) = -lower.solution.d_homog_minus(i, j, k) * lower.optical.streamTransmittance(Location::INSIDE, m, j) -
					lower.solution.homog_minus(i, j) * lower.optical.d_streamTransmittance(Location::INSIDE, m, j, k, input_deriv[lower.layerDerivStart + k]);

				d_mat(i, j_bot) = -sasktran_disco::stokes_negation_factor<NSTOKES>(i) * (lower.solution.d_homog_plus(i, j, k) * lower.optical.streamTransmittance(Location::INSIDE, m, j) +
					lower.solution.homog_plus(i, j) * lower.optical.d_streamTransmittance(Location::INSIDE, m, j, k, input_deriv[lower.layerDerivStart + k]));
			}
			for (uint k = 0; k < upper.numDeriv; ++k) {
				auto& d_mat = d_A[upper.layerDerivStart + k].upper();

				d_mat(i + N*NSTOKES, j_top) = upper.solution.d_homog_minus(i, j, k);
				d_mat(i, j_top) = sasktran_disco::stokes_negation_factor<NSTOKES>(i) * upper.solution.d_homog_plus(i, j, k);
			}

		}
	}
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::RTESolver<NSTOKES, CNSTR>::bvpGroundCondition(AEOrder m, BoundaryIndex p, BVPMatrix& A, VectorDim1<BVPMatrixDenseBlock<NSTOKES>>& d_A) const
{

	// Setup enviroment
	assert(p == this->M_NLYR);
	const uint N = this->M_NSTR / 2;
	typename BVPMatrix::Block mat = A.getBlock(p);
	const auto& solution = m_layers.bottom().solution(m).value;
	const auto& input_deriv = m_layers.inputDerivatives().layerDerivatives();

	const OpticalLayer<NSTOKES, CNSTR>& layer = m_layers[p - 1];
	uint layerDerivStart = (uint)m_layers.inputDerivatives().layerStartIndex(layer.index());
	uint numDeriv = (uint)m_layers.inputDerivatives().numDerivativeLayer(layer.index());

	double sum_vm, sum_vp, d_sum;

	// Build the LHS bottom layer to ground coupling condition
	for(StreamIndex i = 0; i < N*NSTOKES; ++i) {
		for(SolutionIndex j = 0; j < N*NSTOKES; ++j) {
			sum_vm = v_minus(m, layer, i, j);
			mat(i, j) = sasktran_disco::stokes_negation_factor<NSTOKES>(i) * sum_vm * layer.streamTransmittance(Location::INSIDE, m, j);

			sum_vp = v_plus(m, layer, i, j);
			mat(i, j + N*NSTOKES) = sasktran_disco::stokes_negation_factor<NSTOKES>(i) * sum_vp;

			for (uint k = 0; k < numDeriv; ++k) {
				auto& d_mat = d_A[layerDerivStart + k].upper();

				d_sum = d_v_minus(m, layer, i, j, k, input_deriv[layerDerivStart + k]);
				d_mat(i, j) = sasktran_disco::stokes_negation_factor<NSTOKES>(i) * (d_sum * layer.streamTransmittance(Location::INSIDE, m, j) + sum_vm * layer.d_streamTransmittance(Location::INSIDE, m, j, k, input_deriv[layerDerivStart + k]));

				d_sum = d_v_plus(m, layer, i, j, k, input_deriv[layerDerivStart + k]);
				d_mat(i, j + N*NSTOKES) = sasktran_disco::stokes_negation_factor<NSTOKES>(i) * d_sum;
			}
			
		}
	}
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::RTESolver<NSTOKES, CNSTR>::bvpCouplingCondition_BC1(AEOrder m, BoundaryIndex p, uint& loc, Eigen::VectorXd& b, VectorDim1<Eigen::VectorXd>& d_b) const
{

	// Setup enviroment
	assert(p == 0);
	const auto& solution = m_layers.top().solution(m).value;
	const uint N = this->M_NSTR / 2;

	uint numDeriv = (uint)m_layers.inputDerivatives().numDerivative();

	// Build TOA condition. See (29)
	for (StreamIndex i = 0; i < N*NSTOKES; ++i) {
        b[loc] = -solution.dual_Gplus_top().value(i);

        for (uint k = 0; k < numDeriv; ++k) {
            d_b[k][loc] = -solution.dual_Gplus_top().deriv(k, i);
        }
		loc++;
	}
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::RTESolver<NSTOKES, CNSTR>::bvpCouplingCondition_BC2(AEOrder m, BoundaryIndex p, uint& loc, Eigen::VectorXd& b, VectorDim1<Eigen::VectorXd>& d_b) const
{
	// Setup enviroment
	assert(p != 0 && p != this->M_NLYR);
	const uint N = this->M_NSTR / 2;
	const auto& input_deriv = m_layers.inputDerivatives().layerDerivatives();
	const auto numderiv = m_layers.inputDerivatives().numDerivative();

	struct Layer
	{
		const RTEGeneralSolution<NSTOKES, CNSTR>& solution;
		const OpticalLayer<NSTOKES, CNSTR>& optical;
	};
	Layer upper = { m_layers[p - 1].solution(m).value, m_layers[p - 1] };
	Layer lower = { m_layers[p].solution(m).value, m_layers[p] };


	// Build coupling condition. See eq 31
	for (StreamIndex i = 0; i < N*NSTOKES; ++i) {
		double u_p_p = upper.solution.particular_plus(i);
		double u_p_m = upper.solution.particular_minus(i);
		double u_bt = upper.optical.beamTransmittance(Location::FLOOR);

		double l_p_p = lower.solution.particular_plus(i);
		double l_p_m = lower.solution.particular_minus(i);
		double l_bt = lower.optical.beamTransmittance(Location::CEILING);

        b[loc] = -upper.solution.dual_Gminus_bottom().value(i) + lower.solution.dual_Gminus_top().value(i);
        b[loc + N*NSTOKES] = -upper.solution.dual_Gplus_bottom().value(i) + lower.solution.dual_Gplus_top().value(i);

        for (uint j = 0; j < input_deriv.size(); ++j) {
            d_b[j][loc + N*NSTOKES] = lower.solution.dual_Gplus_top().deriv(j, i) - upper.solution.dual_Gplus_bottom().deriv(j, i);
            d_b[j][loc] = lower.solution.dual_Gminus_top().deriv(j, i) - upper.solution.dual_Gminus_bottom().deriv(j, i);
        }

		loc++;
	}
	loc += N*NSTOKES;
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::RTESolver<NSTOKES, CNSTR>::bvpCouplingCondition_BC3(AEOrder m, BoundaryIndex p, uint& loc, Eigen::VectorXd& b, VectorDim1<Eigen::VectorXd>& d_b) const
{
	// Configure
	assert(p == this->M_NLYR);
	RTEGeneralSolution<NSTOKES, CNSTR>& solution = m_layers.bottom().solution(m).value;
	const OpticalLayer<NSTOKES, CNSTR>& layer = m_layers[p - 1];
	const auto& input_deriv = m_layers.inputDerivatives().layerDerivatives();

	// Build ground condition. See eq (36)
	for (StreamIndex i = 0; i < this->M_NSTR / 2 * NSTOKES; ++i) {
		b[loc] = ground_direct_sun(m, layer, i) -u_minus(m, layer, i);

		// We have cross derivatives
		for (uint j = 0; j < input_deriv.size(); ++j) {
			d_b[j][loc] = d_ground_direct_sun(m, layer, i, input_deriv[j], j) - d_u_minus(m, layer, i, j, input_deriv[j]);
		}

		loc++;
	}
}

SASKTRAN_DISCO_INSTANTIATE_TEMPLATE(sasktran_disco::RTESolver);
