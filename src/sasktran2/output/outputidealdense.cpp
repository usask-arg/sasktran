#include <sasktran2/output.h>


namespace sasktran2 {
    template<int NSTOKES>
    void OutputIdealDense<NSTOKES>::resize(int nlos, int nwavel, int nderiv) {
        Output<NSTOKES>::resize(nlos, nwavel, nderiv);
        m_radiance.resize(NSTOKES * nwavel * nlos, nderiv, false);
    }

    template<int NSTOKES>
    void
    OutputIdealDense<NSTOKES>::assign(const sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES> &radiance,
                                      int losidx, int wavelidx) {
        int linear_index = NSTOKES * this->m_nlos * wavelidx + NSTOKES * losidx;

        for(int i = 0; i < NSTOKES; ++i) {
            m_radiance.value(linear_index + i) = radiance.value(i);
            m_radiance.deriv(linear_index + i, Eigen::all) = radiance.deriv(i, Eigen::all);
        }
    }

    template class OutputIdealDense<1>;
    template class OutputIdealDense<3>;

}