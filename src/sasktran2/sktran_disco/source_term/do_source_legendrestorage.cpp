#include "sasktran2/do_source.h"


namespace sasktran2 {
    template<int NSTOKES, int CNSTR>
    LegendrePhaseStorage<NSTOKES, CNSTR>::LegendrePhaseStorage(int nstr) : nstr(nstr) {
        // Storage is linear in stream, azimuth expansion order
        storage.resize(NSTOKES, nstr*nstr);
    }

    template<int NSTOKES, int CNSTR>
    int LegendrePhaseStorage<NSTOKES, CNSTR>::linear_index(int m, int l) const {
        return m*nstr + l;
    }

    template<int NSTOKES, int CNSTR>
    void LegendrePhaseStorage<NSTOKES, CNSTR>::fill(double coszen) {
        double theta = acos(coszen);

        for(int m = 0; m < nstr; ++m) {
            auto calculator = sasktran2::math::WignerDCalculator(m, 0);
            for(int l = 0; l < nstr; ++l) {
                int idx = linear_index(m, l);

                storage(0, idx) = calculator.d(theta, l);
            }

            if constexpr(NSTOKES > 1) {
                auto calculatorneg = sasktran_disco::WignerDCalculator(m, -2);
                auto calculatorpos = sasktran_disco::WignerDCalculator(m, 2);
                for(int l = 0; l < nstr; ++l) {
                    int idx = linear_index(m, l);

                    storage(1, idx) = -0.5 * (calculatorpos.d(theta, l) + calculatorneg.d(theta, l));
                    storage(2, idx) = -0.5 * (calculatorpos.d(theta, l) - calculatorneg.d(theta, l));
                }
            }
        }
    }

    template<int NSTOKES, int CNSTR>
    void LegendrePhaseStorage<NSTOKES, CNSTR>::set_phase_container(
            std::vector<sasktran_disco::LegendrePhaseContainer<NSTOKES>> &container, int m) const {
        for(int l = 0; l < nstr; ++l) {
            int idx = linear_index(m, l);

            if constexpr(NSTOKES == 3) {
                container[l].values = storage(Eigen::all, idx);
            }
            if constexpr(NSTOKES == 1) {
                container[l].value = storage(0, idx);
            }
        }
    }

    SASKTRAN_DISCO_INSTANTIATE_TEMPLATE_STRUCT(LegendrePhaseStorage);
}