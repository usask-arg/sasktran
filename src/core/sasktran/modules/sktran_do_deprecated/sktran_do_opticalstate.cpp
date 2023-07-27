#include "modules/sktran_do_deprecated/include/sktran_do.h"
#include "modules/sktran_do_deprecated/include/sktran_do_opticalstate.h"

template<>
void sktran_do_detail::OpticalState<1>::fill_species_scattering_entry(size_t speciesidx, sktran_do_detail::uint wavelidx, sktran_do_detail::uint heightindex, skOpticalProperties* optprop){
    int numlegendre;

    auto& entry = m_wavelengthentry[wavelidx];

    optprop->LegendreCoefficientsP11(1e7 / m_wavelengths[wavelidx], &entry.speciesa1[speciesidx](0, heightindex), m_numstr + 1, numlegendre);

    // Delta scaling fraction
    double f;
    if (numlegendre == m_numstr + 1) {
        // Have delta scaling
        f = entry.speciesa1[speciesidx](m_numstr, heightindex) / (2 * m_numstr + 1);
    }
    else {
        f = 0.0;
    }
    double ssa = entry.speciesscatxs[speciesidx](heightindex) / entry.speciesxs[speciesidx](heightindex);
    double scaledssa = ssa * (1 - f) / (1 - ssa * f);
    entry.speciesxs[speciesidx](heightindex) *= 1 - ssa * f;
    entry.speciesscatxs[speciesidx](heightindex) = entry.speciesxs[speciesidx](heightindex) * scaledssa;

    for (size_t k = 0; k < numlegendre; ++k) {
        entry.speciesa1[speciesidx](k, heightindex) = (entry.speciesa1[speciesidx](k, heightindex) - f * (2 * k + 1)) / (1 - f);
    }

    /*
    for (size_t k = 0; k < los.size(); ++k) {
        // Scale the phase function so that single scatter is the same
        entry.speciesphase[i](j, k) /= (1 - f);
    }
     */


}

template<int NSTOKES, int CNSTR>
void sktran_do_detail::OpticalState<NSTOKES, CNSTR>::fill_species_scattering_entry(size_t speciesidx, sktran_do_detail::uint wavelidx, sktran_do_detail::uint heightindex, skOpticalProperties* optprop){
    int numlegendre;

    auto& entry = m_wavelengthentry[wavelidx];

    optprop->LegendreCoefficientsPolarized(1e7 / m_wavelengths[wavelidx],
                                           &entry.speciesa1[speciesidx](0, heightindex),
                                           &entry.speciesa2[speciesidx](0, heightindex),
                                           &entry.speciesa3[speciesidx](0, heightindex),
                                           &entry.speciesa4[speciesidx](0, heightindex),
                                           &entry.speciesb1[speciesidx](0, heightindex),
                                           &entry.speciesb2[speciesidx](0, heightindex),
                                           m_numstr + 1,
                                           numlegendre);

    // Delta scaling fraction
    double f;
    if (numlegendre == m_numstr + 1) {
        // Have delta scaling
        f = entry.speciesa1[speciesidx](m_numstr, heightindex) / (2 * m_numstr + 1);
    }
    else {
        f = 0.0;
    }
    double ssa = entry.speciesscatxs[speciesidx](heightindex) / entry.speciesxs[speciesidx](heightindex);
    double scaledssa = ssa * (1 - f) / (1 - ssa * f);
    entry.speciesxs[speciesidx](heightindex) *= 1 - ssa * f;
    entry.speciesscatxs[speciesidx](heightindex) = entry.speciesxs[speciesidx](heightindex) * scaledssa;

    for (size_t k = 0; k < numlegendre; ++k) {
        entry.speciesa1[speciesidx](k, heightindex) = (entry.speciesa1[speciesidx](k, heightindex) - f * (2 * k + 1)) / (1 - f);
        entry.speciesa2[speciesidx](k, heightindex) = (entry.speciesa2[speciesidx](k, heightindex) - f * (2 * k + 1)) / (1 - f);
        entry.speciesa3[speciesidx](k, heightindex) = (entry.speciesa3[speciesidx](k, heightindex) - f * (2 * k + 1)) / (1 - f);
        entry.speciesa4[speciesidx](k, heightindex) = (entry.speciesa4[speciesidx](k, heightindex) - f * (2 * k + 1)) / (1 - f);

        entry.speciesb1[speciesidx](k, heightindex) = (entry.speciesb1[speciesidx](k, heightindex)) / (1 - f);
        entry.speciesb2[speciesidx](k, heightindex) = (entry.speciesb2[speciesidx](k, heightindex)) / (1 - f);
    }
}

template <>
void sktran_do_detail::OpticalState<1>::assign_legendre_quantities(Eigen::Matrix<double,
        Eigen::Dynamic, 6>& bl,
        const OpticalStateWavelengthEntry<1>& entry,
        double altitude_weight, double cs_scat,
        const std::pair<size_t, size_t>& indicies,
        const std::pair<double, double>& weights ) const {

    bl(Eigen::all, 0).noalias() += altitude_weight * cs_scat* ( entry.totala1(Eigen::all, indicies.first) * weights.first + entry.totala1(Eigen::all, indicies.second) * weights.second );
}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::OpticalState<NSTOKES, CNSTR>::assign_legendre_quantities(Eigen::Matrix<double,
        Eigen::Dynamic, 6>& bl,
        const OpticalStateWavelengthEntry<NSTOKES>& entry,
        double altitude_weight, double cs_scat,
        const std::pair<size_t, size_t>& indicies,
        const std::pair<double, double>& weights ) const {

    bl(Eigen::all, 0).noalias() += altitude_weight * cs_scat* ( entry.totala1(Eigen::all, indicies.first) * weights.first + entry.totala1(Eigen::all, indicies.second) * weights.second );
    bl(Eigen::all, 1).noalias() += altitude_weight * cs_scat* ( entry.totala2(Eigen::all, indicies.first) * weights.first + entry.totala2(Eigen::all, indicies.second) * weights.second );
    bl(Eigen::all, 2).noalias() += altitude_weight * cs_scat* ( entry.totala3(Eigen::all, indicies.first) * weights.first + entry.totala3(Eigen::all, indicies.second) * weights.second );
    bl(Eigen::all, 3).noalias() += altitude_weight * cs_scat* ( entry.totala4(Eigen::all, indicies.first) * weights.first + entry.totala4(Eigen::all, indicies.second) * weights.second );
    bl(Eigen::all, 4).noalias() += altitude_weight * cs_scat* ( entry.totalb1(Eigen::all, indicies.first) * weights.first + entry.totalb1(Eigen::all, indicies.second) * weights.second );
    bl(Eigen::all, 5).noalias() += altitude_weight * cs_scat* ( entry.totalb2(Eigen::all, indicies.first) * weights.first + entry.totalb2(Eigen::all, indicies.second) * weights.second );
}

template <>
void sktran_do_detail::OpticalState<1>::assign_species_legendre_quantities(
        Eigen::Matrix<double, Eigen::Dynamic, 6> &bl, const OpticalStateWavelengthEntry<1> &entry,
        int species_idx, const std::pair<size_t, size_t> &indicies, const std::pair<double, double> &weights,
        SpeciesType type) const {
    auto seq = Eigen::seq(0, bl.rows() - 1);

    if(type == ScatteringHeightDependent) {
        bl(Eigen::all, 0).noalias() = ( entry.speciesa1[species_idx](seq, indicies.first) * weights.first + entry.speciesa1[species_idx](seq, indicies.second) * weights.second );
    } else if(type == ScatteringHeightIndependent) {
        bl(Eigen::all, 0).noalias() = entry.speciesa1[species_idx](seq, indicies.first);
    }
}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::OpticalState<NSTOKES, CNSTR>::assign_species_legendre_quantities(
        Eigen::Matrix<double, Eigen::Dynamic, 6> &bl, const OpticalStateWavelengthEntry<NSTOKES> &entry,
        int species_idx, const std::pair<size_t, size_t> &indicies, const std::pair<double, double> &weights,
        SpeciesType type) const {
    auto seq = Eigen::seq(0, bl.rows() - 1);

    if(type == ScatteringHeightDependent) {
        bl(Eigen::all, 0).noalias() = ( entry.speciesa1[species_idx](seq, indicies.first) * weights.first + entry.speciesa1[species_idx](seq, indicies.second) * weights.second );
        bl(Eigen::all, 1).noalias() = ( entry.speciesa2[species_idx](seq, indicies.first) * weights.first + entry.speciesa2[species_idx](seq, indicies.second) * weights.second );
        bl(Eigen::all, 2).noalias() = ( entry.speciesa3[species_idx](seq, indicies.first) * weights.first + entry.speciesa3[species_idx](seq, indicies.second) * weights.second );
        bl(Eigen::all, 3).noalias() = ( entry.speciesa4[species_idx](seq, indicies.first) * weights.first + entry.speciesa4[species_idx](seq, indicies.second) * weights.second );
        bl(Eigen::all, 4).noalias() = ( entry.speciesb1[species_idx](seq, indicies.first) * weights.first + entry.speciesb1[species_idx](seq, indicies.second) * weights.second );
        bl(Eigen::all, 5).noalias() = ( entry.speciesb2[species_idx](seq, indicies.first) * weights.first + entry.speciesb2[species_idx](seq, indicies.second) * weights.second );
    } else if(type == ScatteringHeightIndependent) {
        bl(Eigen::all, 0).noalias() = entry.speciesa1[species_idx](seq, indicies.first);
        bl(Eigen::all, 1).noalias() = entry.speciesa2[species_idx](seq, indicies.first);
        bl(Eigen::all, 2).noalias() = entry.speciesa3[species_idx](seq, indicies.first);
        bl(Eigen::all, 3).noalias() = entry.speciesa4[species_idx](seq, indicies.first);
        bl(Eigen::all, 4).noalias() = entry.speciesb1[species_idx](seq, indicies.first);
        bl(Eigen::all, 5).noalias() = entry.speciesb2[species_idx](seq, indicies.first);
    }
}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::OpticalStateWavelengthEntry<NSTOKES, CNSTR>::copy_height_independent_species(size_t speciesidx) {

    speciesxs[speciesidx].setConstant(speciesxs[speciesidx](0));
    speciesscatxs[speciesidx].setConstant(speciesscatxs[speciesidx](0));

    for(int i = 0; i < speciesa1[speciesidx].cols(); ++i) {
        speciesa1[speciesidx](Eigen::all, i) = speciesa1[speciesidx](Eigen::all, 0);
        speciesa2[speciesidx](Eigen::all, i) = speciesa2[speciesidx](Eigen::all, 0);
        speciesa3[speciesidx](Eigen::all, i) = speciesa3[speciesidx](Eigen::all, 0);
        speciesa4[speciesidx](Eigen::all, i) = speciesa4[speciesidx](Eigen::all, 0);
        speciesb1[speciesidx](Eigen::all, i) = speciesb1[speciesidx](Eigen::all, 0);
        speciesb2[speciesidx](Eigen::all, i) = speciesb2[speciesidx](Eigen::all, 0);
    }

    for(int i = 0; i < speciesphase[speciesidx].rows(); ++i) {
        speciesphase[speciesidx](i, Eigen::all) = speciesphase[speciesidx](0, Eigen::all);
    }
}

void sktran_do_detail::OpticalStateWavelengthEntry<1>::copy_height_independent_species(size_t speciesidx) {
    speciesxs[speciesidx].setConstant(speciesxs[speciesidx](0));
    speciesscatxs[speciesidx].setConstant(speciesscatxs[speciesidx](0));

    for(int i = 0; i < speciesa1[speciesidx].cols(); ++i) {
        speciesa1[speciesidx](Eigen::all, i) = speciesa1[speciesidx](Eigen::all, 0);
    }
    for(int i = 0; i < speciesphase[speciesidx].rows(); ++i) {
        speciesphase[speciesidx](i, Eigen::all) = speciesphase[speciesidx](0, Eigen::all);
    }
}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::OpticalStateWavelengthEntry<NSTOKES, CNSTR>::sum_totals(const nx2dArray<double>& speciesnd, std::vector<SKTRAN_AtmosphericOpticalStateEntry_V21>& allspecies) {
    int numalt = speciesxs[0].rows();
    int numstr = speciesa1[0].rows() - 1;

    for (int i = 0; i < numalt; ++i) {
        for (size_t j = 0; j < allspecies.size(); ++j) {
            totalext(i) += speciesnd.At(j, i) * speciesxs[j](i) * 100; // Convert /cm to /m
            totalscatext(i) += speciesnd.At(j, i) * speciesscatxs[j](i) * 100;

            if (allspecies[j].ParticleOpticalProps()->IsScatterer()) {
                for (size_t k = 0; k < numstr; ++k) {
                    totala1(k, i) += speciesscatxs[j](i) * 100 * speciesa1[j](k, i) * speciesnd.At(j, i);
                    totala2(k, i) += speciesscatxs[j](i) * 100 * speciesa2[j](k, i) * speciesnd.At(j, i);
                    totala3(k, i) += speciesscatxs[j](i) * 100 * speciesa3[j](k, i) * speciesnd.At(j, i);
                    totala4(k, i) += speciesscatxs[j](i) * 100 * speciesa4[j](k, i) * speciesnd.At(j, i);
                    totalb1(k, i) += speciesscatxs[j](i) * 100 * speciesb1[j](k, i) * speciesnd.At(j, i);
                    totalb2(k, i) += speciesscatxs[j](i) * 100 * speciesb2[j](k, i) * speciesnd.At(j, i);
                }

                for (size_t k = 0; k < totalphase.cols(); ++k) {
                    totalphase(i, k) += speciesscatxs[j](i) * 100 * speciesphase[j](i, k) * speciesnd.At(j, i);
                }
            }
        }
        if (totalscatext(i) != 0.0) {
            // If it is 0 then totals are already 0, don't need to normalize
            totala1(Eigen::all, i) /= totalscatext(i);
            totala2( Eigen::all, i) /= totalscatext(i);
            totala3( Eigen::all, i) /= totalscatext(i);
            totala4( Eigen::all, i) /= totalscatext(i);
            totalb1( Eigen::all, i) /= totalscatext(i);
            totalb2( Eigen::all, i) /= totalscatext(i);

            totalphase(i, Eigen::all) /= totalscatext(i);
        }
        totalphase(i, Eigen::all) /= (4 * nxmath::Pi);
    }
}

void sktran_do_detail::OpticalStateWavelengthEntry<1>::sum_totals(const nx2dArray<double>& speciesnd, std::vector<SKTRAN_AtmosphericOpticalStateEntry_V21>& allspecies) {
    int numalt = speciesxs[0].rows();
    int numstr = speciesa1[0].rows() - 1;

    for (int i = 0; i < numalt; ++i) {
        for (size_t j = 0; j < allspecies.size(); ++j) {
            totalext(i) += speciesnd.At(j, i) * speciesxs[j](i) * 100; // Convert /cm to /m
            totalscatext(i) += speciesnd.At(j, i) * speciesscatxs[j](i) * 100;

            if (allspecies[j].ParticleOpticalProps()->IsScatterer()) {
                for (size_t k = 0; k < numstr; ++k) {
                    totala1(k, i) += speciesscatxs[j](i) * 100 * speciesa1[j](k, i) * speciesnd.At(j, i);
                }

                for (size_t k = 0; k < totalphase.cols(); ++k) {
                    totalphase(i, k) += speciesscatxs[j](i) * 100 * speciesphase[j](i, k) * speciesnd.At(j, i);
                }
            }
        }
        if (totalscatext(i) != 0.0) {
            // If it is 0 then totals are already 0, don't need to normalize
            totala1( Eigen::all, i) /= totalscatext(i);
            totalphase(i, Eigen::all) /= totalscatext(i);
        }
        totalphase(i, Eigen::all) /= (4 * nxmath::Pi);
    }
}


namespace sktran_do_detail {
    template class OpticalState<1>;
    template class OpticalState<3>;
    template class OpticalState<4>;
}

