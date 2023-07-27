#pragma once
#include "modules/sktran_do_deprecated/include/sktran_do.h"
#include "omp.h"

// OpticalState is the layer inbetween SASKTRANIF and DO, it stores all of the relevant optical information
// on an altitude grid that is then later translated to layer quantities.  In spherical mode the OpticalState is
// used directly
// Optical information can be stored at multiple wavelengths to use the faster multiwavelength functions in SASKTRANIF

namespace sktran_do_detail {

    // Common interface for components of the code where NSTOKES is not known and not needed, mostly configuration and setup
    class OpticalStateInterface {
    public:
        virtual void configure(const Eigen::VectorXd& internal_altitudes, uint numstr) = 0;

        virtual void fill_tables(double wavelength_nm, const GEODETIC_INSTANT& refpt, SASKTRANAtmosphereInterface* atmo,
                                 const VectorDim1<LineOfSight>& los) = 0;

        virtual void fill_tables(const std::vector<double>& wavelength_nm, const GEODETIC_INSTANT& refpt, SASKTRANAtmosphereInterface* atmo,
                                 const VectorDim1<LineOfSight>& los) = 0;

        virtual void fill_derivative_tables(const std::vector<std::unique_ptr<SKTRAN_DO_UserSpec::WeightingFunctionSpec>>* pert_specs) = 0;
    };

    // Storage for a single wavelength that is common for each NSTOKES
	struct OpticalStateWavelengthEntryCommon {
		std::vector<Eigen::VectorXd> speciesxs;
		std::vector<Eigen::VectorXd> speciesscatxs;
		std::vector<Eigen::MatrixXd> speciesphase;

        std::vector<Eigen::MatrixXd> speciesa1;
        Eigen::MatrixXd totala1;

		Eigen::VectorXd totalext;
		Eigen::VectorXd totalscatext;
		Eigen::MatrixXd totalphase;

		Eigen::VectorXd cumulativeod;

        void resize(size_t num_species, size_t num_altitudes, size_t numstr, size_t numlos) {
            speciesxs.resize(num_species);
            speciesscatxs.resize(num_species);
            speciesa1.resize(num_species);
            speciesphase.resize(num_species);

            for( int i = 0; i < num_species; ++i) {
                speciesxs[i].resize(num_altitudes);
                speciesscatxs[i].resize(num_altitudes);
                speciesa1[i].resize(numstr + 1, num_altitudes); // Add 1 for delta-m scaling later
                speciesa1[i].setZero();

                speciesphase[i].resize(num_altitudes, numlos);
            }

            totalext.resize(num_altitudes);
            totalscatext.resize(num_altitudes);
            totala1.resize(numstr, num_altitudes);  // We delta scale the individual species optical properties then sum for the total
            cumulativeod.resize(num_altitudes);
            totalphase.resize(num_altitudes, numlos);

            totala1.setZero();
            totalext.setZero();
            totalscatext.setZero();
            totalphase.setZero();
        };

	};

    // Additional storage when NSTOKES>1, this could be reduced slightly for NSTOKES=3 but it isn't really
    // a bottleneck
	template <int NSTOKES, int CNSTR=-1>
	struct OpticalStateWavelengthEntry : OpticalStateWavelengthEntryCommon {
        std::vector<Eigen::MatrixXd> speciesa2;
        Eigen::MatrixXd totala2;

        std::vector<Eigen::MatrixXd> speciesa3;
        Eigen::MatrixXd totala3;

        std::vector<Eigen::MatrixXd> speciesb1;
        Eigen::MatrixXd totalb1;

        std::vector<Eigen::MatrixXd> speciesa4;
        Eigen::MatrixXd totala4;

        std::vector<Eigen::MatrixXd> speciesb2;
        Eigen::MatrixXd totalb2;

        void copy_height_independent_species(size_t speciesidx);
        void sum_totals(const nx2dArray<double>& speciesnd, std::vector<SKTRAN_AtmosphericOpticalStateEntry_V21>& allspecies);

        void resize(size_t num_species, size_t num_altitudes, size_t numstr, size_t numlos) {
            OpticalStateWavelengthEntryCommon::resize(num_species, num_altitudes, numstr, numlos);

            speciesa2.resize(num_species);
            speciesa3.resize(num_species);
            speciesa4.resize(num_species);
            speciesb1.resize(num_species);
            speciesb2.resize(num_species);


            for( int i = 0; i < num_species; ++i) {
                speciesa2[i].resize(numstr + 1, num_altitudes); // Add 1 for delta-m scaling later
                speciesa3[i].resize(numstr + 1, num_altitudes); // Add 1 for delta-m scaling later
                speciesa4[i].resize(numstr + 1, num_altitudes); // Add 1 for delta-m scaling later
                speciesb1[i].resize(numstr + 1, num_altitudes); // Add 1 for delta-m scaling later
                speciesb2[i].resize(numstr + 1, num_altitudes); // Add 1 for delta-m scaling later

                speciesa2[i].setZero();
                speciesa3[i].setZero();
                speciesa4[i].setZero();
                speciesb1[i].setZero();
                speciesb2[i].setZero();
            }
            totala2.resize(numstr, num_altitudes);  // We delta scale the individual species optical properties then sum for the total
            totala3.resize(numstr, num_altitudes);  // We delta scale the individual species optical properties then sum for the total
            totala4.resize(numstr, num_altitudes);  // We delta scale the individual species optical properties then sum for the total
            totalb1.resize(numstr, num_altitudes);  // We delta scale the individual species optical properties then sum for the total
            totalb2.resize(numstr, num_altitudes);  // We delta scale the individual species optical properties then sum for the total

            totala2.setZero();
            totala3.setZero();
            totala4.setZero();
            totalb1.setZero();
            totalb2.setZero();
        };
    };

    // Specialized instance for NSTOKES=1, we only need the common storage
	template <>
	struct OpticalStateWavelengthEntry<1> : OpticalStateWavelengthEntryCommon {
        void copy_height_independent_species(size_t speciesidx);

        void sum_totals(const nx2dArray<double> &speciesnd,
                        std::vector<SKTRAN_AtmosphericOpticalStateEntry_V21> &allspecies);

        void resize(size_t num_species, size_t num_altitudes, size_t numstr, size_t numlos) {
            OpticalStateWavelengthEntryCommon::resize(num_species, num_altitudes, numstr, numlos);
        };
    };


    template <int NSTOKES, int CNSTR=-1>
	class OpticalState : public OpticalStateInterface {
    // OpticalState stores the input optical property information on a high resolution altitude grid, this is then
    // integrated over the layers in the DO engine, or used directly in the case of spherical mode
    // The optical state can either be filled on a wavelength by wavelength basis, or it can be prefilled with
    // an array of wavelengths

	public:
        // Define a type for each species. PurelyAbsorbing has no scattering componont, ScatteringHeightDependent
        // is a generic scatter, ScatteringHeightIndependent is a scatterer where the legendre coefficients/cross sections
        // are not a function of altitude.
        enum SpeciesType {
            PurelyAbsorbing,
            ScatteringHeightDependent,
            ScatteringHeightIndependent
        };

		// Mapping between derivatives and values on the altitude grid
		struct DerivEntry {
			int deriv_index;
			double d_totalext;
			double d_totalscatext;
			Eigen::VectorXd d_totalphase;
		};

		OpticalState() {
            m_num_deriv = 0;
		}

		~OpticalState() {

		}

		void configure(const Eigen::VectorXd& internal_altitudes, uint numstr) override {
			m_internal_altitudes = internal_altitudes;
			m_numstr = numstr;
		}

        // Fill method for a single wavelength
        void fill_tables(double wavelength_nm, const GEODETIC_INSTANT& refpt, SASKTRANAtmosphereInterface* atmo,
                         const VectorDim1<LineOfSight>& los) override {
            m_wavelengths.resize(1);
            m_wavelengths[0] = wavelength_nm;
            resize(atmo, los, 1);

            internal_fill_tables(refpt, atmo, los);
        }

        // Fill method for multiple wavelengths
        void fill_tables(const std::vector<double>& wavelength_nm, const GEODETIC_INSTANT& refpt, SASKTRANAtmosphereInterface* atmo,
                         const VectorDim1<LineOfSight>& los) override {
            m_wavelengths = wavelength_nm;
            resize(atmo, los, m_wavelengths.size());

            internal_fill_tables(refpt, atmo, los);
        }

		void internal_fill_tables(const GEODETIC_INSTANT& refpt, SASKTRANAtmosphereInterface* atmo,
			const VectorDim1<LineOfSight>& los) {

			m_brdf = atmo->albedo();
			auto* atmospheric_state = atmo->atmospheric_state();
			auto& all_species = atmo->species();

			GEODETIC_INSTANT loc(refpt.latitude, refpt.longitude, refpt.heightm, refpt.mjd);
			atmospheric_state->UpdateCache(loc);
			bool haschanged;

            // Storage for each thread, Probably not even used anymore
			std::vector<skRTPhaseMatrix> thread_temp_phase(omp_get_num_threads());
			std::vector<int> thread_numlegendre(omp_get_num_threads());
			std::vector<std::vector<double>> thread_legendre(omp_get_num_threads());

            // store m_numstr+1 legendre coefficients for each species to account for delta-m scaling
            for(auto& legendre : thread_legendre) {
                legendre.resize(m_numstr + 1);
            }

            // Temporary wavenumber grid
            std::vector<double> wavenumber(m_wavelengths.size());
            for(size_t idx = 0; idx < m_wavelengths.size(); ++idx) {
                wavenumber[idx] = 1e7 / m_wavelengths[idx];
            }

            // Temporary storage for the output cross sections
            std::vector<double> absxs, scatxs, extxs;
            absxs.resize(m_wavelengths.size());
            scatxs.resize(m_wavelengths.size());
            extxs.resize(m_wavelengths.size());

			// First fill the species specific tables
			for (size_t i = 0; i < all_species.size(); ++i) {
				auto& species = all_species[i];
				auto* clim = species.GetClimatology();
				auto* optprop = species.ParticleOpticalProps();

				optprop->SetAtmosphericState(atmospheric_state);

                // Set the species type
                if(optprop->IsScatterer()) {
                    if(optprop->IsHeightDependent()) {
                        m_species_types[i] = ScatteringHeightDependent;
                    } else {
                        m_species_types[i] = ScatteringHeightIndependent;
                    }
                } else {
                    m_species_types[i] = PurelyAbsorbing;
                }

				for (int j = 0; j < m_internal_altitudes.size(); ++j) {
					loc.heightm = m_internal_altitudes(j);
					clim->UpdateCache(loc);

					clim->GetParameter(species.GetSpecies(), loc, &m_speciesnd.At(i, j), false);

					if (optprop->IsHeightDependent() || j == 0) {
						// If the optical property is height dependent we have to fill the table at every height
						optprop->SetLocation(loc, &haschanged);
                        optprop->CalculateCrossSectionsArray(&wavenumber[0], wavenumber.size(), &absxs[0], &extxs[0], &scatxs[0]);

                        // I think this was causing some problems, CalculatePhaseMatrix is not threadsafe
                        //#pragma omp parallel for schedule(dynamic, 1)
                        for(int wavidx = 0; wavidx < m_wavelengths.size(); ++wavidx) {
                            size_t thread_id = omp_get_thread_num();
                            auto& legendre = thread_legendre[thread_id];
                            auto& numlegendre = thread_numlegendre[thread_id];
                            auto& temp_phase = thread_temp_phase[thread_id];

                            auto& entry = m_wavelengthentry[wavidx];
                            entry.speciesxs[i](j) = extxs[wavidx];
                            entry.speciesscatxs[i](j) = scatxs[wavidx];

                            if (optprop->IsScatterer())
                            {
                                fill_species_scattering_entry(i, wavidx, j, optprop);

                                for (size_t k = 0; k < los.size(); k++) {
                                    optprop->CalculatePhaseMatrix(wavenumber[wavidx], los[k].cos_scattering_angle, &temp_phase);
                                    entry.speciesphase[i](j, k) = temp_phase.At(1, 1);
                                }
                            }

                        }
					}
					else {
                        #pragma omp parallel for schedule(guided, 1)
                        for(int wavidx = 0; wavidx < m_wavelengths.size(); ++wavidx) {
                            auto& entry = m_wavelengthentry[wavidx];
                            entry.copy_height_independent_species(i);
                        }
					}
				}
			}

			// Then combine them together
			// TODO: We need to disable delta scaling for the extinctions when calculating the total average phase function
            // look into this some more
            #pragma omp parallel for schedule(guided, 1)
            for(int wavidx = 0; wavidx < m_wavelengths.size(); ++wavidx) {
                auto& entry = m_wavelengthentry[wavidx];

                entry.sum_totals(m_speciesnd, all_species);
            }
			calculate_optical_depth();
			fill_pressure(refpt, atmo);
		}

        void fill_species_scattering_entry(size_t speciesidx, sktran_do_detail::uint wavelidx, sktran_do_detail::uint heightindex, skOpticalProperties* optprop);

        // Derivative tables are used only in spherical mode
		void fill_derivative_tables(const std::vector<std::unique_ptr<SKTRAN_DO_UserSpec::WeightingFunctionSpec>>* pert_specs) override {
            m_deriv_mapping.resize(m_wavelengths.size()); // I think this size is used for something even if nderiv=0
            for(auto& mapping : m_deriv_mapping) {
                mapping.resize(m_internal_altitudes.size());
            }
            if(pert_specs == nullptr) {
                // don't need to full tables
                return;
            }

			m_num_deriv = pert_specs->size();
			size_t num_los = m_wavelengthentry[0].totalphase.cols();

			auto triangle_function = [](double center, double width, double x) {
				double dist = abs(center - x);
				return 1 - dist / width;
			};

            // Here we calculate the mapping of species perturbations to elements in the optical tables
			for (int i = 0; i < pert_specs->size(); ++i) {
				if (auto spec = dynamic_cast<const SKTRAN_DO_UserSpec::SpeciesWF*>(pert_specs->at(i).get())) {
					size_t species_index = -1;
					for (size_t j = 0; j < m_clim_handles.size(); ++j) {
						if (m_clim_handles[j] == spec->handle) {
							species_index = j;
						}
					}

					double low_alt = spec->altitude - spec->l_width;
					double high_alt = spec->altitude + spec->u_width;

					// Find the altitude grid elements between these two values
					auto lower_bound = std::lower_bound(m_internal_altitudes.cbegin(), m_internal_altitudes.cend(), low_alt);
					auto upper_bound = std::lower_bound(m_internal_altitudes.cbegin(), m_internal_altitudes.cend(), high_alt);

					for (auto it = lower_bound; it != upper_bound; it++) {
						size_t index = std::distance(m_internal_altitudes.cbegin(), it);
						double triangle_val;
						if (*it >= spec->altitude) {
							triangle_val = triangle_function(spec->altitude, spec->u_width, *it);
						}
						else {
							triangle_val = triangle_function(spec->altitude, spec->l_width, *it);
						}
                        #pragma omp parallel for schedule(guided, 1)
                        for(int wavidx = 0; wavidx < m_wavelengths.size(); ++wavidx) {
                            DerivEntry entry;
                            const auto& tableentry = m_wavelengthentry[wavidx];
                            // ext = nd * xs * 100
                            // d_ext = xs * 100
                            entry.d_totalext = tableentry.speciesxs[species_index](index) * triangle_val * 100;
                            entry.d_totalscatext = tableentry.speciesscatxs[species_index](index) * triangle_val * 100;
                            entry.d_totalphase.resize(num_los);
                            if(m_species_types[species_index] != PurelyAbsorbing)
                            {
                                for (uint k = 0; k < num_los; ++k) {
                                    entry.d_totalphase(k) = tableentry.speciesscatxs[species_index](index) * triangle_val * 100 / (tableentry.totalscatext(index)) * (tableentry.speciesphase[species_index](index, k) - tableentry.totalphase(index, k));
                                }
                            } else {
                                for (uint k = 0; k < num_los; ++k) {
                                    entry.d_totalphase(k) = 0;
                                }
                            }

                            entry.deriv_index = i;

                            m_deriv_mapping[wavidx][index].emplace_back(entry);
                        }
					}
				}
			}

		}

		void calculate_optical_depth() {
            // Calculates the vertical OD of the atmospheric state
			// To calculate the vertical OD we just do standard trapezoidal integration along the path
            // This should be exact assuming linear variation of optical properties in altitude
            #pragma omp parallel for schedule(guided, 1)
            for(int wavidx = 0; wavidx < m_wavelengths.size(); ++wavidx) {
                auto& entry = m_wavelengthentry[wavidx];
                double od = 0.0;

                for (auto i = m_internal_altitudes.size() - 1; i > 0; --i) {
                    entry.cumulativeod(i) = od;

                    double width = m_internal_altitudes(i) - m_internal_altitudes(i - 1);
                    double height = (entry.totalext(i) + entry.totalext(i - 1)) / 2.0;

                    od += width * height;
                }

                entry.cumulativeod(0) = od;

                entry.cumulativeod.reverseInPlace();
            }
		}

		void fill_pressure(const GEODETIC_INSTANT& refpt, SASKTRANAtmosphereInterface* atmo) {
            // Fill the pressure tables, needed to later place layers in pressure space
			GEODETIC_INSTANT loc(refpt.latitude, refpt.longitude, refpt.heightm, refpt.mjd);

			atmo->atmospheric_state()->UpdateCache(loc);

			for (int i = 0; i < m_internal_altitudes.size(); ++i) {
				loc.heightm = m_internal_altitudes[i];

				atmo->atmospheric_state()->GetParameter(SKCLIMATOLOGY_PRESSURE_PA, loc, &m_pressure(i), false);
			}
		}

		const Eigen::VectorXd& cumulative_od(size_t wavel_index) const {
			return m_wavelengthentry[wavel_index].cumulativeod;
		}

        void assign_legendre_quantities(Eigen::Matrix<double, Eigen::Dynamic, 6>& bl, const OpticalStateWavelengthEntry<NSTOKES>& entry, double altitude_weight, double cs_scat, const std::pair<size_t, size_t>& indicies, const std::pair<double, double>& weights ) const;
        void assign_species_legendre_quantities(Eigen::Matrix<double, Eigen::Dynamic, 6>& bl, const OpticalStateWavelengthEntry<NSTOKES>& entry, int species_idx,  const std::pair<size_t, size_t>& indicies, const std::pair<double, double>& weights, SpeciesType type ) const;


        void compute_layer_quantities_from_interpolation(double& k, double& kscat, Eigen::Matrix<double, Eigen::Dynamic, 6>& bl, const Eigen::MatrixXd& interpolator, LayerIndex p, size_t wavel_index) const {
            const auto& entry = m_wavelengthentry[wavel_index];

            k = interpolator(p, Eigen::all).dot(entry.totalext);
            kscat = interpolator(p, Eigen::all).dot(entry.totalscatext);

            for(int i = 0; i < m_numstr; ++i) {
                bl(i, 0) = interpolator(p, Eigen::all).dot(entry.totalscatext.cwiseProduct(entry.totala1(i, Eigen::all).transpose()));

                if constexpr(NSTOKES > 1) {
                    bl(i, 1) = interpolator(p, Eigen::all).dot( entry.totalscatext.cwiseProduct(entry.totala2(i, Eigen::all).transpose()) );
                    bl(i, 2) = interpolator(p, Eigen::all).dot( entry.totalscatext.cwiseProduct(entry.totala3(i, Eigen::all).transpose()) );
                    bl(i, 3) = interpolator(p, Eigen::all).dot( entry.totalscatext.cwiseProduct(entry.totala4(i, Eigen::all).transpose()) );
                    bl(i, 4) = interpolator(p, Eigen::all).dot( entry.totalscatext.cwiseProduct(entry.totalb1(i, Eigen::all).transpose()) );
                    bl(i, 5) = interpolator(p, Eigen::all).dot( entry.totalscatext.cwiseProduct(entry.totalb2(i, Eigen::all).transpose()) );
                }
            }
            bl /= kscat;
        }


        void add_discretized_trapezoid(double upper_altitude, double lower_altitude,  double& cs_scat, double& cs_ext, Eigen::Matrix<double, Eigen::Dynamic, 6>& bl, size_t wavel_index) const {
            // We have altitudes lower_altitude and upper_altitude that are guaranteed not to cross opticalstate boundaries,
            // we calculate the integral of the optical state quantities over this altitude range and accumulate them

            // Used in constructing the layer quantities

            const auto& entry = m_wavelengthentry[wavel_index];

            double cs_scat_lower, cs_scat_upper;

            std::pair<size_t, size_t> indicies_lower, indicies_upper;
            std::pair<double, double> weights_lower, weights_upper;

            double altitude_weight = (upper_altitude - lower_altitude) * 0.5;

            interpolation_index(upper_altitude, indicies_upper, weights_upper);
            interpolation_index(lower_altitude, indicies_lower, weights_lower);

            // Scattering extinctions
            cs_scat_lower = entry.totalscatext(indicies_lower.first) * weights_lower.first + entry.totalscatext(indicies_lower.second) * weights_lower.second;
            cs_scat_upper = entry.totalscatext(indicies_upper.first) * weights_upper.first + entry.totalscatext(indicies_upper.second) * weights_upper.second;

            cs_scat += altitude_weight * (cs_scat_lower + cs_scat_upper);

            // Total extinctions
            cs_ext += altitude_weight * ( entry.totalext(indicies_lower.first) * weights_lower.first + entry.totalext(indicies_lower.second) * weights_lower.second );
            cs_ext += altitude_weight * ( entry.totalext(indicies_upper.first) * weights_upper.first + entry.totalext(indicies_upper.second) * weights_upper.second );

            // Legendre quantities
            assign_legendre_quantities(bl, entry, altitude_weight, cs_scat_lower, indicies_lower, weights_lower);
            assign_legendre_quantities(bl, entry, altitude_weight, cs_scat_upper, indicies_upper, weights_upper);
        }

        void fill_interpolation_matrix(double upper_altitude, double lower_altitude,  LayerIndex p, Eigen::MatrixXd& interpolation) const {
            // We have altitudes lower_altitude and upper_altitude that are guaranteed not to cross opticalstate boundaries,
            // we calculate the integral of the optical state quantities over this altitude range and accumulate them

            std::pair<size_t, size_t> indicies_lower, indicies_upper;
            std::pair<double, double> weights_lower, weights_upper;

            double altitude_weight = (upper_altitude - lower_altitude) * 0.5;

            interpolation_index(upper_altitude, indicies_upper, weights_upper);
            interpolation_index(lower_altitude, indicies_lower, weights_lower);

            interpolation(p, indicies_lower.first) += weights_lower.first * altitude_weight;
            interpolation(p, indicies_lower.second) += weights_lower.second * altitude_weight;

            interpolation(p, indicies_upper.first) += weights_upper.first * altitude_weight;
            interpolation(p, indicies_upper.second) += weights_upper.second * altitude_weight;
        }


		SpeciesType interpolate_species(const CLIMATOLOGY_HANDLE& handle, double altitude, double& scatxs, double& totalxs, double& numden, Eigen::Matrix<double, Eigen::Dynamic, 6>& legendre, size_t wavel_index) const {
            // Interpolates a specific species, used to calculate the derivative quantities

            const auto& entry = m_wavelengthentry[wavel_index];

			size_t species_index = -1;
			for (size_t i = 0; i < m_clim_handles.size(); ++i) {
				if (m_clim_handles[i] == handle) {
					species_index = i;
				}
			}
			if (species_index == -1) {

			}
			std::pair<size_t, size_t> indicies;
			std::pair<double, double> weights;

			interpolation_index(altitude, indicies, weights);

			scatxs = entry.speciesscatxs[species_index](indicies.first) * weights.first + entry.speciesscatxs[species_index](indicies.second) * weights.second;
			totalxs = entry.speciesxs[species_index](indicies.first) * weights.first + entry.speciesxs[species_index](indicies.second) * weights.second;
			numden = m_speciesnd.At(species_index, indicies.first) * weights.first + m_speciesnd.At(species_index, indicies.second) * weights.second;

            assign_species_legendre_quantities(legendre, entry, species_index, indicies, weights, m_species_types[species_index]);

            return m_species_types[species_index];
		}

		double altitude_at_pressure(double pressure) const {
			auto lte = std::upper_bound(m_pressure.cbegin(), m_pressure.cend(), pressure, std::greater<double>());
			if (lte == m_pressure.cend()) {
				return m_internal_altitudes(Eigen::last);
			}
			double prop = -1 * (std::log(*(lte - 1)) - std::log(pressure)) / (std::log(*lte) - std::log(*(lte - 1)));
			auto idx = std::distance(m_pressure.cbegin(), lte) - 1;
			double above_alt = m_internal_altitudes(idx + 1);
			double below_alt = m_internal_altitudes(idx);

			return below_alt + prop * (above_alt - below_alt);
		}

		double opticaldepth_at_altitude(double altitude, size_t wavel_index) const {
            const auto& entry = m_wavelengthentry[wavel_index];

			auto lte = std::upper_bound(m_internal_altitudes.cbegin(),
				m_internal_altitudes.cend(),
				altitude
			);
			if (lte == m_internal_altitudes.cend()) {
				return entry.cumulativeod(0);
			}
			// indexes to cumlative_od are idx, idx+1
			auto idx = std::distance(lte, m_internal_altitudes.cend()) - 1;

			// indexes to heights are idx_h and idx_h-1
			auto idx_h = std::distance(m_internal_altitudes.cbegin(), lte);

			double delta_h = m_internal_altitudes[idx_h] - m_internal_altitudes[idx_h - 1];
			double delta_k = entry.totalext[idx_h - 1] - entry.totalext[idx_h];

			double h = m_internal_altitudes[idx_h] - altitude;
			double k0 = entry.totalext[idx_h];

			return entry.cumulativeod[idx] + h * k0 + (delta_k / delta_h) * (h*h) / 2;
		}

		double altitude_at_opticaldepth(double optical_depth, size_t wavel_index) const {
            const auto& entry = m_wavelengthentry[wavel_index];

			auto lte = std::upper_bound(entry.cumulativeod.cbegin(), entry.cumulativeod.cend(), optical_depth);
			if (lte == entry.cumulativeod.cend()) {
				return m_internal_altitudes(0);
			}
			double od_start = *(lte - 1);

			// idx in height
			auto idx_h = std::distance(lte, entry.cumulativeod.cend()) - 1;

			double delta_h = m_internal_altitudes[idx_h + 1] - m_internal_altitudes[idx_h];
			double delta_k = entry.totalext[idx_h] - entry.totalext[idx_h + 1];
			double k0 = entry.totalext[idx_h + 1];

			if (std::abs(od_start - optical_depth) < 1e-10) {
				return m_internal_altitudes[idx_h + 1];
			}
			double a = delta_k / delta_h / 2;
			double b = k0;
			double c = od_start - optical_depth;
			double h;

			h = 2.0 * c / (-b - sqrt(b*b - 4 * a*c));
			return m_internal_altitudes[idx_h + 1] - h;
		}

		const Eigen::VectorXd& pressure() const {
			return m_pressure;
		}

		const Eigen::VectorXd& total_extinction(size_t wavel_index) const {
			return m_wavelengthentry[wavel_index].totalext;
		}

		const Eigen::VectorXd& scattering_extinction(size_t wavel_index) const {
			return m_wavelengthentry[wavel_index].totalscatext;
		}

		const Eigen::MatrixXd& phase_function(size_t wavel_index) const {
			return m_wavelengthentry[wavel_index].totalphase;
		}

		double exact_brdf(const HELIODETIC_UNITVECTOR& incoming_ground, const HELIODETIC_POINT& ground_pt, const SKTRAN_CoordinateTransform_V2& coords, size_t wavel_index) const {
			GEODETIC_INSTANT pt;

			// TODO: Do this calculation and linearize
			double mu_in, mu_out, cosdphi, brdf;
			mu_in = 0.0;
			mu_out = 0.0;
			cosdphi = 0.0;

			m_brdf->BRDF(m_wavelengths[wavel_index], pt, mu_in, mu_out, cosdphi, &brdf);

			return brdf;
		}

		const skBRDF* brdf_object() const {
			return m_brdf;
		}

        size_t wavel_index(double wavelength) const {
            size_t current_index = 0;
            double current_difference = abs(m_wavelengths[current_index] - wavelength);

            for(size_t idx = 1; idx < m_wavelengths.size(); ++idx) {
                if(abs(m_wavelengths[idx] - wavelength) < current_difference) {
                    current_difference = abs(m_wavelengths[idx] - wavelength);
                    current_index = idx;
                }
            }
            return current_index;
        }


	private:
		void resize(SASKTRANAtmosphereInterface* atmo, const VectorDim1<LineOfSight>& los, size_t num_wavel = 1) {
			size_t num_species = atmo->species().size();
			m_clim_handles.resize(num_species);
            m_species_types.resize(num_species);
			for (size_t i = 0; i < num_species; ++i) {
				m_clim_handles[i] = atmo->species()[i].GetSpecies();
			}
            // Wavelength independent quantities
            m_speciesnd.SetSize(num_species, m_internal_altitudes.size());
            m_pressure.resize(m_internal_altitudes.size());
            m_pressure.setZero();

            m_wavelengthentry.resize(num_wavel);

            // Wavelength dependent quantities
            for( auto& entry : m_wavelengthentry) {
                entry.resize(num_species, m_internal_altitudes.size(), m_numstr, los.size());
            }
		}

		void interpolation_index(double altitude, std::pair<size_t, size_t>& indicies, std::pair<double, double>& weights) const {
			nxLinearInterpolate::FindBoundingIndicesAscending(m_internal_altitudes.cbegin(), m_internal_altitudes.cend(), altitude, &indicies.first,
				&indicies.second, &weights.first, &weights.second);

			weights.first = (weights.second - altitude) / (weights.second - weights.first);
			weights.second = (1 - weights.first);
		}

		const InputDerivatives<NSTOKES>* m_input_derivatives;
		const skBRDF* m_brdf;

        nx2dArray<double> m_speciesnd;

        std::vector<double> m_wavelengths;
        std::vector<OpticalStateWavelengthEntry<NSTOKES>> m_wavelengthentry;

        Eigen::VectorXd m_pressure;

		uint m_numstr;
		Eigen::VectorXd m_internal_altitudes;
		std::vector<CLIMATOLOGY_HANDLE> m_clim_handles;
        std::vector<SpeciesType> m_species_types;

		VectorDim3<DerivEntry> m_deriv_mapping;  // [wavel, altitude, deriv]
		size_t m_num_deriv;

		public:
			const VectorDim2<DerivEntry>& deriv_mapping(size_t wavel_index) const {
				return m_deriv_mapping[wavel_index];
			}

			size_t num_deriv() const {
				return m_num_deriv;
			}
	};
}