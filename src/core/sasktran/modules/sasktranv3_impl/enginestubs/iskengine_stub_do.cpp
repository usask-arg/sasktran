#include "../dllimplementation/stdafx.h"
#include "modules/sasktranv3_impl/sktranif_impl_helperclasses.h"

#include <string>
#include <boost/regex.hpp>

// Sasktran-Disco ISK interface

ISKEngine_Stub_DO::ISKEngine_Stub_DO()
{
    m_radiance.SetReuseMemory(true);
    makeScalarSetFunctions();
    makeVectorSetFunctions();
    makeScalarGetFunctions();
    makeVectorGetFunctions();
    makeObjectSetFunctions();
    makeStringSetFunctions();

    // Set default values
    m_userspec.setNumberOfStreams(16);
    m_userspec.setNumberOfLayers(50);
    m_return_los_diagnostics = false;
    m_return_rts_diagnostics = false;
    m_force_single_reference_point = false;

    m_output_optical_depths = false;

    m_num_threads = 0;

    m_wl_batches = 1;
    m_num_stokes = 1;

    m_registered_calls.resize((unsigned int)TrackedCalls::SIZE, false);

    m_required_calls.resize((unsigned int)TrackedCalls::SIZE, true);
    m_required_calls[(unsigned int)TrackedCalls::AtmosphericState] = false;
    m_required_calls[(unsigned int)TrackedCalls::WFAltitudes] = false;
    m_required_calls[(unsigned int)TrackedCalls::SunPosition] = false;

}

bool ISKEngine_Stub_DO::SetBRDF(ISKBrdf_Stub* brdf)
{
    skBRDF* brdf_cast = static_cast<skBRDF*>(brdf->RawObjectPointer());

    if(!brdf_cast) {
        nxLog::Record(NXLOG_ERROR, "Invalid skBRDF object!");
        return false;
    } else {
        m_registered_calls[(unsigned int)TrackedCalls::BRDF] = true;
        m_opticalstate.set_albedo(brdf_cast);
        return true;
    }
}

ISKEngine_Stub_DO::~ISKEngine_Stub_DO()
{
}

bool ISKEngine_Stub_DO::AddLineOfSight(double mjd, const nxVector& obs, const nxVector& look, int* losindex)
{
    bool		ok;

    ok = m_linesofsight.AddLineOfSight(obs, look, mjd);
    if(ok) {
        *losindex = (int) m_linesofsight.NumRays() - 1;
    } else {
        *losindex = -999999;
    }
    m_radiance.erase();
    m_registered_calls[(unsigned int)TrackedCalls::LinesOfSight] = true;
    return ok;
}

bool ISKEngine_Stub_DO::SetPolarizationMode(int specifier)
{
    m_num_stokes = specifier;
    return true;
}

bool ISKEngine_Stub_DO::AddSpecies(const CLIMATOLOGY_HANDLE& species, ISKClimatology_Stub* climatology, ISKOpticalProperty_Stub* opticalproperty)
{
    skClimatology*				climptr;
    nxUnknown*					optbaseptr;
    skOpticalProperties*		optptr;

    optbaseptr = (opticalproperty != nullptr) ? opticalproperty->RawObjectPointer() : nullptr;
    climptr = static_cast<skClimatology*>(climatology->RawObjectPointer());
    optptr = (optbaseptr != nullptr) ? static_cast<skOpticalProperties*>(optbaseptr) : nullptr;
    m_opticalstate.add_species(species, climptr, optptr);
    m_registered_calls[(unsigned int)TrackedCalls::AtmosphericSpecies] = true;
    return true;
}

bool ISKEngine_Stub_DO::AddEmission(const EMISSION_HANDLE& species, ISKEmission_Stub* emissionobject)
{
    nxLog::Record(NXLOG_ERROR, "Sasktran-Disco does not support emissions.");
    return false;
}

bool ISKEngine_Stub_DO::SetAlbedo(double albedo)
{
    m_registered_calls[(unsigned int)TrackedCalls::BRDF] = true;
    m_opticalstate.set_albedo(albedo);
    return true;
}

bool ISKEngine_Stub_DO::SetWavelengths(const double* wavelen, int numwavelen)
{
    m_registered_calls[(unsigned int)TrackedCalls::Wavelengths] = true;
    m_wavelen.assign(wavelen, wavelen + numwavelen);
    return true;
}

bool ISKEngine_Stub_DO::InitializeModel()
{
    try {
        m_userspec.configure();
        m_registered_calls[(unsigned int)TrackedCalls::Initialize] = true;
        return true;
    }
    catch(const std::exception& e) {
        nxLog::Record(NXLOG_ERROR, e.what());
        return false;
    }
}

bool ISKEngine_Stub_DO::CalculateRadiance(const double** radiance, int* numwavelens, int* numlinesofsight)
{
    // The matrices used here are so small that mkl threading hurts more than it helps
#ifdef SKTRAN_DO_USE_MKL
    mkl_set_num_threads(1);
#else
#ifndef SKTRAN_DO_USE_ACCELERATE
    //    LAPACKE_set_nancheck(0);
#endif
#endif
    omp_set_dynamic(0);
    omp_set_num_threads(m_num_threads);

    m_radiance.erase();
    auto numrays = m_linesofsight.NumRays();
    auto numwave = m_wavelen.size();
    m_radiance.SetSize(numrays * m_num_stokes, numwave);
    m_lostransmissionstorage.SetSize(numrays, numwave);
    m_lostransmissionstorage.SetTo(0.0);

    auto handle_calc_except = [&, this](const std::exception& e, int w, int l=-1) {
        // Try to cast the exception to an Disco exception
        bool found_type = false;
        auto internal_except = dynamic_cast<const sktran_do_detail::InternalError*>(&e);
        if(internal_except != nullptr) {
            nxLog::Record(NXLOG_ERROR, internal_except->what());
            found_type = true;
        }
        auto config_except = dynamic_cast<const sktran_do_detail::InvalidConfiguration*>(&e);
        if(config_except != nullptr) {
            nxLog::Record(NXLOG_ERROR, config_except->what());
            found_type = true;
        }
        if(!found_type) {
            nxLog::Record(NXLOG_ERROR, "An unexpected internal exception occured. This is likely a bug! "
                                       "Please submit a issue at: https://arggit.usask.ca/ARGPackages/SasktranDO. ");
        }

        // Set associated m_radiance's to NaN
        if(l != -1) {
            m_radiance.At(l, w) = std::numeric_limits<double>::quiet_NaN();
        } else {
            for(int i = 0; i < numrays; ++i) {
                m_radiance.At(i, w) = std::numeric_limits<double>::quiet_NaN();
            }
        }
    };

    try {
        // check that we are ready
        bool required_calls_not_okay = false;
        for(size_t i = 0; i < m_required_calls.size(); ++i) {
            if(m_required_calls[i] && !m_registered_calls[i]) {
                required_calls_not_okay = true;
            }
        }
        if(required_calls_not_okay) {
            throw sktran_do_detail::InvalidConfiguration(
                    "All required calls have not been called"
            );
        }

        // configure userspecs
        m_userspec.configure();

        if( m_output_optical_depths ) {
            m_ssastorage.SetSize(m_userspec.getNumberOfLayers() , numwave);
            m_odstorage.SetSize(m_userspec.getNumberOfLayers(), numwave);
        }

        // Configure weighting calculations
        // If the user has set the wfaltitudes, good, if not we default to the altitudegrid
        if (m_wf_altitudes.size() == 0) {
            m_wf_altitudes = m_userspec.getAltitudeGrid();
        }

        // If the user has set the wf widths, use them otherwise construct them from the altitudes
        if (m_wf_widths_low.size() == 0) {
            m_wf_widths_low.resize(m_wf_altitudes.size(), 0.0);
            m_wf_widths_high.resize(m_wf_altitudes.size(), 0.0);

            for (int i = 0; i < m_wf_altitudes.size() - 1; ++i) {
                m_wf_widths_high[i] = m_wf_altitudes[i + 1] - m_wf_altitudes[i];
                m_wf_widths_low[i + 1] = m_wf_altitudes[i + 1] - m_wf_altitudes[i];
            }
            // Duplicate the last value to the edges
            m_wf_widths_high[m_wf_widths_high.size() - 1] = m_wf_widths_high[m_wf_widths_high.size() - 2];
            m_wf_widths_low[0] = m_wf_widths_low[1];
        }

        // check that wf specs are sized correctly
        if (m_wf_widths_low.size() > 0 && m_wf_widths_low.size() != m_wf_altitudes.size()) {
            throw sktran_do_detail::InvalidConfiguration("Perturbation widths must be specified for each weighting function!");
        }
        if (m_wf_widths_high.size() > 0 && m_wf_widths_high.size() != m_wf_altitudes.size()) {
            throw sktran_do_detail::InvalidConfiguration("Perturbation widths must be specified for each weighting function!");
        }

        SKTRAN_DO_UserSpec::VectorOfUPtrWeightingFunctionSpecs wfs;
        m_albedo_wf_included  = false;
        wfs.reserve(m_wf_handles.size() * m_wf_altitudes.size());
        for(auto handle : m_wf_handles) {
            if (handle == SKCLIMATOLOGY_ALBEDO) {
                m_albedo_wf_included = true;
            }
            else {
                for (int i = 0; i < m_wf_altitudes.size(); ++i) {
                    wfs.push_back(std::unique_ptr<SKTRAN_DO_UserSpec::WeightingFunctionSpec>(
                            new SKTRAN_DO_UserSpec::SpeciesWF(handle, m_wf_altitudes.at(i), m_wf_widths_low.at(i), m_wf_widths_high.at(i))
                    ));
                }
            }
        }
        if (m_albedo_wf_included) {
            // add albedo to the END of the wf list
            wfs.push_back(std::unique_ptr<SKTRAN_DO_UserSpec::WeightingFunctionSpec>(
                    new SKTRAN_DO_UserSpec::AlbedoWF()
            ));
        }

        // reset the size of our wf storage vector (note that wfs will be length 0 if no weighting functions)
        m_wf_storage.reset(new double[m_wavelen.size() * m_linesofsight.NumRays() * wfs.size() * m_num_stokes]);
        m_userspec.setWFSpecies(wfs);

        // Get ready to do the calculation
        m_los_diagnostics.clear();
        m_rts_diagnostics.clear();
        if(m_return_los_diagnostics) m_los_diagnostics.resize(numrays);
        if(m_return_rts_diagnostics) m_rts_diagnostics.resize(numwave, std::vector<sktran_do_detail::RTSDiagnostics>(numrays));

        // Start the calculation
        std::mutex optical_properties_mutex;

        if(!m_force_single_reference_point) {
#pragma omp parallel for schedule(dynamic, 1)
            for(int wl = 0; wl < numwave * numrays; ++wl) {
                auto w = wl / numrays;
                auto l = wl % numrays;
                // Make a line-of-sight array with the single line of sight
                auto single_los = SKTRAN_LineOfSightArray_V21();
                const SKTRAN_LineOfSightEntry_V2* the_only_los_entry;
                m_linesofsight.GetRay(l, &the_only_los_entry);
                single_los.AddLineOfSight(the_only_los_entry->Observer(),
                                          the_only_los_entry->Look(),
                                          the_only_los_entry->Mjd());

                // Configure line-of-sight diagnostics storage memory
                std::vector<sktran_do_detail::LOSDiagnostics>* los_diag = nullptr;
                if(w == 0 && m_return_los_diagnostics) {
                    // only needs to be done once, wavelength independent so only do it for w=0
                    los_diag = &m_los_diagnostics[l];
                }

                // Get ready for calculation
                sktran_do_detail::RTSDiagnostics* rts_diag = nullptr;
                if(m_return_rts_diagnostics) {
                    rts_diag = &m_rts_diagnostics[w][l];
                }

                // Do the calculation
                try {
                    std::vector<double> rad;
                    std::vector<double> los_transmission;
                    std::vector<std::vector<double>> los_wf;
                    std::vector<double> od;
                    std::vector<double> ssa;

                    auto* input_od = m_output_optical_depths && l == 0 ? &od : nullptr;
                    auto* input_ssa = m_output_optical_depths && l == 0 ? &ssa : nullptr;

                    // Configure a new engine
                    SKTRAN_DO_EngineInterface* engine;
                    SKTRAN_DO_Engine<1> engine1;
                    SKTRAN_DO_Engine<3> engine3;

                    if(m_num_stokes == 1) {
                        engine1.configureModel(m_userspec, single_los, &optical_properties_mutex, los_diag);
                        engine = &engine1;
                    } else if (m_num_stokes == 3) {
                        engine3.configureModel(m_userspec, single_los, &optical_properties_mutex, los_diag);
                        engine = &engine3;
                    }
                    engine->calculateRadiance(rad, m_wavelen[w], &m_opticalstate, &los_wf, rts_diag, &los_transmission, input_od, input_ssa);


                    int wav_offset = static_cast<int>(m_linesofsight.NumRays() * m_num_stokes * wfs.size());
                    int los_offset = static_cast<int>(wfs.size() * m_num_stokes);
                    std::copy(los_wf[0].cbegin(), los_wf[0].cend(), m_wf_storage.get() + w * wav_offset + l * los_offset);
                    for(int s = 0; s < m_num_stokes; ++s) {
                        m_radiance.At(l*m_num_stokes + s, w) = rad[s];
                    }
                    m_lostransmissionstorage.At(l, w) = los_transmission[0];

                    if(m_output_optical_depths && l == 0 ) {
                        for( int k = 0; k < od.size(); ++k ) {
                            m_odstorage.At(k, w) = od[k];
                            m_ssastorage.At(k, w) = ssa[k];
                        }
                    }

                }
                catch(const std::exception& e) {
                    handle_calc_except(e, (int) w, (int) l);
                }
            }
        } else { // force single reference point -> only loop over wavelength

            omp_set_dynamic(0);
            omp_set_num_threads(m_num_threads);

            int wlbatches = m_wl_batches;
            int numinbatch = ceil(double(numwave) / double(wlbatches));
            for (int batch = 0; batch < wlbatches; batch++) {
                int batchstart = batch * numinbatch;
                int batchend = min(int(numwave), int((batch + 1) * numinbatch));

                std::vector<double> batchwv(batchend - batchstart);
                for (int w = 0; w < batchwv.size(); ++w) {
                    batchwv[w] = m_wavelen[w + batchstart];
                }

                SKTRAN_DO_EngineInterface* engine;
                SKTRAN_DO_Engine<1> engine1;
                SKTRAN_DO_Engine<3> engine3;

                if(m_num_stokes == 1) {
                    engine1.configureModel(m_userspec, m_linesofsight, &optical_properties_mutex, nullptr);
                    engine1.prefillWavelengthTables(batchwv, &m_opticalstate, &optical_properties_mutex);
                    engine = &engine1;
                } else if (m_num_stokes == 3) {
                    engine3.configureModel(m_userspec, m_linesofsight, &optical_properties_mutex, nullptr);
                    engine3.prefillWavelengthTables(batchwv, &m_opticalstate, &optical_properties_mutex);
                    engine = &engine3;
                }

#pragma omp parallel for schedule(guided) num_threads(m_num_threads)
                for (int w = batchstart; w < batchend; ++w) {
                    // Temporary space to hold LOS diagnostics
                    std::vector<sktran_do_detail::LOSDiagnostics> los_diag;
                    decltype(los_diag)* los_diag_ptr = m_return_los_diagnostics ? &los_diag : nullptr;

                    // Get ready for calculation
                    sktran_do_detail::RTSDiagnostics* rts_diag_ptr = nullptr;
                    if (m_return_rts_diagnostics) {
                        // Single refpt means single atmosphere means single R.T. solution (different postprocessing though).
                        m_rts_diagnostics[w].resize(1);
                        rts_diag_ptr = &(m_rts_diagnostics[w][0]);
                    }

                    // Try the calculation
                    try {
                        std::vector<double> rad;
                        std::vector<double> los_transmission;
                        std::vector<std::vector<double>> los_wf;

                        std::vector<double> od;
                        std::vector<double> ssa;

                        auto* input_od = m_output_optical_depths ? &od : nullptr;
                        auto* input_ssa = m_output_optical_depths ? &ssa : nullptr;

                        int linearwavelidx = w - batchstart;

                        engine->calculateRadiance(rad, m_wavelen[w], &m_opticalstate, &los_wf, rts_diag_ptr, &los_transmission, input_od, input_ssa, linearwavelidx);

                        if (m_output_optical_depths) {
                            for (int k = 0; k < od.size(); ++k) {
                                m_odstorage.At(k, w) = od[k];
                                m_ssastorage.At(k, w) = ssa[k];
                            }
                        }

                        int wav_offset = static_cast<int>(m_linesofsight.NumRays() * m_num_stokes * wfs.size());
                        int los_offset = static_cast<int>(wfs.size() * m_num_stokes);
                        for (int l = 0; l < numrays; ++l) {
                            std::copy(los_wf[l].cbegin(), los_wf[l].cend(), m_wf_storage.get() + w * wav_offset + l * los_offset);
                            for( int s = 0; s < m_num_stokes; ++s) {
                                m_radiance.At(l*m_num_stokes + s, w) = rad[l*m_num_stokes + s];
                            }
                            m_lostransmissionstorage.At(l, w) = los_transmission[l];

                            if (m_return_los_diagnostics) {
                                m_los_diagnostics[l].clear();
                                m_los_diagnostics[l].reserve(1);
                                m_los_diagnostics[l].push_back(sktran_do_detail::LOSDiagnostics(los_diag[l]));
                            }
                        }
                    }
                    catch (const std::exception& e) { // Mark NaN if an error occured
                        handle_calc_except(e, w);
                    }
                }
            }
        }
    }
    catch(const std::exception& e) {
        auto except = dynamic_cast<const sktran_do_detail::InternalError*>(&e);
        if(except != nullptr) {
            nxLog::Record(NXLOG_ERROR, except->what());
        } else {
            nxLog::Record(NXLOG_ERROR, "An unexpected internal exception occured. This is likely a bug! "
                                       "Please submit a issue at: https://arggit.usask.ca/ARGPackages/SasktranDO. ");
        }
        *radiance = m_radiance.ArrayBasePtr();
        *numwavelens = (int) numwave;
        *numlinesofsight = (int) numrays;
        return false;
    }

    *radiance = m_radiance.ArrayBasePtr();
    *numwavelens = (int) numwave;
    *numlinesofsight = (int) numrays;
    return true;
}

bool ISKEngine_Stub_DO::CalculateStokesVector(const ISKStokesVector** radiance, int* numwavelens, int* numlinesofsight)
{
    const double *radiancetemp;

    bool ok = CalculateRadiance(&radiancetemp, numwavelens, numlinesofsight);

    m_stokesvector.resize((*numwavelens) * (*numlinesofsight));
    const SKTRAN_LineOfSightEntry_V2* entry;

    for( int l = 0; l < *numlinesofsight; ++l) {

        m_linesofsight.GetRay(l, &entry);

        const nxVector propagation = -1 * entry->Look();
        nxVector geo_theta = entry->Observer().UnitVector().Cross(propagation).UnitVector();
        // There is an ambiguity if the observer unit vector is directly in line with the propagation direction
        if(geo_theta.IsZero() || !geo_theta.IsValid()) {
            // in this case just set theta to be north?
            nxGeodetic geo(entry->Observer().Latitude(), entry->Observer().Longitude(), entry->Observer().Z());
            nxVector temp1, temp2;
            geo.GetGeodeticWestSouthUp(&temp1, &geo_theta, &temp2);
        }


        const nxVector geo_phi = propagation.Cross(geo_theta);

        ISKBasisDirection stokes_basis;
        stokes_basis.Assign(propagation, geo_theta, geo_phi);

        for( int w = 0; w < *numwavelens; ++w) {
            int linearindex = w * (*numlinesofsight) + l;

            IQUV stokes;
            stokes.I = m_radiance.at(l*m_num_stokes, w);
            if( m_num_stokes > 1) {
                stokes.Q = m_radiance.at(l*m_num_stokes + 1, w);
                if(m_num_stokes > 2) {
                    stokes.U = m_radiance.at(l*m_num_stokes + 2, w);
                    if(m_num_stokes > 3) {
                        stokes.V = m_radiance.at(l*m_num_stokes + 3, w);
                    }
                }
            }
            m_stokesvector[linearindex].Assign(stokes, stokes_basis);
        }
    }
    *radiance = m_stokesvector.data();

    return ok;
}

bool ISKEngine_Stub_DO::GetWeightingFunctions(const double** wf, int* numwavel, int* numlos, int* numwf)
{
    if(m_wf_handles.size() == 0) {
        nxLog::Record(NXLOG_ERROR, "You have not configured weighting function calculations!");
        return false;
    }
    if (!m_albedo_wf_included) {
        *numwf = static_cast<int>(m_wf_handles.size() *  m_wf_altitudes.size());
    }
    else {
        *numwf = static_cast<int>((m_wf_handles.size() - 1) *  m_wf_altitudes.size() + 1);
    }
    *numlos = static_cast<int>(m_linesofsight.NumRays()*m_num_stokes);
    *numwavel = static_cast<int>(m_wavelen.size());
    *wf = m_wf_storage.get();
    return true;
}

bool ISKEngine_Stub_DO::SetAtmosphericState(ISKClimatology_Stub* climatology)
{
    m_registered_calls[(unsigned int)TrackedCalls::AtmosphericState] = true;
    m_opticalstate.set_atmospheric_state(static_cast<skClimatology*>(climatology->RawObjectPointer()));
    return true;
}

bool ISKEngine_Stub_DO::SetPropertyScalar(const char* raw_property_name, double value)
{
    return setter(raw_property_name, m_scalar_set_functions, value);
}

bool ISKEngine_Stub_DO::SetPropertyArray(const char* raw_property_name, const double* value, int numpoints)
{
    return setter(raw_property_name, m_vector_set_functions, value, numpoints);
}

bool ISKEngine_Stub_DO::SetPropertyObject(const char* raw_property_name, nxUnknown* object)
{
    return setter(raw_property_name, m_object_set_functions, object);
}

bool ISKEngine_Stub_DO::SetPropertyString(const char* raw_property_name, const char* value)
{
    return setter(raw_property_name, m_string_set_functions, value);
}

bool ISKEngine_Stub_DO::GetProperty(const char* raw_property_specifier, const double** value, int* numpoints)
{
    try {
        std::string command(raw_property_specifier);
        int index1 = -1;
        int index2 = -1;
        std::string property_name;

        parseGetProperty(raw_property_specifier, property_name, index1, index2);

        // First check scalars
        auto sit = m_scalar_get_functions.find(property_name);
        if(sit != m_scalar_get_functions.end()) {
            m_getpropertybuffer.resize(1);
            sit->second(index1, index2, &(m_getpropertybuffer[0]));
            *value = &m_getpropertybuffer[0];
            *numpoints = 0;
            return true;
        }

        auto vit = m_vector_get_functions.find(property_name);
        if(vit != m_vector_get_functions.end()) {
            vit->second(index1, index2);
            *numpoints = (int) m_getpropertybuffer.size();
            *value = &m_getpropertybuffer[0];
            return true;
        } else {
            throw sktran_do_detail::InvalidConfiguration("Unknown property: " + property_name);
        }
    }
    catch(const std::exception& e) {
        nxLog::Record(NXLOG_ERROR, e.what());
        return false;
    }
    return true;
}

// DETAILS 

void ISKEngine_Stub_DO::parseGetProperty(const std::string& input,
                                         std::string& prop,
                                         int& index1,
                                         int& index2)
{
    boost::regex prop_regex("^(\\w+)");
    boost::regex index_regex("(\\d+)");
    boost::smatch m;

    // Get property should only be one!
    std::list<std::string> property_list;
    std::string input_search = input;
    while(boost::regex_search(input_search, m, prop_regex)) {
        property_list.push_back(m[0]);
        input_search = m.suffix().str();
        break; // we only take the first match
    }
    if(property_list.size() != 1) throw sktran_do_detail::InvalidConfiguration("Invalid property name syntax: " + input);
    prop = property_list.front();
    // Get indexes
    std::vector<int*> indexes = { &index1, &index2 };
    int i = 0;
    while(boost::regex_search(input_search, m, index_regex)) {
        if(i >= indexes.size()) throw sktran_do_detail::InvalidConfiguration("Too many indexes specified for get: " + prop);
        *indexes[i++] = atoi(m[0].str().c_str());
        input_search = m.suffix().str();
    }
}

void ISKEngine_Stub_DO::makeScalarSetFunctions()
{
    m_scalar_set_functions["averagereferencepoint"] = [&, this](double d) {
        m_force_single_reference_point = std::abs(d) > 1e-6;
    };

    m_scalar_set_functions["numstreams"] = [&, this](double d) {
        m_userspec.setNumberOfStreams(static_cast<sktran_do_detail::uint>(ceil(d - 0.5)));
    };

    m_scalar_set_functions["numlayers"] = [&, this](double d) {
        m_userspec.setNumberOfLayers(static_cast<int>(ceil(d - 0.5)));
    };

    m_scalar_set_functions["numthreads"] = [&, this](double d) {
        m_num_threads = static_cast<int>(ceil(d - 0.5));
    };

    m_scalar_set_functions["numwavelengthbatches"] = [&, this](double d) {
        m_wl_batches = static_cast<int>(ceil(d - 0.5));
    };

    m_scalar_set_functions["szarelseparation"] = [&, this](double d) {
        m_userspec.setSZARelSep(d);
    };

    m_scalar_set_functions["numbrdfexpansionterms"] = [&, this](double d) {
        unsigned int n = (unsigned int) ceil(d - 0.5);
        SKTRAN_DO_UserSpec::NTERMS_IN_EXPANSION nterms;
        switch(n) {
            case 64: 	nterms = SKTRAN_DO_UserSpec::NTERMS_IN_EXPANSION::u64;		break;
            case 128: 	nterms = SKTRAN_DO_UserSpec::NTERMS_IN_EXPANSION::u128;		break;
            case 256: 	nterms = SKTRAN_DO_UserSpec::NTERMS_IN_EXPANSION::u256;		break;
            case 512: 	nterms = SKTRAN_DO_UserSpec::NTERMS_IN_EXPANSION::u512;		break;
            case 1024: 	nterms = SKTRAN_DO_UserSpec::NTERMS_IN_EXPANSION::u1024;	break;
            default:
                throw sktran_do_detail::InvalidConfiguration("Invalid number of quadrature terms for brdf expansion! Must be 64, 128, 256, 512, or 1024.");
        }
        m_userspec.setNumBRDFExpansionTerms(nterms);
    };

    m_scalar_set_functions["diagnostics"] = [&, this](double d) {
        auto diag = static_cast<uint_fast16_t>(d) > 0;
        m_return_los_diagnostics = diag;
        m_return_rts_diagnostics = diag;
    };

    m_scalar_set_functions["losdiagnostics"] = [&, this](double d) {
        m_return_los_diagnostics = static_cast<uint_fast16_t>(d) > 0;
    };

    m_scalar_set_functions["rtsdiagnostics"] = [&, this](double d) {
        m_return_rts_diagnostics = static_cast<uint_fast16_t>(d) > 0;
    };

    m_scalar_set_functions["outputopticaldepth"] = [&, this](double d) {
        m_output_optical_depths = true;
    };

    m_scalar_set_functions["wfform"] = [&, this](double d) {
        bool didx = static_cast<uint_fast16_t>(d) > 0;
        SKTRAN_DO_UserSpec::WeightingFunctionForm form;
        if(didx) form = SKTRAN_DO_UserSpec::WeightingFunctionForm::dI_dX;
        else form = SKTRAN_DO_UserSpec::WeightingFunctionForm::dI_dLogX;
        m_userspec.setWFReturnForm(form);
    };

    m_scalar_set_functions["convergencecriteria"] = [&, this](double d) {
        m_userspec.setCauchyCriterion(d);
    };

    m_scalar_set_functions["forcenumazimuthterms"] = [&, this](double d) {
        m_userspec.setForceNumberAzimuthTerms(static_cast<sktran_do_detail::uint>(d));
    };

    m_scalar_set_functions["ssadither"] = [&, this](double d) {
        m_userspec.setSSAEqual1Dither(d);
    };

    m_scalar_set_functions["usepseudospherical"] = [&, this](double d) {
        bool use_ps = d > 0.5;

        m_userspec.setUsePsuedoSpherical(use_ps);
    };

    m_scalar_set_functions["usegreensfunction"] = [&, this](double d) {
        bool use = d > 0.5;

        m_userspec.setUseGreensFunction(use);
    };

    m_scalar_set_functions["uselosspherical"] = [&, this](double d) {
        bool use_ps = d > 0.5;

        m_userspec.setUseLOSSpherical(use_ps);
    };

    m_scalar_set_functions["useupwellingsphericalcorrection"] = [&, this](double d) {
        bool use_upwelling_spher = d > 0.5;

        m_userspec.setUseUpwellingSpherical(use_upwelling_spher);
    };

    m_scalar_set_functions["singlescatteronly"] = [&, this](double d) {
        bool ss_only = d > 0.5;

        m_userspec.setSSOnly(ss_only);
    };

    m_scalar_set_functions["numsphericalsza"] = [&, this](double d) {
        unsigned int n = (unsigned int)ceil(d - 0.5);

        m_userspec.setNumSZA(n);
    };

    m_scalar_set_functions["layerconstructionmethod"] = [&, this](double d) {
        unsigned int n = (unsigned int)ceil(d - 0.5);

        if (n == 0) {
            m_userspec.setLayerConstructionMethod(SKTRAN_DO_UserSpec::LayerConstructionMethod::uniform_pressure);
        }
        else if (n == 2) {
            m_userspec.setLayerConstructionMethod(SKTRAN_DO_UserSpec::LayerConstructionMethod::uniform_height);
        }
        else if (n == 3) {
            m_userspec.setLayerConstructionMethod(SKTRAN_DO_UserSpec::LayerConstructionMethod::match_altitudegrid);
        }
        else {
            throw sktran_do_detail::InvalidConfiguration("layerconstructionmethod must be either 0, 1, 2, or 3");
        }
    };

    m_scalar_set_functions["losgeoidadjustmentmethod"] = [&, this](double d) {
        unsigned int n = (unsigned int)ceil(d - 0.5);

        if (n == 0) {
            m_userspec.setLineOfSightAdjustment(SKTRAN_DO_UserSpec::LineOfSightAdjustment::translate_observer);
        }
        else if (n == 1) {
            m_userspec.setLineOfSightAdjustment(SKTRAN_DO_UserSpec::LineOfSightAdjustment::match_target_angles);
        }
        else {
            throw sktran_do_detail::InvalidConfiguration("layerconstructionmethod must be either 0 or 1");
        }
    };
}


void ISKEngine_Stub_DO::makeVectorSetFunctions()
{
    m_vector_set_functions["sun"] = [&, this](const double* sun, int n) {
        m_solar_position.FromSequence(sun);
        m_userspec.setSolarPosition(m_solar_position);
        m_registered_calls[(unsigned int)TrackedCalls::SunPosition] = true;
    };

    m_vector_set_functions["altitudegrid"] = [&, this](const double* altitudes, int n) {
        std::vector<double> alt_grid(altitudes, altitudes + n);
        m_userspec.setAltitudeGrid(alt_grid);
    };

    m_vector_set_functions["manuallayeraltitudes"] = [&, this](const double* altitudes, int n) {
        std::vector<double> alt_grid(altitudes, altitudes + n);
        m_userspec.setManualLayerAltitudes(alt_grid);
        m_userspec.setLayerConstructionMethod(SKTRAN_DO_UserSpec::LayerConstructionMethod::manual);
    };

    m_vector_set_functions["radiancetoa"] = [&, this](const double* sun, int n) {
        m_userspec.setTOAIntensities(sun[0]);
    };
    m_vector_set_functions["wfaltitudes"] = [&, this](const double* altitudes, int n) {
        m_wf_altitudes.resize(n);
        std::copy(altitudes, altitudes + n, m_wf_altitudes.begin());
        m_registered_calls[(unsigned int)TrackedCalls::WFAltitudes] = true;
    };
    m_vector_set_functions["wfwidths"] = [&, this](const double* widths, int n) {
        m_wf_widths_low.resize(n);
        m_wf_widths_high.resize(n);
        std::copy(widths, widths + n, m_wf_widths_low.begin());
        std::copy(widths, widths + n, m_wf_widths_high.begin());
    };
    m_vector_set_functions["wfwidthslow"] = [&, this](const double* widths, int n) {
        m_wf_widths_low.resize(n);
        std::copy(widths, widths + n, m_wf_widths_low.begin());
    };
    m_vector_set_functions["wfwidthshigh"] = [&, this](const double* widths, int n) {
        m_wf_widths_high.resize(n);
        std::copy(widths, widths + n, m_wf_widths_high.begin());
    };
    m_vector_set_functions["surfaceemissionwavelengths"] = [&, this](const double* wavel, int n) {
        std::vector<double> wavel_grid(wavel, wavel+n);
        m_userspec.setSurfaceEmissionWavelengths(wavel_grid);
    };
    m_vector_set_functions["surfaceemissionvalues"] = [&, this](const double* value, int n) {
        std::vector<double> value_grid(value, value + n);
        m_userspec.setSurfaceEmissionValues(value_grid);
    };
}

#define SKDO_EXCEPT_LOS_DIAG_NOT_SET sktran_do_detail::InvalidConfiguration("Property not availible. The 'losdiagnostics' property must have been set!")
#define SKDO_EXCEPT_RTS_DIAG_NOT_SET sktran_do_detail::InvalidConfiguration("Property not availible. The 'rtsdiagnostics' property must have been set!")
#define SKDO_EXCEPT_TOO_FEW_INDEXES sktran_do_detail::InvalidConfiguration("Not enough indexes were given!")
#define SKDO_EXCEPT_TOO_MANY_INDEXES sktran_do_detail::InvalidConfiguration("Too many indexes given!")
#define SKDO_EXCEPT_INDEX_OUT_OF_RANGE sktran_do_detail::InvalidConfiguration("An index was out-of-range!")

void ISKEngine_Stub_DO::makeScalarGetFunctions()
{
    m_scalar_get_functions["averagereferencepoint"] = [&, this](int, int, double* value) {
        *value = m_force_single_reference_point ? 1.0 : 0.0;
    };

    m_scalar_get_functions["numstreams"] = [&, this](int, int, double* value) {
        *value = static_cast<double>(m_userspec.getNumberOfStreams());
    };
    m_scalar_get_functions["numlayers"] = [&, this](int, int, double* value) {
        *value = static_cast<double>(m_userspec.getNumberOfLayers());
    };
    m_scalar_get_functions["numbrdfexpansionterms"] = [&, this](int, int, double* value) {
        *value = static_cast<double>(m_userspec.getNumBRDFQuadratureTerms());
    };

    m_scalar_get_functions["localsolarzenith"] = [&, this](int los_index, int not_used, double* value) {
        if(!m_return_los_diagnostics) throw SKDO_EXCEPT_LOS_DIAG_NOT_SET;
        if(not_used != -1) throw SKDO_EXCEPT_TOO_MANY_INDEXES;
        if(los_index == -1) throw SKDO_EXCEPT_TOO_FEW_INDEXES;
        if(los_index >= m_linesofsight.NumRays()) throw SKDO_EXCEPT_INDEX_OUT_OF_RANGE;
        *value = m_los_diagnostics.at(los_index).at(0).local_solar_zenith;
    };
    m_scalar_get_functions["localloszenith"] = [&, this](int los_index, int not_used, double* value) {
        if(!m_return_los_diagnostics) throw SKDO_EXCEPT_LOS_DIAG_NOT_SET;
        if(not_used != -1) throw SKDO_EXCEPT_TOO_MANY_INDEXES;
        if(los_index == -1) throw SKDO_EXCEPT_TOO_FEW_INDEXES;
        if(los_index >= m_linesofsight.NumRays()) throw SKDO_EXCEPT_INDEX_OUT_OF_RANGE;
        if(m_los_diagnostics.at(los_index).size() == 1) {
            *value = m_los_diagnostics.at(los_index).at(0).local_look_zentih;
        } else {
            *value = std::nan("1");
        }
    };
    m_scalar_get_functions["observerod"] = [&, this](int wav_index, int los_index, double* value) {
        *value = m_rts_diagnostics.at(wav_index).at(los_index).los_diagnostic.at(0).observer_od;
    };
    m_scalar_get_functions["localrelativeazimuth"] = [&, this](int los_index, int not_used, double* value) {
        if(!m_return_los_diagnostics) throw SKDO_EXCEPT_LOS_DIAG_NOT_SET;
        if(not_used != -1) throw SKDO_EXCEPT_TOO_MANY_INDEXES;
        if(los_index == -1) throw SKDO_EXCEPT_TOO_FEW_INDEXES;
        if(los_index >= m_linesofsight.NumRays()) throw SKDO_EXCEPT_INDEX_OUT_OF_RANGE;
        *value = m_los_diagnostics.at(los_index).at(0).local_relative_azimuth;
    };

    return;
}

void ISKEngine_Stub_DO::makeVectorGetFunctions()
{
    m_vector_get_functions["sun"] = [&, this](int, int) {
        m_getpropertybuffer.resize(3);
        const auto* sun = m_userspec.getSolarPosition(m_linesofsight.MeanMJD());
        if(sun == nullptr) throw sktran_do_detail::InvalidConfiguration("The sun has not been set!");
        m_getpropertybuffer[0] = sun->X();
        m_getpropertybuffer[1] = sun->Y();
        m_getpropertybuffer[2] = sun->Z();
    };
    // The following function require a line-of-sight index but no others
    m_vector_get_functions["referencepoint"] = [&, this](int los_index, int not_used) {
        if(!m_return_los_diagnostics) throw SKDO_EXCEPT_LOS_DIAG_NOT_SET;
        if(not_used != -1) throw SKDO_EXCEPT_TOO_MANY_INDEXES;
        if(los_index == -1) throw SKDO_EXCEPT_TOO_FEW_INDEXES;
        if(los_index >= m_linesofsight.NumRays()) throw SKDO_EXCEPT_INDEX_OUT_OF_RANGE;
        const auto& t = m_los_diagnostics.at(los_index);
        m_getpropertybuffer.resize(4);
        if(t.size() == 0) {
            std::fill_n(m_getpropertybuffer.begin(), 4, std::nan("1"));
        } else {
            // TODO: calculate the reference point
            GEODETIC_INSTANT pt;
            m_getpropertybuffer[0] = pt.latitude;
            m_getpropertybuffer[1] = pt.longitude;
            m_getpropertybuffer[2] = pt.heightm;
            m_getpropertybuffer[3] = pt.mjd;
        }
    };

    m_vector_get_functions["lostransmission"] = [&, this](int notused1, int notused2) {
        m_getpropertybuffer.resize(m_lostransmissionstorage.size());
        for (unsigned int i = 0; i < m_lostransmissionstorage.size(); ++i) {
            m_getpropertybuffer[i] = *(m_lostransmissionstorage.ArrayBasePtr() + i);
        }
    };

    m_vector_get_functions["layeropticaldepth"] = [&, this](int notused1, int notused2) {
        if(!m_output_optical_depths) throw SKDO_EXCEPT_RTS_DIAG_NOT_SET;
        m_getpropertybuffer.resize(m_odstorage.size());

        for (unsigned int i = 0; i < m_odstorage.size(); ++i) {
            m_getpropertybuffer[i] = *(m_odstorage.ArrayBasePtr() + i);
        }

    };

    m_vector_get_functions["layerssa"] = [&, this](int notused1, int notused2) {
        if(!m_output_optical_depths) throw SKDO_EXCEPT_RTS_DIAG_NOT_SET;
        m_getpropertybuffer.resize(m_ssastorage.size());

        for (unsigned int i = 0; i < m_ssastorage.size(); ++i) {
            m_getpropertybuffer[i] = *(m_ssastorage.ArrayBasePtr() + i);
        }
    };

    // The following function require both a wavelength index and a line-of-sight index
    m_vector_get_functions["daopticaldepths"] = [&, this](int wav_index, int los_index) {
        if(!m_return_rts_diagnostics) throw SKDO_EXCEPT_RTS_DIAG_NOT_SET;
        if(wav_index == -1 || los_index == -1) throw SKDO_EXCEPT_TOO_FEW_INDEXES;
        if(wav_index >= m_wavelen.size() || los_index >= m_linesofsight.NumRays()) throw SKDO_EXCEPT_INDEX_OUT_OF_RANGE;
        if(m_force_single_reference_point) los_index = 0;

        auto num_layers = m_userspec.getNumberOfLayers();
        m_getpropertybuffer.resize(num_layers);
        for(unsigned int i = 0; i < num_layers; ++i) {
            m_getpropertybuffer[i] = m_rts_diagnostics.at(wav_index).at(los_index).atmo_diagnostics.layer_optical_depths.at(i);
        }
    };
    m_vector_get_functions["dassa"] = [&, this](int wav_index, int los_index) {
        if(!m_return_rts_diagnostics) throw SKDO_EXCEPT_RTS_DIAG_NOT_SET;
        if(wav_index == -1 || los_index == -1) throw SKDO_EXCEPT_TOO_FEW_INDEXES;
        if(wav_index >= m_wavelen.size() || los_index >= m_linesofsight.NumRays()) throw SKDO_EXCEPT_INDEX_OUT_OF_RANGE;
        if(m_force_single_reference_point) los_index = 0;

        auto num_layers = m_userspec.getNumberOfLayers();
        m_getpropertybuffer.resize(num_layers);
        for(unsigned int i = 0; i < num_layers; ++i) {
            m_getpropertybuffer[i] = m_rts_diagnostics.at(wav_index).at(los_index).atmo_diagnostics.layer_ssa.at(i);
        }
    };
    m_vector_get_functions["daboundaryaltitudes"] = [&, this](int wav_index, int los_index) {
        if(!m_return_rts_diagnostics) throw SKDO_EXCEPT_RTS_DIAG_NOT_SET;
        if(wav_index == -1 || los_index == -1) throw SKDO_EXCEPT_TOO_FEW_INDEXES;
        if(wav_index >= m_wavelen.size() || los_index >= m_linesofsight.NumRays()) throw SKDO_EXCEPT_INDEX_OUT_OF_RANGE;
        if(m_force_single_reference_point) los_index = 0;

        auto num_layers = m_userspec.getNumberOfLayers();
        m_getpropertybuffer.resize(num_layers + 1);
        for(unsigned int i = 0; i < num_layers + 1; ++i) {
            m_getpropertybuffer[i] = m_rts_diagnostics.at(wav_index).at(los_index).atmo_diagnostics.layer_boundary_altitudes.at(i);
        }
    };
    m_vector_get_functions["daphasef"] = [&, this](int wav_index, int los_index) {
        if(!m_return_rts_diagnostics) throw SKDO_EXCEPT_RTS_DIAG_NOT_SET;
        if(wav_index == -1 || los_index == -1) throw SKDO_EXCEPT_TOO_FEW_INDEXES;
        if(wav_index >= m_wavelen.size() || los_index >= m_linesofsight.NumRays()) throw SKDO_EXCEPT_INDEX_OUT_OF_RANGE;
        if(m_force_single_reference_point) los_index = 0;

        auto num_layers = m_userspec.getNumberOfLayers();
        auto num_streams = m_userspec.getNumberOfStreams();
        m_getpropertybuffer.resize(num_layers * num_streams * num_streams);
        for(unsigned int i = 0; i < num_layers * num_streams * num_streams; ++i) {
            unsigned int layer_index = i / (num_streams*num_streams);
            unsigned int m = (i % (num_streams*num_streams)) / num_streams;
            unsigned int l = i % num_streams;
            m_getpropertybuffer[i] = m_rts_diagnostics.at(wav_index).at(los_index).atmo_diagnostics.layer_phasef_expansion.at(layer_index).at(m).at(l);
        }
    };
    m_vector_get_functions["radiancecomponents"] = [&, this](int wav_index, int los_index) {
        if(!m_return_rts_diagnostics) throw SKDO_EXCEPT_RTS_DIAG_NOT_SET;
        if(wav_index == -1 || los_index == -1) throw SKDO_EXCEPT_TOO_FEW_INDEXES;
        if(wav_index >= m_wavelen.size() || los_index >= m_linesofsight.NumRays()) throw SKDO_EXCEPT_INDEX_OUT_OF_RANGE;
        if(m_force_single_reference_point) los_index = 0;

        auto num_streams = m_userspec.getNumberOfStreams();
        m_getpropertybuffer.resize(num_streams);
        for(unsigned int i = 0; i < num_streams; ++i) {
            m_getpropertybuffer[i] = m_rts_diagnostics.at(wav_index).at(los_index).los_diagnostic.at(0).intensity_components.at(i);
        }
    };
    m_vector_get_functions["reflectedcomponents"] = [&, this](int wav_index, int los_index) {
        if(!m_return_rts_diagnostics) throw SKDO_EXCEPT_RTS_DIAG_NOT_SET;
        if(wav_index == -1 || los_index == -1) throw SKDO_EXCEPT_TOO_FEW_INDEXES;
        if(wav_index >= m_wavelen.size() || los_index >= m_linesofsight.NumRays()) throw SKDO_EXCEPT_INDEX_OUT_OF_RANGE;
        if(m_force_single_reference_point) los_index = 0;

        auto num_streams = m_userspec.getNumberOfStreams();
        m_getpropertybuffer.resize(num_streams);
        for(unsigned int i = 0; i < num_streams; ++i) {
            m_getpropertybuffer[i] = m_rts_diagnostics.at(wav_index).at(los_index).los_diagnostic.at(0).reflection_components.at(i);
        }
    };
    m_vector_get_functions["participatingsourceterms"] = [&, this](int wav_index, int los_index) {
        if(!m_return_rts_diagnostics) throw SKDO_EXCEPT_RTS_DIAG_NOT_SET;
        if(wav_index == -1 || los_index == -1) throw SKDO_EXCEPT_TOO_FEW_INDEXES;
        if(wav_index >= m_wavelen.size() || los_index >= m_linesofsight.NumRays()) throw SKDO_EXCEPT_INDEX_OUT_OF_RANGE;
        if(m_force_single_reference_point) los_index = 0;

        auto num_layers = m_userspec.getNumberOfLayers();
        auto num_streams = m_userspec.getNumberOfStreams();
        m_getpropertybuffer.resize(num_layers * num_streams);
        for(unsigned int i = 0; i < num_layers * num_streams; ++i) {
            unsigned int layer_index = i / num_streams;
            unsigned int m = i % num_streams;
            m_getpropertybuffer[i] = m_rts_diagnostics.at(wav_index).at(los_index).los_diagnostic.at(0).layer_participating_source_term.at(layer_index).at(m);
        }
    };
    m_vector_get_functions["wfaltitudes"] = [&, this](int wav_index, int los_index) {
        if(!m_return_rts_diagnostics) throw SKDO_EXCEPT_RTS_DIAG_NOT_SET;
        if(wav_index == -1 || los_index == -1) throw SKDO_EXCEPT_TOO_FEW_INDEXES;
        if(wav_index >= m_wavelen.size() || los_index >= m_linesofsight.NumRays()) throw SKDO_EXCEPT_INDEX_OUT_OF_RANGE;
        if(m_force_single_reference_point) los_index = 0;

        m_getpropertybuffer.resize(m_rts_diagnostics.at(wav_index).at(los_index).ptrb_altitudes.size());
        std::copy(m_rts_diagnostics.at(wav_index).at(los_index).ptrb_altitudes.cbegin(),
                  m_rts_diagnostics.at(wav_index).at(los_index).ptrb_altitudes.cend(),
                  m_getpropertybuffer.begin()
        );
    };
    m_vector_get_functions["wfwidths"] = [&, this](int wav_index, int los_index) {
        if(!m_return_rts_diagnostics) throw SKDO_EXCEPT_RTS_DIAG_NOT_SET;
        if(wav_index == -1 || los_index == -1) throw SKDO_EXCEPT_TOO_FEW_INDEXES;
        if(wav_index >= m_wavelen.size() || los_index >= m_linesofsight.NumRays()) throw SKDO_EXCEPT_INDEX_OUT_OF_RANGE;
        if(m_force_single_reference_point) los_index = 0;

        m_getpropertybuffer.resize(m_rts_diagnostics.at(wav_index).at(los_index).ptrb_heights.size());
        std::copy(m_rts_diagnostics.at(wav_index).at(los_index).ptrb_heights.cbegin(),
                  m_rts_diagnostics.at(wav_index).at(los_index).ptrb_heights.cend(),
                  m_getpropertybuffer.begin()
        );
    };
    m_vector_get_functions["ptrbssaqty"] = [&, this](int wav_index, int los_index) {
        if(!m_return_rts_diagnostics) throw SKDO_EXCEPT_RTS_DIAG_NOT_SET;
        if(wav_index == -1 || los_index == -1) throw SKDO_EXCEPT_TOO_FEW_INDEXES;
        if(wav_index >= m_wavelen.size() || los_index >= m_linesofsight.NumRays()) throw SKDO_EXCEPT_INDEX_OUT_OF_RANGE;
        if(m_force_single_reference_point) los_index = 0;

        m_getpropertybuffer.resize(m_rts_diagnostics.at(wav_index).at(los_index).ptrb_ssa_qty.size());
        std::copy(m_rts_diagnostics.at(wav_index).at(los_index).ptrb_ssa_qty.cbegin(),
                  m_rts_diagnostics.at(wav_index).at(los_index).ptrb_ssa_qty.cend(),
                  m_getpropertybuffer.begin()
        );
    };
    m_vector_get_functions["ptrboptdqty"] = [&, this](int wav_index, int los_index) {
        if(!m_return_rts_diagnostics) throw SKDO_EXCEPT_RTS_DIAG_NOT_SET;
        if(wav_index == -1 || los_index == -1) throw SKDO_EXCEPT_TOO_FEW_INDEXES;
        if(wav_index >= m_wavelen.size() || los_index >= m_linesofsight.NumRays()) throw SKDO_EXCEPT_INDEX_OUT_OF_RANGE;
        if(m_force_single_reference_point) los_index = 0;

        m_getpropertybuffer.resize(m_rts_diagnostics.at(wav_index).at(los_index).ptrb_optd_qty.size());
        std::copy(m_rts_diagnostics.at(wav_index).at(los_index).ptrb_optd_qty.cbegin(),
                  m_rts_diagnostics.at(wav_index).at(los_index).ptrb_optd_qty.cend(),
                  m_getpropertybuffer.begin()
        );
    };
}

void ISKEngine_Stub_DO::makeObjectSetFunctions()
{
    return; /* ** NOT TESTED **
			m_objectsetfunctions["albedoclim"] = [&, this](nxUnknown* obj) {
			skClimatology* clim = static_cast<skClimatology*>(obj);
			if(clim != nullptr) {
			skBRDF_AlbedoPlane* alb = new skBRDF_AlbedoPlane(clim);
			m_opticalstate.SetAlbedoObject(alb);
			m_registered_calls[TrackedCalls::BRDF] = true;
			} else {
			throw sasktran_disco::InvalidConfiguration("Invalid climatology object!");
			}
			};
			return true;*/
}

void ISKEngine_Stub_DO::makeStringSetFunctions() {
    m_string_set_functions["wfclimatologies"] = [&, this](const char* cstr) {
        m_wf_handles.clear();
        std::string s(cstr);
        boost::regex e("(\\w+)");
        boost::smatch m;

        while(boost::regex_search(s, m, e)) {
            std::string handle_str = m[0];
            CLIMATOLOGY_HANDLE* clim;
            nxString nxs(handle_str.c_str());
            nxs.MakeLower();
            // Special cases for albedo
            if (nxs == "albedo" || nxs == "brdf")
            {
                clim = &SKCLIMATOLOGY_ALBEDO;
            }
            else {
                clim = FindGlobalClimatologyHandle(handle_str.c_str());
            }

            if(*clim == SKCLIMATOLOGY_UNDEFINED) {
                std::string err = "Unknown climatology when configuring weighting functions: " + handle_str;
                throw sktran_do_detail::InvalidConfiguration(err);
            }
            m_wf_handles.push_back(*clim);
            s = m.suffix().str();
        }
    };
}