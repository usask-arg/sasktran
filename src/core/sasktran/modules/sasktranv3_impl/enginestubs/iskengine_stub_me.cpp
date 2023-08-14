#include "../dllimplementation/stdafx.h"
#include "modules/sasktranv3_impl/sktranif_impl_helperclasses.h"
#include <boost/regex.hpp>

ISKEngine_Stub_ME::ISKEngine_Stub_ME() {
    m_nstokes = 1;

    m_altitude_grid.resize(101);
    for (int i = 0; i < 101; ++i) {
        m_altitude_grid(i) = i * 1000;
    }

    m_geodetic.SelectGeoid(m_geodetic.WGS84);

    MakeScalarSetFunctions();
    MakeVectorSetFunctions();
    MakeObjectSetFunctions();
    MakeVectorGetFunctions();
    MakeScalarGetFunctions();
    MakeStringSetFunctions();
}

bool ISKEngine_Stub_ME::MakeScalarSetFunctions() {
    AddSetScalarFunction("msmode",
                         [&, this](double d ) {
                             int specifier = (int) ceil(d - 0.5);
                             bool ok = true;

                             if(specifier == 0) {
                                 m_config.set_multiple_scatter_source(sasktran2::Config::MultipleScatterSource::none);
                             }

                             if(specifier == 1) {
                                 m_config.set_multiple_scatter_source(sasktran2::Config::MultipleScatterSource::discrete_ordinates);
                             }

                             if(specifier == 2) {
                                 m_config.set_multiple_scatter_source(sasktran2::Config::MultipleScatterSource::hr);
                             }

                             return ok;
                         }
                         );

    AddSetScalarFunction("BASE",
                         [&, this](double d ) {
                             int specifier = (int) ceil(d - 0.5);
                             bool ok = true;

                             return ok;
                         }
    );

    AddSetScalarFunction("nstokes",
                         [&, this](double d ) {
                             int specifier = (int) ceil(d - 0.5);
                             bool ok = true;

                             m_nstokes = specifier;

                             return ok;
                         }
    );

    AddSetScalarFunction("numdostreams",
                         [&, this](double d ) {
                             int specifier = (int) ceil(d - 0.5);
                             bool ok = true;

                             m_config.set_num_do_streams(specifier);

                             return ok;
                         }
    );

    AddSetScalarFunction("numssmoments",
                         [&, this](double d ) {
                             int specifier = (int) ceil(d - 0.5);
                             bool ok = true;

                             m_config.set_num_singlescatter_moments(specifier);

                             return ok;
                         }
    );

    AddSetScalarFunction("numthreads",
                         [&, this](double d ) {
                             int specifier = (int) ceil(d - 0.5);
                             bool ok = true;

                             m_config.set_num_threads(specifier);

                             return ok;
                         }
    );

    AddSetScalarFunction("initializehrwithdo",
                         [&, this](double d ) {
                             int specifier = (int) ceil(d - 0.5);
                             bool ok = true;

                             if(specifier == 1) {
                                 m_config.set_initialize_hr_with_do(true);
                             } else {
                                 m_config.set_initialize_hr_with_do(false);
                             }

                             return ok;
                         }
    );

    AddSetScalarFunction("applydeltascaling",
                         [&, this](double d ) {
                             int specifier = (int) ceil(d - 0.5);
                             bool ok = true;

                             if(specifier == 1) {
                                 m_config.set_apply_delta_scaling(true);
                             } else {
                                 m_config.set_apply_delta_scaling(false);
                             }

                             return ok;
                         }
    );

    AddSetScalarFunction("numhriterations",
                         [&, this](double d ) {
                             int specifier = (int) ceil(d - 0.5);
                             bool ok = true;

                             m_config.set_num_hr_spherical_iterations(specifier);

                             return ok;
                         }
    );

    AddSetScalarFunction("numhrincoming",
                         [&, this](double d ) {
                             int specifier = (int) ceil(d - 0.5);
                             bool ok = true;

                             m_config.set_num_hr_incoming(specifier);

                             return ok;
                         }
    );

    AddSetScalarFunction("numhroutgoing",
                         [&, this](double d ) {
                             int specifier = (int) ceil(d - 0.5);
                             bool ok = true;

                             m_config.set_num_hr_outgoing(specifier);

                             return ok;
                         }
    );

    return true;
}

bool ISKEngine_Stub_ME::MakeVectorSetFunctions() {
    AddSetVectorFunction( "setsun",
                          [&, this](const double* sun, int n)
                          {
                              bool ok;
                              ok = ( n == 3 );
                              if (!ok)
                              {
                                  BOOST_LOG_TRIVIAL(error) << "ISKEngine CO, SetSun is not the correct size";
                              }

                              m_manual_sun = std::make_unique<nxVector>(sun[0], sun[1], sun[2]);

                              return ok;
                          }
    );

    AddSetVectorFunction( "altitudegrid",
                          [&, this](const double* alts, int n)
                          {
                              bool ok = true;

                              m_altitude_grid.resize(n);
                              for(int i = 0; i < n; ++i) {
                                  m_altitude_grid(i) = alts[i];
                              }

                              return ok;
                          }
    );

    return true;
}

bool ISKEngine_Stub_ME::MakeObjectSetFunctions() {
    return true;

}

bool ISKEngine_Stub_ME::MakeVectorGetFunctions() {
    AddGetVectorFunction( "referencepoint",
                          [&, this]( int index )
                          {
                              GEODETIC_INSTANT pt = m_geometry_constructor.reference_point();
                              m_getpropertybuffer.resize(4);
                              m_getpropertybuffer[0] = pt.latitude;
                              m_getpropertybuffer[1] = pt.longitude;
                              m_getpropertybuffer[2] = pt.heightm;
                              m_getpropertybuffer[3] = pt.mjd;
                              return true;
                          }
    );

    AddGetVectorFunction("ssa",
                          [&, this]( int wavelindex )
                          {
                              auto& atmosphere = dynamic_cast<sktran_me::AtmosphereConstructor<1>*>(m_atmosphere_constructor.get())->atmosphere();

                              int num_grid = atmosphere.storage().ssa.rows();

                              m_getpropertybuffer.resize(num_grid);

                              for(int i = 0; i < num_grid; ++i) {
                                  m_getpropertybuffer[i] = atmosphere.storage().ssa(i, wavelindex);
                              }

                            return true;
                          }
    );

    AddGetVectorFunction("extinction",
                         [&, this]( int wavelindex )
                         {
                             auto& atmosphere = dynamic_cast<sktran_me::AtmosphereConstructor<1>*>(m_atmosphere_constructor.get())->atmosphere();

                             int num_grid = atmosphere.storage().total_extinction.rows();

                             m_getpropertybuffer.resize(num_grid);

                             for(int i = 0; i < num_grid; ++i) {
                                 m_getpropertybuffer[i] = atmosphere.storage().total_extinction(i, wavelindex);
                             }

                             return true;
                         }
    );

    return true;

}

bool ISKEngine_Stub_ME::MakeScalarGetFunctions() {
    return true;

}

bool ISKEngine_Stub_ME::MakeStringSetFunctions() {
    AddSetStringFunction("wfspecies",
                         [&, this](const char* cstr)
                         {
                             std::vector<std::string> wf_handles;
                             std::string s(cstr);
                             boost::regex e("(\\w+)");
                             boost::smatch m;

                             while (boost::regex_search(s, m, e)) {
                                 std::string handle_str = m[0];
                                 wf_handles.push_back(handle_str);
                                 s = m.suffix().str();
                             }

                             m_wfhandler.set_from_strings(wf_handles);

                             return true;
                         }
    );
    return true;
}

bool ISKEngine_Stub_ME::AddLineOfSight(double mjd, const nxVector& observer, const nxVector& lookvector, int* losindex) {
    bool		ok;

    ok = m_linesofsight.AddLineOfSight( observer, lookvector, mjd );
    if (ok)
    {
        *losindex = (int)m_linesofsight.NumRays() - 1;
    }
    else
    {
        *losindex = -999999;
    }
    m_geometry_is_configured = false;

    return ok;
}

bool ISKEngine_Stub_ME::AddSpecies(const CLIMATOLOGY_HANDLE& species, ISKClimatology_Stub* climatology, ISKOpticalProperty_Stub* opticalproperty) {
    // Add to atmosphere inteface
    m_atmosphere_interface.add_species(species, dynamic_cast<skClimatology*>(climatology->RawObjectPointer()), dynamic_cast<skOpticalProperties*>(opticalproperty->RawObjectPointer()));

	return true;
}

bool ISKEngine_Stub_ME::AddEmission(const EMISSION_HANDLE& species, ISKEmission_Stub* emission) {
    // TODO: Implement when emissions work

	return true;

}

bool ISKEngine_Stub_ME::SetAtmosphericState(ISKClimatology_Stub* climatology) {
    m_atmosphere_interface.set_atmospheric_state(dynamic_cast<skClimatology*>(climatology->RawObjectPointer()));

	return true;
}

bool ISKEngine_Stub_ME::SetAlbedo(double albedo) {
    // Add to atmosphere inteface
    m_atmosphere_interface.set_albedo(albedo);

	return true;
}

bool ISKEngine_Stub_ME::SetBRDF(ISKBrdf_Stub* brdf) {
    // Add to atmosphere inteface
    m_atmosphere_interface.set_albedo(dynamic_cast<skBRDF*>(brdf->RawObjectPointer()));

	return true;
}

bool ISKEngine_Stub_ME::SetPolarizationMode(int polarizationmode) {
    m_nstokes = polarizationmode;

    if(m_nstokes == 1 || m_nstokes == 3) {
        return true;
    }

    BOOST_LOG_TRIVIAL(error) << "EngineME only supports polarization modes 1 and 3, which correspond to the number of stokes vector elements included";

	return false;
}

bool ISKEngine_Stub_ME::SetWavelengths(const double* wavelen, int numwavelen) {
    m_wavelengths.assign(wavelen, wavelen + numwavelen);

	return true;
}

bool ISKEngine_Stub_ME::InitializeModel() {
    if(m_geometry_is_configured) {
        // Geometry hasn't changed, and we have already been here, so just return
        return true;
    }

    // TODO: Get the sun unit vector
    if(!m_manual_sun) {
        BOOST_LOG_TRIVIAL(error) << "Error, sun must be manually specified with SetSun for EngineCO";
    }

    nxVector sun_unit = *m_manual_sun;

    // Construct the viewing and geometry objects
    m_geometry_constructor.construct_from_config(m_config, m_altitude_grid, m_geodetic, m_linesofsight, m_viewing_rays, m_geometry, sun_unit);
    m_wfhandler.set_geometry(*m_geometry);


    // Create the internal engine object and initialize it
    if(m_nstokes == 1) {
        m_engine = std::make_unique<Sasktran2<1>>(m_config, m_geometry.get(), *m_viewing_rays);
    } else if (m_nstokes == 3) {
        m_engine = std::make_unique<Sasktran2<3>>(m_config, m_geometry.get(), *m_viewing_rays);
    } else {
        // Should be unreachable
        BOOST_LOG_TRIVIAL(error) << "m_nstokes is not 1 or 3";
    }

    m_geometry_is_configured = true;

	return true;

}


bool ISKEngine_Stub_ME::CalculateRadiance(const double** radiance, int* numwavelens, int* numlinesofsight) {
    // Create the internal engine object, atmosphere, and do the calculation
    if(m_nstokes == 1) {
        m_atmosphere_constructor = std::make_unique<sktran_me::AtmosphereConstructor<1>>(m_atmosphere_interface, m_wfhandler);

        // Construction of atmosphere/output
        dynamic_cast<sktran_me::AtmosphereConstructor<1>*>(m_atmosphere_constructor.get())->construct_atmospheric_state(*m_geometry,
                                                                                                                        m_config,
                                                                                                                        *m_viewing_rays,
                                                                                                                        m_geometry_constructor.reference_point(),
                                                                                                                        m_wavelengths);


        dynamic_cast<Sasktran2<1>*>(m_engine.get())->calculate_radiance(dynamic_cast<sktran_me::AtmosphereConstructor<1>*>(m_atmosphere_constructor.get())->atmosphere(),
                                                                        dynamic_cast<sktran_me::AtmosphereConstructor<1>*>(m_atmosphere_constructor.get())->output());
    } else if (m_nstokes == 3) {
        m_atmosphere_constructor = std::make_unique<sktran_me::AtmosphereConstructor<3>>(m_atmosphere_interface, m_wfhandler);

        // Construction of atmosphere/output
        dynamic_cast<sktran_me::AtmosphereConstructor<3>*>(m_atmosphere_constructor.get())->construct_atmospheric_state(*m_geometry,
                                                                                                                        m_config,
                                                                                                                        *m_viewing_rays,
                                                                                                                        m_geometry_constructor.reference_point(),
                                                                                                                        m_wavelengths);

        dynamic_cast<Sasktran2<3>*>(m_engine.get())->calculate_radiance(dynamic_cast<sktran_me::AtmosphereConstructor<3>*>(m_atmosphere_constructor.get())->atmosphere(),
                                                                        dynamic_cast<sktran_me::AtmosphereConstructor<3>*>(m_atmosphere_constructor.get())->output());
    } else {
        // Should be unreachable
        BOOST_LOG_TRIVIAL(error) << "m_nstokes is not 1 or 3";
    }

    m_atmosphere_constructor->assign_output_radiance(radiance, numwavelens, numlinesofsight);

	return true;
}

bool ISKEngine_Stub_ME::CalculateStokesVector(const ISKStokesVector** radiancep, int* numwavelens, int* numlinesofsight) {
    BOOST_LOG_TRIVIAL(error) << "CalculateStokesVector should not be used for SasktranME";

	return false;
}

bool ISKEngine_Stub_ME::GetWeightingFunctions(const double** wf, int* numwavel, int* numlinesofsight, int* numwf) {
    *numwavel = m_wavelengths.size();
    *numlinesofsight = m_viewing_rays->observer_rays().size();
    *numwf = m_wfhandler.num_output_wf();

    m_atmosphere_constructor->assign_wf_buffer(wf);

	return true;
}
