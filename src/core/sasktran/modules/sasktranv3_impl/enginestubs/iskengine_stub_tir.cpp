#include "../dllimplementation/stdafx.h"
#include "modules/sasktranv3_impl/sktranif_impl_helperclasses.h"

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::ISKEngine_Stub_TIR			2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

ISKEngine_Stub_TIR::ISKEngine_Stub_TIR()
{
    m_radiance.SetReuseMemory(true);
    m_numthreads = 0;
    m_modelisconfigured = false;
    m_usecache = false;
    m_cacheisvalid = false;
    MakeScalarSetFunctions();
    MakeVectorSetFunctions();
    MakeScalarGetFunctions();
    MakeVectorGetFunctions();
    MakeStringSetFunctions();
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::~ISKEngine_Stub_TIR			2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

ISKEngine_Stub_TIR::~ISKEngine_Stub_TIR()
{
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::AddLineOfSight				2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::AddLineOfSight(double mjd, const nxVector& obs, const nxVector& look, int* losindex)
{
    bool ok;

    ok = m_linesofsight.AddLineOfSight(obs, look, mjd);
    if (ok)
    {
        *losindex = (int)m_linesofsight.NumRays() - 1;
    }
    else
    {
        *losindex = -999999;
    }
    m_modelisconfigured = false;
    m_cacheisvalid = false;
    m_radiance.erase();
    return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::CheckModelNotInitialized	2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool  ISKEngine_Stub_TIR::CheckModelNotInitialized(const char* propertystr) const
{
    if (m_modelisconfigured)  nxLog::Record(NXLOG_WARNING, "ISKEngine TIR, Cannot change property %s after the geometry has been initialized", (const char*)propertystr);
    return !m_modelisconfigured;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::MakeScalarSetFunctions		2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::MakeScalarSetFunctions()
{
    m_scalarsetfunctions[nxString("usecache")] = [&, this](double d)
    {
        int specifier = (int)ceil(d - 0.5);
        bool ok = true;

        if (specifier == 0)
        {
            m_usecache = false;
        }
        else if (specifier == 1)
        {
            m_usecache = true;
        }
        else
        {
            ok = false;
            nxLog::Record(NXLOG_WARNING, "ISKEngine TIR, Unknown specifier (%d) for property usecache", (int)specifier);
        }

        return ok;
    };

    m_scalarsetfunctions[nxString("calctemperaturewf")] = [&, this](double d)
    {
        int specifier = (int)ceil(d - 0.5);
        bool ok = true;

        ok = CheckModelNotInitialized("calctemperaturewf");
        if (ok)
        {
            if (specifier == 0)
            {
                m_specs.WeightingFunctionSpecs().SetDoTemperatureWF(false);
            }
            else if (specifier == 1)
            {
                m_specs.WeightingFunctionSpecs().SetDoTemperatureWF(true);
            }
            else
            {
                ok = false;
                nxLog::Record(NXLOG_WARNING, "ISKEngine TIR, Unknown specifier (%d) for property calctemperaturewf", (int)specifier);
            }
        }

        return ok;
    };

    m_scalarsetfunctions[nxString("userefraction")] = [&, this](double d)
    {
        int specifier = (int)ceil(d - 0.5);
        bool ok = true;

        ok = CheckModelNotInitialized("userefraction");
        if (ok)
        {
            if (specifier == 0)
            {
                ok = ok && m_specs.RayTracingSpecs().SetLinesOfSightType(RayTracerTypeTIR::straight);
            }
            else if (specifier == 1)
            {
                ok = ok && m_specs.RayTracingSpecs().SetLinesOfSightType(RayTracerTypeTIR::curved);
            }
            else
            {
                ok = false;
                nxLog::Record(NXLOG_WARNING, "ISKEngine TIR, Unknown specifier (%d) for property userefraction", (int)specifier);
            }
        }
        return ok;
    };

    m_scalarsetfunctions[nxString("opticalpropertiesheightres")] = [&, this](double d)
    {
        bool ok;

        ok = CheckModelNotInitialized("opticalpropertiesheightres");
        if (ok)
        {
            m_specs.OpticalPropertiesSpecs().SetHeightRes(d);
            if (!ok)
            {
                nxLog::Record(NXLOG_WARNING, "ISKEngine TIR, Error setting optical properties height spacing to %e", (double)d);
            }
        }
        return ok;
    };

    m_scalarsetfunctions[nxString("useadaptiveintegration")] = [&, this](double d)
    {
        int specifier = (int)ceil(d - 0.5);
        bool ok = true;

        ok = CheckModelNotInitialized("useadaptiveintegration");
        if (ok)
        {
            if (specifier == 0)
            {
                m_specs.IntegratorSpecs().SetOpticalPropertiesType(OpticalPropertiesIntegratorTypeTIR::straight);
            }
            else if (specifier == 1)
            {
                m_specs.IntegratorSpecs().SetOpticalPropertiesType(OpticalPropertiesIntegratorTypeTIR::adaptive);
            }
            else
            {
                ok = false;
                nxLog::Record(NXLOG_WARNING, "ISKEngine TIR, Unknown specifier (%d) for property useadaptiveintegration", (int)specifier);
            }
        }
        return ok;
    };

    m_scalarsetfunctions[nxString("uselinearextinction")] = [&, this](double d)
    {
        int specifier = (int)ceil(d - 0.5);
        bool ok = true;

        ok = CheckModelNotInitialized("uselinearextinction");
        if (ok)
        {
            if (specifier == 0)
            {
                m_specs.IntegratorSpecs().SetLayerExtinctionType(LayerExtinctionTypeTIR::constant);
            }
            else if (specifier == 1)
            {
                m_specs.IntegratorSpecs().SetLayerExtinctionType(LayerExtinctionTypeTIR::linearwithheight);
            }
            else
            {
                ok = false;
                nxLog::Record(NXLOG_WARNING, "ISKEngine TIR, Unknown specifier (%d) for property uselinearextinction", (int)specifier);
            }
        }
        return ok;
    };

    m_scalarsetfunctions[nxString("usevmrwfunit")] = [&, this](double d)
    {
        int specifier = (int)ceil(d - 0.5);
        bool ok = true;

        ok = CheckModelNotInitialized("usevmrwfunit");
        if (ok)
        {
            if (specifier == 0)
            {
                m_specs.WeightingFunctionSpecs().SetSpeciesWFUnit(SpeciesWFUnitTIR::numberdensity);
            }
            else if (specifier == 1)
            {
                m_specs.WeightingFunctionSpecs().SetSpeciesWFUnit(SpeciesWFUnitTIR::vmr);
            }
            else
            {
                ok = false;
                nxLog::Record(NXLOG_WARNING, "ISKEngine TIR, Unknown specifier (%d) for property usevmrwfunit", (int)specifier);
            }
        }
        return ok;
    };

    m_scalarsetfunctions[nxString("sourcetermorder")] = [&, this](double d)
    {
        int specifier = (int)ceil(d - 0.5);
        bool ok = true;

        ok = CheckModelNotInitialized("sourcetermorder");
        if (ok)
        {
            if (specifier == 0)
            {
                m_specs.IntegratorSpecs().SetSourceTermType(SourceTermOrderTIR::order0);
            }
            else if (specifier == 2)
            {
                m_specs.IntegratorSpecs().SetSourceTermType(SourceTermOrderTIR::order2);
            }
            else
            {
                ok = false;
                nxLog::Record(NXLOG_WARNING, "ISKEngine TIR, Source term order (%d) is not supported; avaible options are 0 or 2", (int)specifier);
            }
        }
        return ok;
    };

    m_scalarsetfunctions[nxString("surfaceheight")] = [&, this](double d)
    {
        bool ok;

        ok = CheckModelNotInitialized("surfaceheight");
        if (ok)
        {
            m_specs.RayTracingSpecs().SetGroundShiftAlt(d);
            if (!ok) nxLog::Record(NXLOG_WARNING, "ISKEngine TIR, Error setting property surfaceheight to %e", (double)d);
        }
        return ok;
    };

    m_scalarsetfunctions[nxString("opticaltabledimensions")] = [&, this](double d)
    {
        int specifier = (int)ceil(d - 0.5);
        bool ok;

        ok = CheckModelNotInitialized("opticaltabledimensions");
        if (ok)
        {
            if (specifier == 1)
            {
                m_specs.OpticalPropertiesSpecs().SetOpticalPropertiesTableDimension(AtmosphereDimensionTIR::dim1);
            }
            else if (specifier == 2)
            {
                m_specs.OpticalPropertiesSpecs().SetOpticalPropertiesTableDimension(AtmosphereDimensionTIR::dim2);
            }
            else
            {
                ok = false;
                nxLog::Record(NXLOG_WARNING, "ISKEngine TIR, Invalid number of dimensions (%d) for property opticaltabledimensions", (int)specifier);
            }
        }
        return ok;
    };

    m_scalarsetfunctions[nxString("raytracingshells")] = [&, this](double d)
    {
        bool ok;

        ok = CheckModelNotInitialized("raytracingshells");
        if (ok)
        {
            ok = m_specs.RayTracingSpecs().SetShellSpacing(d);
            if (!ok)
            {
                nxLog::Record(NXLOG_WARNING, "ISKEngine TIR, Error setting shell spacing to %e", (double)d);
            }
        }
        return ok;
    };

    m_scalarsetfunctions[nxString("maxopticaldepthofcell")] = [&, this](double d)
    {
        bool ok;

        ok = CheckModelNotInitialized("maxopticaldepthofcell");
        if (ok)
        {
            ok = (d > 0);
            if (ok)
            {
                m_specs.IntegratorSpecs().SetMaxOpticalDepth(d);
            }
            else
            {
                nxLog::Record(NXLOG_WARNING, "ISKEngine TIR, Invalid Maximum Optical Depth of Cell [%e] entered for property maxopticaldepthofcell", d);
                ok = false;
            }
        }
        return ok;
    };

    m_scalarsetfunctions[nxString("minextinctionratioofcell")] = [&, this](double d)
    {
        bool ok;

        ok = CheckModelNotInitialized("minextinctionratioofcell");
        if (ok)
        {
            ok = (d > 0 && d <= 1.0);
            if (ok)
            {
                m_specs.IntegratorSpecs().SetMinExtinctionRatio(d);
            }
            else
            {
                nxLog::Record(NXLOG_WARNING, "ISKEngine TIR, Invalid minimum extinction ratio of cell [%e] entered", (double)d);
            }
        }
        return ok;
    };

    m_scalarsetfunctions[nxString("groundemissivity")] = [&, this](double d)
    {
        bool ok;

        ok = CheckModelNotInitialized("groundemissivity");
        if (ok)
        {
            ok = (d >= 0);
            if (ok)
            {
                m_specs.OpticalPropertiesSpecs().SetGroundEmissivity(d);
            }
            else
            {
                nxLog::Record(NXLOG_WARNING, "ISKEngine TIR, Invalid surface emissivity [%e] entered; can not be negative", (double)d);
            }
        }
        return ok;
    };

    m_scalarsetfunctions[nxString("numthreads")] = [&, this](double d)
    {
        int numthreads = (int)ceil(d - 0.5);
        bool ok;


        ok = CheckModelNotInitialized("numthreads");
        if (ok)
        {
            ok = (numthreads >= 0 && numthreads < 1000);
            if (!ok)
            {
                nxLog::Record(NXLOG_WARNING, "ISKEngine TIR, Invalid number of threads (%d) for property numthreads. Enter a number between 0 and 1000 inclusive", (int)numthreads);
            }
            else
            {
                m_numthreads = numthreads;
            }
        }
        return ok;
    };

    m_scalarsetfunctions[nxString("toaheight")] = [&, this](double d)
    {
        bool ok;

        ok = CheckModelNotInitialized("toaheight");
        if (ok)
        {
            ok = m_specs.RayTracingSpecs().SetTOAHeight(d);
            if (!ok) nxLog::Record(NXLOG_WARNING, "ISKEngine TIR, Error setting property toaheight to %e", (double)d);
        }
        return ok;
    };

    return true;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::MakeVectorSetFunctions		2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::MakeVectorSetFunctions()
{
    m_vectorsetfunctions[nxString("wavelengths")] = [&, this](const double* h, int n)
    {
        bool ok;

        ok = CheckModelNotInitialized("wavelengths");
        if (ok)
        {
            ok = ok && SetWavelengths(h, n);
        }
        return ok;
    };

    m_vectorsetfunctions[nxString("wfheights")] = [&, this](const double* h, int n)
    {
        bool ok;

        ok = CheckModelNotInitialized("wfheights");
        if (ok)
        {
            m_specs.WeightingFunctionSpecs().SetWeightingFunctionHeight(std::vector<double>(h, h + n));
        }
        return ok;
    };

    m_vectorsetfunctions[nxString("wfwidths")] = [&, this](const double* h, int n)
    {
        bool ok;

        ok = CheckModelNotInitialized("wfwidths");
        if (ok)
        {
            m_specs.WeightingFunctionSpecs().SetWeightingFunctionWidth(std::vector<double>(h, h + n));
        }
        return ok;
    };

    return true;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::MakeScalarGetFunctions		2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::MakeScalarGetFunctions()
{
    m_scalargetfunctions[nxString("numlinesofsight")] = [&, this](double* val)
    {
        *val = (double)m_linesofsight.NumRays();
        return true;
    };

    return true;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::MakeVectorGetFunctions		2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::MakeVectorGetFunctions()
{
    m_vectorgetfunctions[nxString("referencepoint")] = [&, this](int index)
    {
        GEODETIC_INSTANT pt = m_engine.ReferencePoint();
        m_getpropertybuffer.resize(4);
        m_getpropertybuffer[0] = pt.latitude;
        m_getpropertybuffer[1] = pt.longitude;
        m_getpropertybuffer[2] = pt.heightm;
        m_getpropertybuffer[3] = pt.mjd;
        return true;
    };

    m_vectorgetfunctions[nxString("wavel")] = [&, this](int index)
    {
        m_getpropertybuffer = m_wavelen;
        return true;
    };

    m_vectorgetfunctions[nxString("observer")] = [&, this](int index)
    {
        bool ok = true;

        ok = ok && (int)m_linesofsight.NumRays() > index;
        ok = ok && index >= 0;
        if (!ok)
        {
            nxLog::Record(NXLOG_WARNING, "ISKEngine_Stub_TIR::GetPropertyArray, The observer array getproperty specifies an index (%d) which is out of the line of sight range [0..%d]", (int)index, (int)m_linesofsight.NumRays());
        }
        /* else */
        const SKTRAN_LineOfSightEntry_V2* entry;
        ok = ok && m_linesofsight.GetRay(index, &entry);
        nxVector pt = entry->Observer();
        m_getpropertybuffer.resize(3);
        m_getpropertybuffer[0] = pt.X();
        m_getpropertybuffer[1] = pt.Y();
        m_getpropertybuffer[2] = pt.Z();
        return ok;
    };

    m_vectorgetfunctions[nxString("look")] = [&, this](int index)
    {
        bool ok = true;

        ok = ok && (int)m_linesofsight.NumRays() > index;
        ok = ok && index >= 0;
        if (!ok)
        {
            nxLog::Record(NXLOG_WARNING, "ISKEngine_Stub_TIR::GetPropertyArray, The observer array getproperty specifies an index (%d) which is out of the line of sight range [0..%d]", (int)index, (int)m_linesofsight.NumRays());
        }
        /* else */
        const SKTRAN_LineOfSightEntry_V2* entry;
        ok = ok && m_linesofsight.GetRay(index, &entry);
        nxVector pt = entry->Look();
        m_getpropertybuffer.resize(3);
        m_getpropertybuffer[0] = pt.X();
        m_getpropertybuffer[1] = pt.Y();
        m_getpropertybuffer[2] = pt.Z();
        return ok;
    };

    m_vectorgetfunctions[nxString("tangentalts")] = [&](int index)
    {
        bool ok = true;
        m_getpropertybuffer.resize(m_linesofsight.NumRays());
        double ground = m_engine.Coordinates()->AltitudeToRadius(0.0);
        for (int rayidx = 0; rayidx < m_linesofsight.NumRays(); rayidx++)
        {
            auto ray = m_engine.LinesOfSight().RayAt(rayidx);
            double Robs, Tobs, Rt;
            ray->CalculateBaseLineTangentPointDetails(0, &Robs, &Tobs, &Rt);

            m_getpropertybuffer[rayidx] = Rt - ground;
        }

        return ok;
    };

    return true;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::MakeStringSetFunctions		2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::MakeStringSetFunctions()
{
    m_stringsetfunctions["geoidmodel"] = [&, this](const char* cstr)
    {
        bool ok;

        ok = CheckModelNotInitialized("geoidmodel");
        if (ok)
        {
            nxString nxstr(cstr);
            nxstr.MakeUpper();

            if (nxstr == "GEOID_SPHERE")
            {
                m_specs.RayTracingSpecs().SetGeoidModel(nxGeodetic::GEOID_SPHERE);
            }
            else if (nxstr == "IAU1976")
            {
                m_specs.RayTracingSpecs().SetGeoidModel(nxGeodetic::IAU1976);
            }
            else if (nxstr == "GRS80")
            {
                m_specs.RayTracingSpecs().SetGeoidModel(nxGeodetic::GRS80);
            }
            else if (nxstr == "MERIT83")
            {
                m_specs.RayTracingSpecs().SetGeoidModel(nxGeodetic::MERIT83);
            }
            else if (nxstr == "WGS84")
            {
                m_specs.RayTracingSpecs().SetGeoidModel(nxGeodetic::WGS84);
            }
            else
            {
                ok = false;
                nxLog::Record(NXLOG_WARNING, "ISKEngine TIR, Invalid option (%s) for property geoidmodel", cstr);
            }
        }

        return ok;
    };

    m_stringsetfunctions["wfspecies"] = [&, this](const char* cstr)
    {
        bool ok;

        ok = CheckModelNotInitialized("wfspecies");

        if (ok)
        {
            // cstr must contain a list of climatologies, separated by spaces or commas
            m_specs.WeightingFunctionSpecs().ClearWFSpecies();

            std::string handle_list(cstr);
            std::istringstream comma_stream(handle_list);
            std::string handle_str;
            while (std::getline(comma_stream, handle_str, ','))
            {
                std::istringstream space_stream(handle_str);
                while (std::getline(space_stream, handle_str, ' '))
                {
                    if (handle_str.size() > 0)
                    {
                        auto* clim = FindGlobalClimatologyHandle(handle_str.c_str());
                        if (*clim == SKCLIMATOLOGY_UNDEFINED)
                        {
                            nxLog::Record(NXLOG_WARNING, "ISKEngine_Stub_TIR::MakeStringSetFunctions, Unknown climatology when configuring weighting functions: %s", handle_str.c_str());
                            return false;
                        }
                        m_specs.WeightingFunctionSpecs().AddWeightingFunctionSpecies(*clim);
                    }
                }
            }
        }

        return ok;
    };

    m_stringsetfunctions["addwfspecies"] = [&, this](const char* cstr)
    {
        bool ok;

        ok = CheckModelNotInitialized("addwfspecies");

        if (ok)
        {
            auto* handle = FindGlobalClimatologyHandle(cstr);
            if (*handle == SKCLIMATOLOGY_UNDEFINED)
            {
                nxLog::Record(NXLOG_WARNING, "ISKEngine_Stub_TIR::MakeStringSetFunctions, Unknown climatology handle [%s] when adding weighting function", cstr);
                return false;
            }
            m_specs.WeightingFunctionSpecs().AddWeightingFunctionSpecies(*handle);
        }

        return ok;
    };

    return true;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::AddWeightingFunctionSpecies	2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::AddWeightingFunctionSpecies(CLIMATOLOGY_HANDLE& species)
{
    m_specs.WeightingFunctionSpecs().AddWeightingFunctionSpecies(species);
    m_modelisconfigured = false;
    return true;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::SetPolarizationMode			2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::SetPolarizationMode(int specifier)
{
    if (specifier > 0)
    {
        nxLog::Record(NXLOG_ERROR, "Sasktran-TIR does not support polarization!", specifier);
        return false;
    }
    return true;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::ParseCommandAndIndex		2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::ParseCommandAndIndex(const nxString& input, nxString& cmd, int& index)
{
    nxStringArray	tokens;
    int				numtoken;
    bool			ok;

    numtoken = nxStrtok(input, &tokens, "([]) ,:;");
    ok = (numtoken == 2);
    if (ok)
    {
        cmd = tokens.GetAt(0);
        index = atoi(tokens.GetAt(1));
        cmd.MakeLower();
    }

    if (!ok)
    {
        cmd = input;
        cmd.MakeLower();
        index = -1;
    }
    return true;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::AddSpecies					2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::AddSpecies(const CLIMATOLOGY_HANDLE& species, ISKClimatology_Stub* climatology, ISKOpticalProperty_Stub* opticalproperty)
{
    skClimatology*				climptr;
    nxUnknown*					optbaseptr;
    skOpticalProperties*		optptr;
    bool						ok;
    bool						iswfspecies = false;

    optbaseptr = (opticalproperty != nullptr) ? opticalproperty->RawObjectPointer() : nullptr;
    climptr = dynamic_cast<skClimatology*>(climatology->RawObjectPointer());
    optptr = (optbaseptr != nullptr) ? dynamic_cast<skOpticalProperties*>(optbaseptr) : nullptr;

    m_cacheisvalid = m_cacheisvalid && m_opticalstate.ContainsSpecies(species);		// invalidate cached absorption if new species is added

    // invalidate cached cross sections if this is not a weighting function species
    for (const CLIMATOLOGY_HANDLE& handle : m_engine.InternalSpecs().WeightingFunctionSpecs().WFSpecies())
    {
        if (handle == species)
        {
            iswfspecies = true;
        }
    }
    m_cacheisvalid = m_cacheisvalid && iswfspecies;

    ok = m_opticalstate.AddSpecies(species, climptr, optptr);

    return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::AddEmission					2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::AddEmission(const EMISSION_HANDLE& species, ISKEmission_Stub* emissionobject)
{
    nxLog::Record(NXLOG_WARNING, "ISKEngine_Stub_TIR::AddEmission, SASKTRAN TIR automatically includes thermal emissions and no other emissions are supported; use engine option 'groundemissivity' to change surface emissivity from the default value of 1.0");
    return true;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::SetAlbedo					2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::SetAlbedo(double albedo)
{
    return true;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::SetBRDF						2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::SetBRDF(ISKBrdf_Stub* brdf)
{
    return true;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::SetWavelengths				2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::SetWavelengths(const double* wavelen, int numwavelen)
{
    if (!m_modelisconfigured)	// can't change wavelengths after configuring model
    {
        m_wavelen.assign(wavelen, wavelen + numwavelen);
        m_radiance.erase();
        m_cacheisvalid = false;
    }
    return true;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::InitializeModel				2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::InitializeModel()
{
    bool ok = true;
    if (!m_modelisconfigured)
    {
        ok = ok && m_engine.ConfigureModel(m_specs, m_wavelen, m_linesofsight, m_numthreads);
        m_modelisconfigured = ok;
        m_cacheisvalid = false;
        if (!ok)
        {
            nxLog::Record(NXLOG_WARNING, "ISKEngine TIR, Error initializing the model. Thats not good");
        }
    }
    return ok;
}

/*-----------------------------------------------------------------------------
 *		ISKEngine_Stub_TIR::CalculateRadiance						2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::CalculateRadiance(const double** radiance, int* numwavelens, int* numlinesofsight)
{
    std::vector<std::vector<double>> r;
    bool ok = true;
    bool ok1;
    size_t numrays = m_linesofsight.NumRays();
    size_t numwave = m_wavelen.size();

    if (!m_modelisconfigured)
    {
        ok = ok && InitializeModel();
    }

    r.resize(numwave);
    for (size_t i = 0; i < numwave; i++)
    {
        r[i].resize(numrays);
    }
    m_radiance.SetSize(numrays, numwave);
    int numwf = (int)m_engine.NumWF();
    int numwfspecies = (int)m_specs.WeightingFunctionSpecs().NumWFSpecies();
    numwfspecies += m_specs.WeightingFunctionSpecs().GetDoTemperatureWF() ? 1 : 0;
    m_wfbuffer.erase();
    if (numwf != 0)
    {
        m_wfbuffer.SetSize(numwf * numwfspecies, numrays, numwave);
    }
    if (m_usecache && !m_cacheisvalid)
    {
        nxLog::Record(NXLOG_INFO, "ISKEngine_Stub_TIR::CalculateRadiance, Option to use cached cross sections was specified but the cache is not valid.\n"
                                  "Cross sections will be recomputed. This could happen if this is the first call to CalculateRadiance, if a species NOT\n"
                                  "included in the weighting functions was altered, or the if the weighting function species themselves were changed.");
    }
    ok1 = m_engine.CalculateRadiance(&r, &m_opticalstate, m_usecache && m_cacheisvalid);
    std::vector<std::vector<std::vector<std::vector<double>>>> wf;
    if ((numwf != 0) && ok1)
    {
        wf = m_engine.GetWF();
    }
    for (size_t waveidx = 0; waveidx < numwave; waveidx++)
    {
        for (size_t rayidx = 0; rayidx < numrays; rayidx++)
        {
            m_radiance.At(rayidx, waveidx) = ok1 ? r.at(waveidx).at(rayidx) : std::numeric_limits<double>::quiet_NaN();
        }
        if ((numwf != 0) && ok1)
        {
            for (size_t rayidx = 0; rayidx < numrays; rayidx++)
            {
                for (size_t speciesidx = 0; speciesidx < numwfspecies; speciesidx++)
                {
                    for (size_t wfidx = 0; wfidx < numwf; wfidx++)
                    {
                        m_wfbuffer.At(speciesidx * numwf + wfidx, rayidx, waveidx) = wf[waveidx][rayidx][speciesidx][wfidx];
                    }
                }
            }
        }
        ok = ok && ok1;
    }
    *radiance = m_radiance.ArrayBasePtr();
    *numwavelens = (int)numwave;
    *numlinesofsight = (int)numrays;
    if (ok) m_cacheisvalid = true;

    return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::CalculateStokesVector		2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::CalculateStokesVector(const ISKStokesVector** radiance, int* numwavelens, int* numlinesofsight)
{
    nxLog::Record(NXLOG_ERROR, "Sasktran-TIR does not support polarized calculations. Use radiances as opposed to stokes vectors.");
    return false;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::GetWeightingFunctions		2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::GetWeightingFunctions(const double** wf, int* numwavel, int* numlos, int* numwf)
{
    int numcalcwf = (int)m_engine.NumWF();
    if (numcalcwf == 0)
    {
        *wf = nullptr;
        *numwavel = 0;
        *numlos = 0;
        *numwf = 0;
    }
    else
    {
        *numwavel = (int)m_wfbuffer.ZSize();
        *numlos = (int)m_wfbuffer.YSize();
        *numwf = (int)m_wfbuffer.XSize();

        *wf = m_wfbuffer.ArrayBasePtr();
    }

    return true;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::SetAtmosphericState			2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::SetAtmosphericState(ISKClimatology_Stub* climatology)
{
    bool ok = true;

    ok = ok && m_opticalstate.SetAtmosphericStateModel(dynamic_cast<skClimatology*>(climatology->RawObjectPointer()));
    m_cacheisvalid = false;		// invalidate cross sections b/c temperature and pressure used to compute them may have changed

    return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::SetPropertyScalar			2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::SetPropertyScalar(const char* propertyname, double value)
{
    nxString str(propertyname);
    str.MakeLower();
    auto funciterator = m_scalarsetfunctions.find(str);
    if (funciterator == std::end(m_scalarsetfunctions))
    {
        return false;
    }
    else
    {
        return funciterator->second(value);
    }
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::SetPropertyArray			2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::SetPropertyArray(const char* propertyname, const double* value, int numpoints)
{
    bool	ok;
    nxString str(propertyname);
    str.MakeLower();

    auto funciterator = m_vectorsetfunctions.find(str);
    if (funciterator == std::end(m_vectorsetfunctions))
    {
        nxLog::Record(NXLOG_WARNING, "ISKEngine_Stub_TIR::Set, this object does not support array property [%s]\n", (const char*)propertyname);
        return false;
    }
    else
    {
        return funciterator->second(value, numpoints);
    }
    return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::SetPropertyObject			2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::SetPropertyObject(const char* propertyname, nxUnknown* object)
{
    nxString str(propertyname);
    str.MakeLower();
    auto funciterator = m_objectsetfunctions.find(str);
    if (funciterator == std::end(m_objectsetfunctions))
    {
        nxLog::Record(NXLOG_WARNING, "ISKEngine_Stub_TIR::Set, this object does not support any object properties including [%s]\n", (const char*)propertyname);
        return false;
    }
    else
    {
        return funciterator->second(object);
    };
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::SetPropertyString			2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::SetPropertyString(const char* propertyname, const char* value)
{
    nxString str(propertyname);
    str.MakeLower();
    auto funciterator = m_stringsetfunctions.find(str);
    if (funciterator == std::end(m_stringsetfunctions))
    {
        nxLog::Record(NXLOG_WARNING, "ISKEngine_Stub_TIR::Set, this object does not support any string properties including [%s]\n", (const char*)propertyname);
        return false;
    }
    else
    {
        return funciterator->second(value);
    }
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::GetPropertyScalar			2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::GetPropertyScalar(const char* propertyname, double* value)
{
    nxString str(propertyname);
    str.MakeLower();
    auto funciterator = m_scalargetfunctions.find(str);
    if (funciterator == std::end(m_scalargetfunctions))
    {
        return false;
    }
    else
    {
        return funciterator->second(value);
    }
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_TIR::Get Property				2018-06-18
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_TIR::GetProperty(const char* propertyname, const double** value, int* numpoints)
{
    nxString name(propertyname);
    nxString cmd;
    int index;
    double	scalarvalue;
    bool ok;

    ok = GetPropertyScalar(propertyname, &scalarvalue);
    if (ok)
    {
        m_getpropertybuffer.resize(1, scalarvalue);
        *numpoints = 0;
        *value = &m_getpropertybuffer[0];
    }
    else
    {

        ok = ParseCommandAndIndex(name, cmd, index);
        auto funcinterator = m_vectorgetfunctions.find(cmd);
        if (funcinterator == std::end(m_vectorgetfunctions))
        {
            nxLog::Record(NXLOG_WARNING, "ISKEngine(TIR), The TIR engine not support property <%s>", (const char*)cmd);
            *numpoints = 0;
            *value = NULL;
            ok = false;
        }
        else
        {
            ok = funcinterator->second(index);
            *numpoints = (int)m_getpropertybuffer.size();
            *value = &m_getpropertybuffer[0];
        }
    }
    return ok;
}
