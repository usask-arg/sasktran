#include "include/sktran_tir_internals.h"

/*-----------------------------------------------------------------------------
 *			skRTRefractiveIndex_Profile::skRTRefractiveIndex_Profile		2012-8-27*/
/** **/
/*---------------------------------------------------------------------------*/

skRTRefractiveIndex_Profile_TIR::skRTRefractiveIndex_Profile_TIR()
{
}


/*-----------------------------------------------------------------------------
 *			skRTRefractiveIndex_Profile::~skRTRefractiveIndex_Profile		2012-8-27*/
/** **/
/*---------------------------------------------------------------------------*/

skRTRefractiveIndex_Profile_TIR::~skRTRefractiveIndex_Profile_TIR()
{
}


/*-----------------------------------------------------------------------------
 *			skRTRefractiveIndex_Profile::CalculateProfileTIR		2017-10-11*/
 /** **/
 /*---------------------------------------------------------------------------*/

bool skRTRefractiveIndex_Profile_TIR::CalculateProfileTIR(SKTRAN_TIR_AtmosphericOpticalState *opticalstate, const SKTRAN_GridDefRayTracingShells_V21 *raytracingspecs, double wavelen_nm, GEODETIC_INSTANT referencepoint)
{
	bool					ok;
	//double					*shellheight;
	skClimatology			*species_temperature_pressure;
	size_t					heightidx;
	size_t					numshells;
	double					wavenum = (1 / wavelen_nm) * 1.0E07;
	size_t					numbad;


	// (re-)calculate the index of refraction profile based on the opticalstate and ray tracing grid

	ok = (referencepoint.heightm >= 0.0) &&
		(referencepoint.latitude >= -90.0 && referencepoint.latitude <= 90.0) &&
		(referencepoint.longitude >= -180.0 && referencepoint.longitude <= 360.0) &&
		(referencepoint.mjd >= 10000.0);	// mjd is invalid for less than 10000 for some reason

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skRTRefractiveIndex_Profile::CalculateProfile, Reference point given not valid.");
	}

	if (ok)
	{
		m_referencepoint = referencepoint;

		ok = ok && opticalstate->SetTimeAndLocation(referencepoint, true);

		numshells = raytracingspecs->NumShells();

		if (numshells == 0)
		{
			nxLog::Record(NXLOG_WARNING, "skRTRefractiveIndex_Profile::CalculateProfile, Please set the shell altitudes for the refractivity grid.");
		}
		if (ok &&  numshells != 0)
		{
			m_heights.resize(numshells);
			m_temperature.resize(numshells);
			m_temperature.resize(numshells);
			m_refractiveindex.resize(numshells);

			//shellheight = raytracingspecs->begin();
			heightidx = 0;

			// may need to add in check that h <= 120 km, since pressure/temperature profiles only valid till 120 km
			/*while ( shellheight != raytracingspecs->end() )
			{
				m_heights.at(heightidx) = *shellheight;	// midpoint of cell
				++heightidx;
				++shellheight;
			}*/

			for (size_t heightidx = 0; heightidx < raytracingspecs->NumGridPoints(); heightidx++)
				m_heights.at(heightidx) = raytracingspecs->At(heightidx);

			opticalstate->GetAtmosphericStateModel(&species_temperature_pressure);

			if (ok)
			{
				m_pressure.resize(m_heights.size());
				m_temperature.resize(m_heights.size());

				species_temperature_pressure->GetHeightProfile(SKCLIMATOLOGY_PRESSURE_PA, m_referencepoint, &m_heights.front(), (int)m_heights.size(), &m_pressure.front(), true, &numbad);
				species_temperature_pressure->GetHeightProfile(SKCLIMATOLOGY_TEMPERATURE_K, m_referencepoint, &m_heights.front(), (int)m_heights.size(), &m_temperature.front(), true, &numbad);

				UpdateRefractiveIndexTIR(opticalstate, wavenum);
			}
		}
	}

	InitializeCubicSplineInterpolation();	// once the index of refraction profile is found

	return ok;
}


void skRTRefractiveIndex_Profile_TIR::UpdateRefractiveIndexTIR(SKTRAN_TIR_AtmosphericOpticalState *opticalstate, double wavenum)
{
	bool					ok_h2o, ok_co2;
	skClimatology			*species_h2o, *species_co2, *species_air;	// search for other species, for now, only other one is water vapour and co2
	std::vector<double>		wvp;
	std::vector<double>		co2p;
	std::vector<double>		airp;
	size_t					numbad;
	size_t					numshells;
	double					Rv = 461.495; // specific gas constant for water vapour (J/(kg K))
	std::complex<double>	n;

	numshells = m_heights.size();

	const CLIMATOLOGY_HANDLE waterName = SKCLIMATOLOGY_H2O_CM3;
	ok_h2o = opticalstate->GetSpeciesClimatology(waterName, &species_h2o);	// get humidity profile, if it's there

	const CLIMATOLOGY_HANDLE co2Name = SKCLIMATOLOGY_CO2_CM3;
	ok_co2 = opticalstate->GetSpeciesClimatology(co2Name, &species_co2);
	const CLIMATOLOGY_HANDLE airName = SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3;
	ok_co2 = ok_co2 && opticalstate->GetSpeciesClimatology(airName, &species_air);  // need air density to determine CO2 VMR

	if (ok_h2o && ok_co2)
	{
		wvp.resize(m_heights.size());
		species_h2o->GetHeightProfile(SKCLIMATOLOGY_H2O_CM3, m_referencepoint, &m_heights.front(), (int)m_heights.size(), &wvp.front(), true, &numbad);
		co2p.resize(m_heights.size());
		species_co2->GetHeightProfile(SKCLIMATOLOGY_CO2_CM3, m_referencepoint, &m_heights.front(), (int)m_heights.size(), &co2p.front(), true, &numbad);
		airp.resize(m_heights.size());
		species_air->GetHeightProfile(SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3, m_referencepoint, &m_heights.front(), (int)m_heights.size(), &airp.front(), true, &numbad);

		for (size_t i = 0; i < numshells; i++)
		{
			// wvp is not yet correct, it is still a density:
			wvp[i] *= 1000000; // in per cm^3, change to per m^3
			wvp[i] *= Rv * m_temperature[i];
			wvp[i] *= 2.99151e-26; //kg/molecule
			// wvp is now correct for this height
			// convert co2 from /cm3 to ppm
			co2p[i] = co2p[i] / airp[i] * 1.0e6;

			m_calculator->Set_WaterVapourPressure(wvp[i]);
			m_calculator->Set_TotalPressure(m_pressure[i]);
			m_calculator->Set_Temperature(m_temperature[i]);
			m_calculator->Set_CO2ppm(co2p[i]);

			n = m_calculator->RefractiveIndex(wavenum);
			m_refractiveindex[i] = n.real();
		}
	}
	else if (ok_h2o && !ok_co2)
	{
		wvp.resize(m_heights.size());
		species_h2o->GetHeightProfile(SKCLIMATOLOGY_H2O_CM3, m_referencepoint, &m_heights.front(), (int)m_heights.size(), &wvp.front(), true, &numbad);

		for (size_t i = 0; i < numshells; i++)
		{
			// wvp is not yet correct, it is still a density:
			wvp[i] *= 1000000; // in per cm^3, change to per m^3
			wvp[i] *= Rv * m_temperature[i];
			wvp[i] *= 2.99151e-26; //kg/molecule
			// wvp is now correct for this height
			m_calculator->Set_WaterVapourPressure(wvp[i]);
			m_calculator->Set_TotalPressure(m_pressure[i]);
			m_calculator->Set_Temperature(m_temperature[i]);

			n = m_calculator->RefractiveIndex(wavenum);
			m_refractiveindex[i] = n.real();
		}
	}
	else if (ok_co2 && !ok_h2o)
	{
		co2p.resize(m_heights.size());
		species_co2->GetHeightProfile(SKCLIMATOLOGY_CO2_CM3, m_referencepoint, &m_heights.front(), (int)m_heights.size(), &co2p.front(), true, &numbad);
		airp.resize(m_heights.size());
		species_air->GetHeightProfile(SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3, m_referencepoint, &m_heights.front(), (int)m_heights.size(), &airp.front(), true, &numbad);

		for (size_t i = 0; i < numshells; i++)
		{
			// convert co2 from /cm3 to ppm
			co2p[i] = co2p[i] / airp[i] * 1.0e6;

			m_calculator->Set_TotalPressure(m_pressure[i]);
			m_calculator->Set_Temperature(m_temperature[i]);
			m_calculator->Set_CO2ppm(co2p[i]);

			n = m_calculator->RefractiveIndex(wavenum);
			m_refractiveindex[i] = n.real();
		}
	}
	else
	{
		for (size_t i = 0; i < numshells; i++)
		{
			m_calculator->Set_TotalPressure(m_pressure[i]);
			m_calculator->Set_Temperature(m_temperature[i]);

			n = m_calculator->RefractiveIndex(wavenum);
			m_refractiveindex[i] = n.real();
		}
	}

}

