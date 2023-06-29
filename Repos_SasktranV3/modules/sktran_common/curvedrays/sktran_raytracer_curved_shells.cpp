#include "../sktran_common.h"
#include <boost/math/quadrature/trapezoidal.hpp>

SKTRAN_RayTracer_Curved_Shells::SKTRAN_RayTracer_Curved_Shells(std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords) : SKTRAN_RayTracer_Base(coords)
{
}

SKTRAN_RayTracer_Curved_Shells::~SKTRAN_RayTracer_Curved_Shells()
{
}

void SKTRAN_RayTracer_Curved_Shells::Initialize(std::shared_ptr<SKTRAN_GridDefRayTracingShells_V21> shells, std::unique_ptr<skRTRefractiveIndex_Profile> refracprofile)
{
	m_shells = shells;
	m_refractiveindex = std::move(refracprofile);
}

bool SKTRAN_RayTracer_Curved_Shells::ConfigureOptical(SKTRAN_AtmosphericOpticalState_V21 * opticalstate,  double wavelen_nm, GEODETIC_INSTANT referencepoint)
{
	return m_refractiveindex->CalculateProfile(opticalstate, m_shells.get(), wavelen_nm, referencepoint);
}

double SKTRAN_RayTracer_Curved_Shells::TangentRadius(double firstguessrt) const
{
	const size_t maxiter = 500;			// Takes no time, may as well do a bunch
	const double tolerance = 1e-6;		// 1 micron

	size_t currentiter = 1;

	double currentrt = firstguessrt;
	double nextrt = 0.0;

	while (currentiter < maxiter)
	{
		double n = IndexOfRefraction(CoordsObject()->RadiusToAltitude(currentrt));

		nextrt = firstguessrt / n;

		if (fabs(nextrt - currentrt) < tolerance)
		{
			break;
		}
		currentiter++;
		if (currentiter != maxiter)
		{
			currentrt = nextrt;
		}
	}

	if (currentiter == maxiter)
	{
		// Usually not a problem, but this can happen at low tangent altitudes
		if (fabs(nextrt - currentrt) < 1)
		{
			// Somewhat close to convergence, just use the average of the last two iterations
			currentrt = (currentrt + nextrt) / 2.0;
		}
		else
		{
			// Actually failed to converge, issue a warning
			nxLog::Record(NXLOG_WARNING, "SKTRAN_RayTracer_Curved TangentRadius failed to converge within the maximum number of iterations");
		}
	}

	return currentrt;
}

bool SKTRAN_RayTracer_Curved_Shells::IntegratePath(double rt, double nt, double r1, double r2, double* ds, double* dphi) const
{
	// If r1 and r2 are almost the same then we don't need to bother doing the curved calculation
	if(std::abs(r1-r2) < 0.1)
	{
		*ds = (sqrt(r2*r2 - rt * rt) - sqrt(r1*r1 - rt * rt));
		*dphi = 0;
		return true;
	}

	double integral, deflection;
	// Eq 21 from Thompson 1982, Ray tracing in a refracting spherically symmetric atmosphere
	// Integral to find the path length of the cell
	auto integrand = [&, this](double x) {
		double r = x*x + rt;
		double n = this->IndexOfRefraction(CoordsObject()->RadiusToAltitude(r));

		double sqf = sqrt(1.0 + (n - nt) / n * rt / (x*x)) * sqrt(x*x + (n + nt) / n * rt);
		double F = sqrt(x*x + 2.0*rt) + sqf;
		double G = sqrt(x*x + 2.0*rt)*sqf;

		return 2.0*rt*rt*(nt + n) / n * (nt - n) / n * (x*x + rt) / (x*x*F*G);

	};

	try
	{
		integral = boost::math::quadrature::trapezoidal(integrand, sqrt(r1 - rt), sqrt(r2 - rt), 1e-6);
	}
	catch(const std::exception& e)
	{
		nxLog::Record(NXLOG_WARNING, "Failed to calculate path length");
		nxLog::Record(NXLOG_WARNING, e.what());
		nxLog::Record(NXLOG_WARNING, "rt: %f, nt: %f, r1: %f, r2: %f", rt, nt, r1, r2);
	}
	

	if (integral != integral)
	{
		// Error
		nxLog::Record(NXLOG_ERROR, "Error integrating curved ray path");
	}

	double P = integral + (sqrt(r2*r2 - rt * rt) - sqrt(r1*r1 - rt * rt)); // Cell path length

	// Also need to find the deflection angle, use a modified version of eq 12
	auto integrand_angle = [&, this](double x) {
		double r = x * x + rt;
		double n = this->IndexOfRefraction(CoordsObject()->RadiusToAltitude(r));

		double t1 = sqrt((n - nt) / n * r / (x*x) + nt / n);
		double t2 = sqrt(r + rt*nt/n);

		return 2.0*nt*rt / r * (1.0 / (n*t1*t2));
	};

	try
	{
		deflection = boost::math::quadrature::trapezoidal(integrand_angle, sqrt(r1 - rt), sqrt(r2 - rt), 1e-10);

	}
	catch(const std::exception& e)
	{
		nxLog::Record(NXLOG_WARNING, "Failed to calculate deflection");
		nxLog::Record(NXLOG_WARNING, e.what());
		nxLog::Record(NXLOG_WARNING, "rt: %f, nt: %f, r1: %f, r2: %f", rt, nt, r1, r2);
	}
	

	// output
	*ds = P;
	*dphi = deflection;

	return true;
}

bool SKTRAN_RayTracer_Curved_Shells::TraceRay(SKTRAN_RayOptical_Curved * aray) const
{
	// A little messy since we handle all ray tracing cases in the same function
	// Set the default case to observer outside the atmosphere, looking through, not hitting the ground
	bool groundishit = false;
	bool lookingdown = true;
	bool observeroutside = true;

	double Robs, Tobsbase, Rtbase;

	aray->CalculateBaseLineTangentPointDetails(0, &Robs, &Tobsbase, &Rtbase);
	aray->StorageVar()->ClearStorage();

	double tangentradius = TangentRadius(Rtbase);
	double tangentaltitude = CoordsObject()->RadiusToAltitude(tangentradius);

	// Check if the observer is inside the atmosphere
	if (Robs < CoordsObject()->AltitudeToRadius(m_shells->At(m_shells->NumCells())))
	{
		observeroutside = false;
		
		// Check to see if we are looking up
		lookingdown = (aray->GetObserver().UnitVector() & aray->LookVector()) < 0.0;
	}


	// Check to see if the refracted tangent radius is below the surface
	if (tangentaltitude < 0.0)
	{
		// Only true if we are looking up, modify this later
		groundishit = true;
	}

	double nt = IndexOfRefraction(tangentaltitude);

	// Will have an intersection with every shell above the tangent radius
	SKTRAN_GridIndex abovetangentindex;
	m_shells->IndexOfPointEqualOrAbove(tangentaltitude, &abovetangentindex);

	if (fabs(CoordsObject()->AltitudeToRadius(m_shells->At(abovetangentindex)) - tangentradius) < 1e-3)
	{
		// We will get errors trying to integrate over a small domain, so
		// Ignore the last shell and just use the tangent radius
		abovetangentindex++;
	}
	
	double ds, dphi;
	std::vector<double> pathlengths;
	std::vector<double> angles;
	std::vector<double> lowerradius;
	std::vector<double> upperradius;

	SKTRAN_GridIndex upperintegrationindex;
	if (observeroutside)
	{
		upperintegrationindex = m_shells->NumCells();
	}

	// Integrate through each shell to get the horizontal angle and curved path length
	// Some of these might not be used depending on the observers geometry
	for (SKTRAN_GridIndex idx = m_shells->NumCells(); idx > abovetangentindex; idx--)
	{
		IntegratePath(tangentradius, nt, CoordsObject()->AltitudeToRadius(m_shells->At(idx-1)), CoordsObject()->AltitudeToRadius(m_shells->At(idx)), &ds, &dphi);
		angles.push_back(dphi);
		pathlengths.push_back(ds);
		lowerradius.push_back(CoordsObject()->AltitudeToRadius(m_shells->At(idx - 1)));
		upperradius.push_back(CoordsObject()->AltitudeToRadius(m_shells->At(idx)));
	}

	if (!groundishit && angles.size() > 0)
	{
		// Have to do an additional integral for the tangent shell
		// Integrate the final tangent shell, Shift the domain of integration by 1 mm to avoid the singularity
		IntegratePath(tangentradius, nt, tangentradius + 1e-3, CoordsObject()->AltitudeToRadius(m_shells->At(abovetangentindex)) + 1e-3, &ds, &dphi);
		angles.push_back(dphi);
		pathlengths.push_back(ds);
		lowerradius.push_back(tangentradius);
		upperradius.push_back(CoordsObject()->AltitudeToRadius(m_shells->At(abovetangentindex)));
	}
	
	// Now calculate the three dimensional locations of all the intersections

	HELIODETIC_VECTOR start;
	HELIODETIC_POINT startpt;
	HELIODETIC_UNITVECTOR look = aray->LookVector();
	HELIODETIC_VECTOR observer = aray->GetObserver();
	// If we are outside the atmosphere we have to first calculate the geometric distance to the top of the atmosphere
	if (observeroutside)
	{
		// Assume no refraction to the top shell, calculate the length geometrically
		double topradius = CoordsObject()->AltitudeToRadius(m_shells->At(m_shells->NumCells())); // Radius of the top shell
		double disttotangent = sqrt(topradius*topradius - Rtbase * Rtbase);                          // Geometric distance of the top shell to the geometric TP
		double disttotop = Tobsbase - disttotangent;                                                 // Distance to the top shell from observer

		start = observer + HELIODETIC_VECTOR(look, disttotop);
		startpt.FromVector(start, CoordsObject().get());
	}
	else
	{
		// We start at the observer location
		start = observer;
		startpt.FromVector(start, CoordsObject().get());
	}
	HELIODETIC_POINT nextpt;
	HELIODETIC_UNITVECTOR uv;

	// Construct the x, y basis for the observer/look vector plane
	nxVector x_basis(start.X(), start.Y(), start.Z());
	x_basis /= x_basis.Magnitude();

	nxVector lookv(look.X(), look.Y(), look.Z());

	// If we are looking straight up then there is a problem, but if we are
	// looking straight up then the ray doesn't curve at all anyways so we can
	// set the y_basis to just be 0
	nxVector y_basis;

	if (std::abs(x_basis.AngleTo(lookv)) < 0.00001 || std::abs(x_basis.AngleTo(-1.0*lookv)) < 0.00001)
	{
		y_basis = nxVector(0, 0, 0);
	}
	else
	{
		nxVector normal = x_basis.Cross(lookv);
		normal /= normal.Magnitude();

		y_basis = normal.Cross(x_basis);
	}

	double total_phi = 0;

	// Put the angles, pathlengths in the correct order for each ray tracing case
	std::vector<double> traceangles;
	std::vector<double> tracepathlengths;
	std::vector<double> tracenextradius;

	if (observeroutside)
	{
		if (groundishit)
		{
			// Observer outside, hitting ground

			traceangles.reserve(angles.size());
			tracepathlengths.reserve(pathlengths.size());
			// We hit all of the shells once
			for (int idx = 0; idx < angles.size(); idx++)
			{
				traceangles.push_back(angles[idx]);
				tracepathlengths.push_back(pathlengths[idx]);
				tracenextradius.push_back(lowerradius[idx]);
			}
		}
		else
		{
			// Observer outside, standard limb view

			// We hit all of the shells twice
			traceangles.reserve(angles.size() * 2);
			tracepathlengths.reserve(pathlengths.size() * 2);
			for (int idx = 0; idx < angles.size(); idx++)
			{
				traceangles.push_back(angles[idx]);
				tracepathlengths.push_back(pathlengths[idx]);
				tracenextradius.push_back(lowerradius[idx]);
			}
			for (int idx = (int)(angles.size()) - 1; idx >= 0; idx--)
			{
				traceangles.push_back(angles[idx]);
				tracepathlengths.push_back(pathlengths[idx]);
				tracenextradius.push_back(upperradius[idx]);
			}
		}
	}
	else
	{
		SKTRAN_GridIndex aboveobserverindex, belowobserverindex;
		m_shells->IndexOfPointEqualOrAbove(CoordsObject()->RadiusToAltitude(Robs), &aboveobserverindex);
		m_shells->IndexOfPointBelowOrEqual(CoordsObject()->RadiusToAltitude(Robs), &belowobserverindex);

		if (lookingdown)
		{
			// We have to integrate an extra shell from the observer to the lower shell if the lower shell
			// is not below the tangent shell
			if (CoordsObject()->AltitudeToRadius(m_shells->At(belowobserverindex)) > tangentradius)
			{
				IntegratePath(tangentradius, nt, CoordsObject()->AltitudeToRadius(m_shells->At(belowobserverindex)), Robs, &ds, &dphi);

				traceangles.push_back(dphi);
				tracepathlengths.push_back(ds);
				tracenextradius.push_back(CoordsObject()->AltitudeToRadius(m_shells->At(belowobserverindex)));
			}

			if (groundishit)
			{
				// Observer inside, looking down hitting ground

				// Just need the observer shell and all below it
				for (int idx = (int)(angles.size() - belowobserverindex); idx < angles.size(); idx++)
				{
					traceangles.push_back(angles[idx]);
					tracepathlengths.push_back(pathlengths[idx]);
					tracenextradius.push_back(lowerradius[idx]);
				}
			}
			else
			{
				// Observer inside, looking down limb view

				int index = (int)std::distance(std::begin(lowerradius), std::lower_bound(std::begin(lowerradius), std::end(lowerradius), Robs, [](const double& a, const double& b) { return b < a; }));

				// Observer shell and all below it
				for (int idx = index + 1; idx < angles.size(); idx++)
				{
					traceangles.push_back(angles[idx]);
					tracepathlengths.push_back(pathlengths[idx]);
					tracenextradius.push_back(lowerradius[idx]);
				}
				// And then all of them on the way out
				for (int idx = (int)(angles.size()) - 1; idx >= 0; idx--)
				{
					traceangles.push_back(angles[idx]);
					tracepathlengths.push_back(pathlengths[idx]);
					tracenextradius.push_back(upperradius[idx]);
				}
			}
		}
		else
		{
			// Need to integrate an extra shell, the observer to the shell just above
			IntegratePath(tangentradius, nt, Robs, CoordsObject()->AltitudeToRadius(m_shells->At(aboveobserverindex)), &ds, &dphi);
			// Observer inside, looking up
			groundishit = false;

			traceangles.push_back(dphi);
			tracepathlengths.push_back(ds);
			tracenextradius.push_back(CoordsObject()->AltitudeToRadius(m_shells->At(aboveobserverindex)));

			int index = (int)std::distance(std::begin(upperradius), std::lower_bound(std::begin(upperradius), std::end(upperradius), Robs, [](const double& a, const double& b) { return b < a; }));

			// -1 since one shell is rt
			for (int idx = index - 2; idx >= 0; idx--)
			{
				traceangles.push_back(angles[idx]);
				tracepathlengths.push_back(pathlengths[idx]);
				tracenextradius.push_back(upperradius[idx]);
			}
		}
	}


	// Calculate the intersection for every shell
	for (size_t idx = 0; idx < traceangles.size(); idx++)
	{
		double phi = traceangles[idx];
		double ds = tracepathlengths[idx];
		double nextradius = tracenextradius[idx];

		if (ds < 1e-6)
		{
			// Small shell that we can just ignore, sometimes we get rounding errors
			// that cause ds to be negative even
			continue;
		}

		total_phi += phi;

		nxVector nextv =  nextradius*(cos(total_phi) * x_basis + sin(total_phi) * y_basis);

		HELIODETIC_VECTOR nexthv;
		nexthv.SetCoords(nextv.X(), nextv.Y(), nextv.Z());

		nextpt.FromVector(nexthv, CoordsObject().get());

		nxVector avglook(nextpt.Vector().X() - startpt.Vector().X(),
			nextpt.Vector().Y() - startpt.Vector().Y(), nextpt.Vector().Z() - startpt.Vector().Z());
		avglook /= avglook.Magnitude();
		uv.SetCoords(avglook.X(), avglook.Y(), avglook.Z());

		aray->StorageVar()->PushBack(&uv, &startpt, ds);
		
		startpt = nextpt;

	}
	// And add the last point
	if (traceangles.size() > 0)
	{
		aray->StorageVar()->PushBack(&uv, &nextpt, -9999999); // ds shouldn't be used here so enter a large negative number to detect errors
	}

	aray->StorageVar()->SetGroundIsHit(groundishit);

	return true;
}
