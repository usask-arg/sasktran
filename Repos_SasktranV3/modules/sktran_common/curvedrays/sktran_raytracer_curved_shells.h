

class SKTRAN_RayTracer_Curved_Shells : public SKTRAN_RayTracer_Base
{
public:
	std::shared_ptr<SKTRAN_GridDefRayTracingShells_V21> m_shells;
	std::unique_ptr< skRTRefractiveIndex_Profile> m_refractiveindex;

	// double IndexOfRefraction(double altitude) const { return 1.0 + (100000.0 - altitude) / (100000.0) * 0.0003; }
	double IndexOfRefraction(double altitude) const { return m_refractiveindex->ExponentialLinearInterp(altitude); }
	double TangentRadius(double firstguessrt) const;
	bool IntegratePath(double rt, double nt, double r1, double r2, double* ds, double* dphi) const;

public:
	SKTRAN_RayTracer_Curved_Shells(std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords);
	virtual ~SKTRAN_RayTracer_Curved_Shells();

	void Initialize(std::shared_ptr<SKTRAN_GridDefRayTracingShells_V21> shells, std::unique_ptr< skRTRefractiveIndex_Profile> refracprofile);
	bool ConfigureOptical(SKTRAN_AtmosphericOpticalState_V21 *opticalstate, double wavelen_nm, GEODETIC_INSTANT referencepoint) override;

	virtual bool TraceRay(SKTRAN_RayOptical_Curved* aray) const;


};