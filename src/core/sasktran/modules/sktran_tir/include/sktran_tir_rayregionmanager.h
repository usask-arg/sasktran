/**
 * SASKTRAN TIR Ray Tracing Region Manager
 */

/**
 * SKTRAN_TIR_RayTracingRegionManager
 * 2018-09-13
 *
 * Similar to the SKTRAN_HR_RayTracingRegionManager, this class allows the
 * ray manager to determine and save the locations where the line of sight
 * enters and leaves the atmosphere.
 */
class SKTRAN_TIR_RayTracingRegionManager : public SKTRAN_RayTracingRegionManager
{
private:
	nxVector m_inreferencepoint;
	nxVector m_outreferencepoint;

	double m_fixedearthradius;
	double m_minalt;

private:
	virtual bool UpdateBoundingReferences(const SKTRAN_LineOfSightArray_V21& linesofsight);
	virtual bool GetRayEndpointsObserverOutside(const nxVector& observer,
												const nxVector& look,
												nxVector* startpt,
												nxVector* endpt);
public:
	SKTRAN_TIR_RayTracingRegionManager() : SKTRAN_RayTracingRegionManager() { m_fixedearthradius = std::numeric_limits<double>::quiet_NaN(); m_minalt = 0.0; }
	~SKTRAN_TIR_RayTracingRegionManager() {}
	bool GetBoundingReferences(nxVector& in,
							   nxVector& out) const;
	virtual bool UpdateUndefinedParametersFromLinesOfSight(const SKTRAN_LineOfSightArray_V21& linesofsight);
	bool SetEarthRadius(double earthradius_meters);
	bool MakeCoordinateSystem(std::shared_ptr<const SKTRAN_CoordinateTransform_V2>* usercoordinates,
							  double groundaltitude,
							  double toa_altitude,
							  nxGeodetic::GEOID_MODEL geoidmodel,
							  bool userdefinedgeoidmodel) const;
	bool CheckParameters() const;
	void Clear();
	bool SetLowerBoundAltitude(double lowerheight_meters);
};
