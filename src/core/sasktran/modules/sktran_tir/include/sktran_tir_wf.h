/**
 * SASKTRAN TIR Weighting Functions
 * 2018-09-12
 */

/**
 * SKTRAN_TIR_Perturbation_Base
 *
 * Defines perturbations for weighting function calculations. Based on the
 * perturbation schemes used by the HR engine.
 */
class SKTRAN_TIR_Perturbation_Base
{
public:
	SKTRAN_TIR_Perturbation_Base() { };
	virtual ~SKTRAN_TIR_Perturbation_Base() { };

	virtual HELIODETIC_POINT PerturbationLocation(const SKTRAN_CoordinateTransform_V2& coords) const = 0;
	virtual double PerturbationAltitudeWidth() const = 0;
	virtual size_t NumBoundingGeometry() const = 0;
	virtual std::unique_ptr<SKTRAN_GeometryObject> BoundingGeometryObject(std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords,
																		  size_t idx) const = 0;
};

/**
 * SKTRAN_TIR_Perturbation_Absorption_Linear
 *
 * Defines a perturbation centered at an altitude which decreases linearly in
 * radius from the center value to a value of 0.
 */
class SKTRAN_TIR_Perturbation_Absorption_Linear : public SKTRAN_TIR_Perturbation_Base
{
private:
	double m_center;				// center altitude of perturbation
	double m_distancetozeroabove;	// from m_center to altitude where perturbation is zero
	double m_distancetozerobelow;	// typically same value as m_distancetozeroabove
public:
	// Initialize perturbation before calling other methods
	bool Initialize(double center,
					double distancetozerobelow,
					double distancetozeroabove);

	// Get the relative weighting for a perturbation at the altitude given by location.
	virtual bool PerturbationWeight(const HELIODETIC_POINT& location,
									bool* isperturbation,
									double* value) const;
	virtual bool PerturbationWeight(const double altitude,
									bool* isperturbation,
									double* value) const;

	// Get location of the center point of the perturbation
	virtual HELIODETIC_POINT PerturbationLocation(const SKTRAN_CoordinateTransform_V2& coords) const;

	// Bounding geometry objects are spherical shells corresponding to m_center, m_distancetozeroabove, and m_distancetozerobelow
	virtual size_t NumBoundingGeometry() const override { return 3; };	
	virtual std::unique_ptr<SKTRAN_GeometryObject> BoundingGeometryObject(std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords,
																		  size_t idx) const override;

	virtual double PerturbationAltitudeWidth() const { return m_distancetozeroabove + m_distancetozerobelow; }
};
