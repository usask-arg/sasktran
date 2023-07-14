/**
 * SASKTRAN TIR Ray Storage
 *
 * Extends the SASKTRAN raystorage classes to allow saving of parameters used
 * for calculating weighting functions in the TIR engine.
 */

/**
 * SKTRAN_RayStorage_Straight_TIR
 * 2018-09-17
 */
class SKTRAN_RayStorage_Straight_TIR : public SKTRAN_RayStorage_Straight
{
private:

	mutable std::vector<double> m_extinction;

	// weighting functions, maps each retrieval quantity to its weighting function, defined at each quadrature point
	mutable std::map<CLIMATOLOGY_HANDLE, std::vector<double>> m_wf;

public:
	SKTRAN_RayStorage_Straight_TIR(std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords) : SKTRAN_RayStorage_Straight(coords) {}
	
	// SKTRAN_RayStorage_Straight Overrides
	virtual bool Reserve(size_t numquadraturepoints) override;
	virtual bool Resize(size_t numquadraturepoints) override;
	virtual void TruncateToNumElements(size_t numels) override;
	virtual void ClearStorage() override;
	virtual bool PushBack(SKTRAN_Distance r,
						  SKTRAN_Distance distFromTan,
						  SKTRAN_Distance s) override;
	virtual bool Insert(SKTRAN_Distance r,
						SKTRAN_Distance distFromTan,
						SKTRAN_Distance s,
						size_t index) override;
	virtual double CellLength(size_t cellindex) const override;
	
	// TIR Functionality
	virtual double WFAtPoint(const CLIMATOLOGY_HANDLE& species, size_t quadindex) const { return m_wf[species][quadindex]; }
	virtual void SetWFAtPoint(const CLIMATOLOGY_HANDLE& species, size_t quadindex, double value) const { m_wf[species][quadindex] = value; }
	virtual void AddWFSpecies(const CLIMATOLOGY_HANDLE& species) const;
	virtual double ExtinctionAtCellStart(size_t cellindex) const { return m_extinction[cellindex]; }
	virtual void SetExtinction(size_t cellindex, double ext) const { m_extinction[cellindex] = ext; }
	
	std::map<CLIMATOLOGY_HANDLE, std::vector<double>>& WFStorageObject() const { return m_wf; }
};

/**
 * SKTRAN_RayStorage_CurvedPiecewise_TIR
 * 2018-10-10
 */
class SKTRAN_RayStorage_CurvedPiecewise_TIR : public SKTRAN_RayStorage_CurvedPiecewise
{
private:
	mutable std::vector<double> m_extinction;

	mutable std::map<CLIMATOLOGY_HANDLE, std::vector<double>> m_wf;

public:
	// SKTRAN_RayStorage_CurvedPiecewise Overrides
	SKTRAN_RayStorage_CurvedPiecewise_TIR(std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords) : SKTRAN_RayStorage_CurvedPiecewise(coords) {}

	virtual bool Reserve(size_t numquadraturepoints) override;
	virtual bool Resize(size_t numquadraturepoints) override;
	virtual void TruncateToNumElements(size_t numels) override;
	virtual void ClearStorage() override;
	virtual bool PushBack(HELIODETIC_UNITVECTOR* uv,
						  HELIODETIC_POINT* pt,
						  SKTRAN_Distance celllength) override;
	virtual bool Insert(HELIODETIC_UNITVECTOR* uv,
						HELIODETIC_POINT* pt,
						SKTRAN_Distance celllength,
						size_t index) override;

	// TIR Functionality
	virtual double WFAtPoint(const CLIMATOLOGY_HANDLE& species, size_t quadindex) const { return m_wf[species][quadindex]; }
	virtual void SetWFAtPoint(const CLIMATOLOGY_HANDLE& species, size_t quadindex, double value) const { m_wf[species][quadindex] = value; }
	virtual void AddWFSpecies(const CLIMATOLOGY_HANDLE& species) const;
	virtual double ExtinctionAtCellStart(size_t cellindex) const override { return m_extinction[cellindex]; }
	virtual void SetExtinction(size_t cellindex, double ext) const override { m_extinction[cellindex] = ext; }

	std::map<CLIMATOLOGY_HANDLE, std::vector<double>>& WFStorageObject() const { return m_wf; }
};
