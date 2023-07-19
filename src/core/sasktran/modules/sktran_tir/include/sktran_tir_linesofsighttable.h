/**
 * SASKTRAN TIR Lines of Sight Table
 */

/**
 * SKTRAN_TIR_LinesOfSightTable
 * 2019-04-22
 */
class SKTRAN_TIR_LinesOfSightTable
{
private:
	std::vector<std::unique_ptr<SKTRAN_RayOptical_Base>>		m_opticalrays;
	SKTRAN_LineOfSightArray_V21									m_observerlinesofsight;

private:
	bool											ReleaseResources();

public:
	SKTRAN_TIR_LinesOfSightTable() {};
	virtual										   ~SKTRAN_TIR_LinesOfSightTable();
	bool											SetLinesOfSight(const SKTRAN_LineOfSightArray_V21& linesofsight,
																	const SKTRAN_CoordinateTransform_V2& coords);
	bool											CreateRays(const SKTRAN_RayFactory_Base* rayfactory);
	SKTRAN_RayOptical_Base*							RayAt(size_t idx);
	std::unique_ptr<SKTRAN_RayOptical_Base>&		RayEntryAt(size_t idx);
	const SKTRAN_RayOptical_Base*					RayAt(size_t idx) const;
	SKTRAN_LineOfSightArray_V21*					LinesOfSightArray() { return &m_observerlinesofsight; }
	size_t											NumRays() { return m_opticalrays.size(); }
	double											MeanMJD() const;
};
