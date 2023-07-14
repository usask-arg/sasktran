/**
 * SASKTRAN TIR Thread Manager
 */

class SKTRAN_TIR_Thread_Storage
{
private:
	SKTRAN_TIR_LinesOfSightTable m_linesofsighttable;

private:
	bool ReleaseResources();

public:
	SKTRAN_TIR_Thread_Storage();
	~SKTRAN_TIR_Thread_Storage();
	// Won't compile if these aren't defined
	SKTRAN_TIR_Thread_Storage(const SKTRAN_TIR_Thread_Storage& other) {}
	SKTRAN_TIR_Thread_Storage& operator= (SKTRAN_TIR_Thread_Storage& other) { return *this; }
	bool Initialize(const SKTRAN_LineOfSightArray_V21& linesofsight,
					std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords);
	bool Configure(const SKTRAN_RayFactory_Base* rayfactory);
	const size_t NumRays() { return m_linesofsighttable.NumRays(); }
	bool TraceRayAt(size_t losidx);
	SKTRAN_RayOptical_Base* RayAt(size_t losidx) { return m_linesofsighttable.RayAt(losidx); }
};

class SKTRAN_TIR_Thread_Manager
{
private:
	size_t m_numthreads;
	std::vector<SKTRAN_TIR_Thread_Storage> m_threadstore;

private:
	bool ReleaseResources();

public:
	SKTRAN_TIR_Thread_Manager();
	~SKTRAN_TIR_Thread_Manager();

	bool Initialize(size_t numthreads,
					const SKTRAN_LineOfSightArray_V21& linesofsight,
					std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords,
					const SKTRAN_RayFactory_Base* rayfactory);
	bool CalculateRadianceMultiWavel(std::vector<std::vector<SKTRAN_StokesScalar>>* losradiance,
									 std::vector<double>& wavelen,
									 SKTRAN_TIR_AtmosphericOpticalState* opticalstate,
									 std::vector<std::vector<skRTStokesVector>>* losvector);

	std::vector<SKTRAN_TIR_Thread_Storage>& ThreadStorage() { return m_threadstore; }
};
