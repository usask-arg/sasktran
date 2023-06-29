//#include "sktran_hr_internals.h"



/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Thread_Storage		2014-10-30*/
/** Class which every thread owns.  Currently every thread only need's its own
 *  ray
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_HR_Thread_Storage
{
	private:
		std::unique_ptr<SKTRAN_RayOptical_Base>			m_rayopt;

	private:
		bool											ReleaseResources();

	public:
														SKTRAN_HR_Thread_Storage		();
													   ~SKTRAN_HR_Thread_Storage		();
													   SKTRAN_HR_Thread_Storage		( const SKTRAN_HR_Thread_Storage& other )		{	if (other.m_rayopt.get() != nullptr)
																																		{
																																			nxLog::Record(NXLOG_ERROR,"SKTRAN_HR_Thread_Storage, Copy Constructor. Cannot copy existing valid thread storage. Code needs reworking");
																																			throw("SKTRAN_HR_Thread_Storage, Copy Constructor Error");
																																		} 
																																		//m_rayopt = std::move(other.m_rayopt);
																																	}
		SKTRAN_HR_Thread_Storage&						operator=						( SKTRAN_HR_Thread_Storage& other )		{ m_rayopt = std::move(other.m_rayopt); return *this;}
		bool											Initialize						( const SKTRAN_RayFactory_Base*   rayfactory);
		SKTRAN_RayOptical_Base*							Ray								(); 
};

//class SKTRAN_HR_Diffuse_Table;
class SKTRAN_HR_Diffuse_Table_CPU;




/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Thread_Manager		2014-10-30*/
/** Workhorse class for the multithreaded portions of the HR engine.
 *  The main operations are
 *    CreateDiffuseFirstOrder - Calculates the first order incoming radiance
 *       to all of the diffuse points
 *    ScatterCPU - Scatters the incoming radiance at every diffuse points
 *       to get the next order outgoing. #order is the order of diffuse scatter
 *       ie #order=1 for the first incoming-to-outgoing scatter, which actually
 *       actually the observer's s
 *    NextOrderCPU - Calculates the next order incoming radiance at every
 *       diffuse point
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_HR_Thread_Manager
{
	private:
		size_t														m_numthreads;
		std::vector<SKTRAN_HR_Thread_Storage >						m_threadstore;
		std::shared_ptr< const SKTRAN_CoordinateTransform_V2>		m_coords;
		const SKTRAN_OpticalPropertiesIntegrator_Base*				m_optint;
		const SKTRAN_SourceTermIntegrator_Base*						m_srcint;
		const SKTRAN_TableOpticalProperties_Base*					m_opttable;
//		const SKTRAN_RayTracer_Base*								m_raytracer;
		std::shared_ptr<const SKTRAN_RayFactory_Base>				m_rayfactory;
		const SKTRAN_SolarTransmission_Base*						m_solartrans;

	private:
		bool										ReleaseResources();

	public:
													SKTRAN_HR_Thread_Manager();
		virtual									   ~SKTRAN_HR_Thread_Manager();

		bool										SetNumThreads				( size_t numthreads );
		bool										CreateDiffuseFirstOrder		( SKTRAN_HR_Diffuse_Table_CPU* diffusetable, double wlen );
		bool										CreateScatteringMatrix		( SKTRAN_HR_Diffuse_Table_CPU* diffusetable );
		//bool										CreateNextOrderMatrix		( SKTRAN_HR_Diffuse_Table& diffusetable );
		bool										ComputeFieldCPU				( SKTRAN_HR_Diffuse_Table_CPU* diffusetable, size_t numorders, double wlen );
		bool										ScatterCPU                  ( SKTRAN_HR_Diffuse_Table_CPU* diffusetable );
		bool										NextOrderCPU				( SKTRAN_HR_Diffuse_Table_CPU* diffusetable );

		virtual bool								Initialize					( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>&	coords,
																				  const SKTRAN_OpticalPropertiesIntegrator_Base&		optintegrator,
																				  const SKTRAN_SourceTermIntegrator_Base&				srcintegrator,
																				  std::shared_ptr< const SKTRAN_RayFactory_Base>		rayfactory,
																				  const SKTRAN_SolarTransmission_Base&					solartrans,
																				  const SKTRAN_TableOpticalProperties_Base&				opttable);

};

