//#include "sktran_hr_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_LinesOfSightTable		2014-10-30*/
/**  Class to store lines of sight and their associated rays for the
 *   engine.  Input lines of sight are translated to the osculating
 *   sphere grid
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_HR_LinesOfSightTable
{
	private:
		std::vector< std::unique_ptr<SKTRAN_RayOptical_Base> >		m_opticalrays;
		SKTRAN_LineOfSightArray_V21									m_observerlinesofsight;

	private:
		bool											ReleaseResources();

	public:
														SKTRAN_HR_LinesOfSightTable		() {};
		virtual										   ~SKTRAN_HR_LinesOfSightTable		();
		bool											SetLinesOfSight					( const SKTRAN_LineOfSightArray_V21& linesofsight, const SKTRAN_CoordinateTransform_V2& coords );
		bool											CreateRays						( const SKTRAN_RayFactory_Base* rayfactory );
		SKTRAN_RayOptical_Base*							RayAt							( size_t idx );
		std::unique_ptr<SKTRAN_RayOptical_Base>&		RayEntryAt						( size_t idx );
		const SKTRAN_RayOptical_Base*					RayAt							( size_t idx ) const;
		SKTRAN_LineOfSightArray_V21*					LinesOfSightArray				() { return &m_observerlinesofsight; }
		size_t											NumRays							() { return m_opticalrays.size(); }
		double											MeanMJD							() const;

};
