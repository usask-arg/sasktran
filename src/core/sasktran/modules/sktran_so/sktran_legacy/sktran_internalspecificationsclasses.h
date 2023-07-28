/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_RayTracing_V21_General	2010-4-28*/
/**   @ingroup specs
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_SpecsInternal_RayTracing_V21 : public nxUnknown
{
	private:
		std::shared_ptr<SKTRAN_GridDefRayTracingShells_V21>	m_raytracingshells;

	public:
															SKTRAN_SpecsInternal_RayTracing_V21();
		virtual											   ~SKTRAN_SpecsInternal_RayTracing_V21();
		bool												ConfigureRayTracingShellAlts	( const double*  alts_meters, size_t npts );

	public:
		virtual std::shared_ptr<SKTRAN_GridDefRayTracingShells_V21>	RayTracingShells () const { return m_raytracingshells;}
		virtual       size_t										MaxShellsAlongRay() const;

};



/*-----------------------------------------------------------------------------
 *					class SKTRAN_SpecsInternal_OpticalPropertiesGrid_V21		2010-4-14*/
/**   @ingroup specs
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_SpecsInternal_OpticalPropertiesGrid_V21: public nxUnknown
{
	public:
																	SKTRAN_SpecsInternal_OpticalPropertiesGrid_V21() {}
		virtual													   ~SKTRAN_SpecsInternal_OpticalPropertiesGrid_V21() {}

		virtual const SKTRAN_GridDefOpticalPropertiesRadii_V21*		OpticalPropertiesGrid						  ()  const  = 0;
		virtual bool										CreateEmptyOpticalPropertiesTable			  ( SKTRAN_TableOpticalProperties_V21**	  opticalpropertiestable ) const = 0; 
};
