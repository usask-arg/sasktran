#if !defined(SASKTRANV21_LEGACY_H_INCLUDED)
#define SASKTRANV21_LEGACY_H_INCLUDED


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_OpticalPropertiesGrid_1D_Height		2010-4-29*/
/** **/
/*---------------------------------------------------------------------------*/

class SKTRAN_SpecsUser_OpticalPropertiesGrid_1D_Height : public SKTRAN_SpecsUser_OpticalPropertiesGrid_V21
{
	private:
		std::vector<double>								m_heights;

	public:
														SKTRAN_SpecsUser_OpticalPropertiesGrid_1D_Height(){};
		virtual										   ~SKTRAN_SpecsUser_OpticalPropertiesGrid_1D_Height(){};
		bool											ConfigureOpticalPropertyShells	( const double* alts, size_t npts );

	public:
		virtual bool									CreateInternalSpecs( const SKTRAN_SpecsInternal_OpticalPropertiesGrid_V21** internalspecs ) const;
};

class SKTRAN_TableRayLOSFactory_Legacy;

#endif
