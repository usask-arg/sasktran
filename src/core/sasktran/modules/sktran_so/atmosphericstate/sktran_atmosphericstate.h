



/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_V21		2007-11-20*/
/** A class used to cache the optical properties of the atmosphere for
 *	a range of altitudes and scattering angles.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_TableOpticalProperties_V21 : public SKTRAN_TableOpticalProperties_Base
{
	public:
		 											SKTRAN_TableOpticalProperties_V21	() {}
		virtual 								   ~SKTRAN_TableOpticalProperties_V21	() {}

	//	virtual bool								IsOptionTrue						( SKTRAN_TableOpticalProperties_V21::OPTIONSENUM options) const = 0;
	//	virtual double								TotalExtinctionPerCM				( const HELIODETIC_POINT& point ) const = 0;
	//	virtual double								ScatteringExtinctionPerCM			( const HELIODETIC_POINT& point ) const = 0;
	//	virtual bool								GetScatteringCoefficientCM2			( const HELIODETIC_POINT& point, double cosangle,  SKTRAN_PhaseMatrixScalar* phasematrix ) const = 0;
	//	virtual bool								GetAlbedo							( const HELIODETIC_POINT& point, double* albedo ) const = 0;
	//	virtual bool								GetEffectiveExtinctionPerCMWithHeight( const HELIODETIC_POINT& point, double h0, double r0, double h1, double r1, double* sigmak, double* sigmaf ) const = 0;
	//	virtual bool								ConfigureOptical					( double wavelen, SKTRAN_AtmosphericOpticalState_V21& opticalstate ) = 0;

		virtual bool								ConfigureGeometry					( const SKTRAN_SpecsInternal_Base* specs       ) = 0;
		virtual bool								GetScatteringCoefficientCM2ForJindex( const SKTRANSO_JIndex* aboundingpoints, size_t anumpoints,  double acosangle, SKTRAN_PhaseMatrixScalar* acoefficient ) const = 0;
		virtual bool								GetBoundingSpatialPoints			( const HELIODETIC_POINT& location, SKTRANSO_JIndex* interppoints, size_t* numpoints ) const = 0;
		virtual bool								GetBoundingScatteringPoints			( double cosscatteringangle, SKTRANSO_JIndex* interppoints, size_t* numpoints ) const = 0;
		virtual bool								GetEffectiveExtinctionPerCMWithHeight1( const HELIODETIC_POINT& point, double h0, double r0, double h1, double r1, double* sigma0, double* sigma1 ) const = 0;
		virtual bool								GetEffectiveExtinctionPerCMWithHeight1( const HELIODETIC_POINT& startpoint, const HELIODETIC_POINT& endpoint, double* sigma0, double* sigma1 ) const { nxLog::Record(NXLOG_ERROR,"SKTRAN_TableOpticalProperties_V21::GetEffectiveExtinctionPerCMWithHeight This version of the function was not implemented for this version of code"); return false;}

	//public:
	//	virtual double								TotalExtinctionPerCM(double wavelength, const HELIODETIC_POINT& point) const override { return false; }
	//	virtual bool								GetBRDF(double wavelength, const HELIODETIC_POINT& point, double mu_in, double mu_out, double cosdphi, double* brdf) const override { return false; }
	//	virtual bool								GetBRDFGeodetic(double wavelength, const GEODETIC_INSTANT& point, double mu_in, double mu_out, double cosdphi, double* brdf) const override { return false; }
	//	virtual bool								GetEffectiveExtinctionPerCMWithHeight1(double wavelength, const HELIODETIC_POINT& startpoint, const HELIODETIC_POINT& endpoint, double* sigma0, double* sigma1) const override { return false; }
	//	virtual bool								GetScatteringCoefficientCM2(double wavelength, const HELIODETIC_POINT& point, double cosangle, SKTRAN_PhaseMatrixScalar* scatcoeff) const override { return false; }
	//	virtual double								ScatteringExtinctionPerCM(double wavelength, const HELIODETIC_POINT& point) const override { return false; }
	//	virtual bool								CreateInterpolationForPoint(double wavelength, const HELIODETIC_POINT& point, SKTRAN_TableOpticalProperties_PointCache& interpolator) const override { return false; }

};



/*-----------------------------------------------------------------------------
 *					SKTRANSO_InternalEmissionPropertiesTable		 2015- 3- 3*/
/** **/
/*---------------------------------------------------------------------------*/

class SKTRANSO_InternalEmissionPropertiesTable : public SKTRAN_TableEmission_Base, public SKTRANSO_JindexTableBase
{
	public:
		 											SKTRANSO_InternalEmissionPropertiesTable	() {}
		virtual 								   ~SKTRANSO_InternalEmissionPropertiesTable	() {}
		virtual bool								ConfigureGeometry			( const SKTRAN_SpecsInternal_Base* specs       ) = 0;
		virtual bool								ConfigureOptical			( double wavelen, const SKTRAN_CoordinateTransform_V2* coords, SKTRAN_AtmosphericEmission* emissionstate ) = 0;


};
