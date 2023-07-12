//#include "sktran_common_internals.h"

class SKTRAN_UnitSphere_V2;

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_3D_UnitSphere		2014-2-7*/
/** @ingroup optprop
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_TableOpticalProperties_3D_UnitSphere : public SKTRAN_TableOpticalProperties_Base //public SKTRAN_TableOpticalProperties_V21
{
	protected:
		const SKTRAN_UnitSphere_V2*							m_unitsphere;
		const SKTRAN_GridDefOpticalPropertiesRadii_V21*		m_alts;
		std::vector< std::vector< std::vector<double> > >					m_totalextinction;
		std::vector< std::vector< std::vector<double> > >					m_scatextinction;
		//std::vector<double>                                 m_singlescatt;
		const skBRDF*										m_albedo;
		double												m_wavelen;
		double												m_mjd;
		const SKTRAN_GridDefScatterAngle_V21*				m_scatteranglegrid;
		const SKTRAN_GridDefWavelength*						m_wavelengthgrid;
		mutable std::vector<size_t>							m_speedhelper;						//!< Helps delaunay lookup speed

		size_t												m_wavelindex;
		std::vector<double>									m_wavelengtharray;
		std::vector<double>									m_wavenumber;
		std::vector<double>									m_scatanglevector;
		bool												m_firsttime;
		bool												m_forcecacheupdates;
		bool												m_inonedimmode;

	private:
		virtual double										InterpTable							( const std::vector< std::vector<double> >& table, const HELIODETIC_POINT& loc) const;
		virtual bool										FillTablesAtIndex					( size_t altidx, size_t locidx, SKTRAN_AtmosphericOpticalState_V21& opticalstate );
		virtual bool										FillTablesAtIndexMultiWavel			( size_t altidx, size_t locidx, SKTRAN_AtmosphericOpticalState_V21& opticalstate );
		virtual void										ReleaseResources					();
		bool												DumpTotalExtinction					() const;
		bool                                                GetUniquePointWeights               ( const HELIODETIC_POINT& point, double cosAngle, SKTRAN_GridIndex gridindices[12], double gridweights[12], size_t& numNonZero ) const;
		SKTRAN_GridIndex                                    TableSubToInd                       ( size_t wavidx, SKTRAN_GridIndex locidx, SKTRAN_GridIndex altidx, SKTRAN_GridIndex angidx ) const;

		inline bool                                         GetScatteringMatrixCM2_FromWeights      ( const SKTRAN_GridIndex gridindex[], const double gridweight[], size_t numels, SKTRAN_ScatMat_MIMSNC& pmatrix ) const;
		inline bool                                         GetScatteringCoefficientCM2_FromWeights ( const SKTRAN_GridIndex gridindex[], const double gridweight[], size_t numels, SKTRAN_PhaseMatrixScalar& scatcoeff ) const;

	protected:
		virtual bool										CalcSphereIndices					( const HELIODETIC_POINT& loc, double* weights, size_t* indices, size_t& numindex ) const;
		virtual bool										CalcAltIndices						( const HELIODETIC_POINT& loc, double* weights, size_t* indices, size_t& numindex ) const;
		virtual bool										CalcScatterIndices					( double cosangle, double* weights, SKTRAN_GridIndex* indices, size_t& numindex ) const;
		virtual bool										CalcWavelengthIndices				( double wavelength, double* weights, size_t* indices, size_t& numindex ) const;

	public:
															SKTRAN_TableOpticalProperties_3D_UnitSphere();
		virtual											   ~SKTRAN_TableOpticalProperties_3D_UnitSphere();

		virtual bool										IsOptionTrue						( SKTRAN_TableOpticalProperties_Base::OPTIONSENUM options) const;
		virtual bool										ConfigureGeometry					( const SKTRAN_Specifications_Base* specs       );

		virtual bool										SetAltitudes						( const SKTRAN_GridDefOpticalPropertiesRadii_V21& alts );
		virtual bool										SetUnitSphere						( const SKTRAN_UnitSphere_V2& unitsphere );
		virtual bool										SetScatterGrid						( const SKTRAN_GridDefScatterAngle_V21& scatgrid );
		virtual bool										SetWavelengthGrid					( const SKTRAN_GridDefWavelength& wavgrid );

		virtual bool										ConfigureOptical					( double wavelen, SKTRAN_AtmosphericOpticalState_V21& opticalstate );
//		virtual bool										ConfigureOpticalMT					( double wavelen, SKTRAN_AtmosphericOpticalState_V21& opticalstate );

		// These are easy to implement, just haven't gotten around to it yet
		virtual double										TotalExtinctionPerCM				( const SKTRAN_RayStorage_Base*  r, size_t index ) const {return false;}
		virtual bool										GetAlbedo							( const SKTRAN_RayStorage_Base*  r, size_t index, double* albedo ) const {return false;}
		virtual bool										GetEffectiveExtinctionPerCMWithHeight1( const SKTRAN_RayStorage_Base* r, size_t startPtIndex, double* sigma0, double* sigma1 ) const;
		virtual bool										GetScatteringCoefficientCM2			( const SKTRAN_RayStorage_Base*  r, size_t index, double cosangle,  SKTRAN_PhaseMatrixScalar* phasematrix ) const {return false;}
		virtual double										ScatteringExtinctionPerCM			( const SKTRAN_RayStorage_Base*  r, size_t index ) const {return -9999;}

		virtual double										TotalExtinctionPerCM				( const HELIODETIC_POINT& point ) const;
		virtual bool										GetBRDF								( const HELIODETIC_POINT& point, double mu_in, double mu_out, double cosdphi, double* brdf ) const override;
		virtual bool										GetBRDFGeodetic						( const GEODETIC_INSTANT& point, double mu_in, double mu_out, double cosdphi, double* brdf ) const override;

		virtual bool										GetEffectiveExtinctionPerCMWithHeight1( const HELIODETIC_POINT& startpoint, const HELIODETIC_POINT& endpoint, double* sigma0, double* sigma1 ) const;
		virtual bool										GetScatteringCoefficientCM2			( const HELIODETIC_POINT& point, double cosangle,  SKTRAN_PhaseMatrixScalar* scatcoeff ) const;
		virtual double										ScatteringExtinctionPerCM			( const HELIODETIC_POINT& point ) const;

		bool												SetWavelengths						(std::vector<double>& wavelengths);

//		virtual bool										GetScatteringCoefficientCM2			( const SKTRANSO_JIndex* boundingpoints, size_t numpoints, double cosangle, SKTRAN_PhaseMatrixScalar* coefficient ) const;
//		virtual bool										GetBoundingSpatialPoints			( const HELIODETIC_POINT& location, SKTRANSO_JIndex* interppoints, size_t* numpoints ) const;
//		virtual bool										GetBoundingScatteringPoints			( double cosscatteringangle, SKTRANSO_JIndex* interppoints, size_t* numpoints ) const;
		void												SetForceCacheUpdates				( bool force ) { m_forcecacheupdates = force; }
		
		virtual bool                                GetScatteringMatrixCM2                      ( const HELIODETIC_POINT& point, double cosAngle, SKTRAN_ScatMat_MIMSNC&  pmatrix   ) const override;
		virtual bool                                GetResultOfUnpolarizedScatterCM2            ( const HELIODETIC_POINT& point, double cosAngle, SKTRAN_Stokes_NC& stokesvec ) const override;

        virtual bool                                CreateInterpolationForPoint                 ( const HELIODETIC_POINT& point, SKTRAN_TableOpticalProperties_PointCache& interpolator ) const override;

	// overloads for interpolating on the wavelength grid without setting an internal variable to represent wavelength (thread safe)
	private:
		virtual double										InterpTable(const std::vector< std::vector< std::vector<double> > >& table, double wavelength, const HELIODETIC_POINT& loc) const;
		bool                                                GetUniquePointWeights(double wavelength, const HELIODETIC_POINT& point, double cosAngle, SKTRAN_GridIndex gridindices[24], double gridweights[24], size_t& numNonZero) const;

	public:
		virtual double										TotalExtinctionPerCM(double wavelength, const HELIODETIC_POINT& point) const;
		virtual bool										GetBRDF(double wavelength, const HELIODETIC_POINT& point, double mu_in, double mu_out, double cosdphi, double* brdf) const;
		virtual bool										GetBRDFGeodetic(double wavelength, const GEODETIC_INSTANT& point, double mu_in, double mu_out, double cosdphi, double* brdf) const;
		virtual bool										GetEffectiveExtinctionPerCMWithHeight1(double wavelength, const HELIODETIC_POINT& startpoint, const HELIODETIC_POINT& endpoint, double* sigma0, double* sigma1) const;
		virtual bool										GetEffectiveExtinctionPerCMWithHeight1(double wavelength, const SKTRAN_RayStorage_Base* r, size_t startPtIndex, double* sigma0, double* sigma1) const;
		virtual bool										GetScatteringCoefficientCM2(double wavelength, const HELIODETIC_POINT& point, double cosangle, SKTRAN_PhaseMatrixScalar* scatcoeff) const;
		virtual double										ScatteringExtinctionPerCM(double wavelength, const HELIODETIC_POINT& point) const;
		virtual bool										GetScatteringMatrixCM2(double wavelength, const HELIODETIC_POINT& point, double cosAngle, SKTRAN_ScatMat_MIMSNC&  pmatrix) const;
		virtual bool										GetResultOfUnpolarizedScatterCM2(double wavelength, const HELIODETIC_POINT& point, double cosAngle, SKTRAN_Stokes_NC& stokesvec) const;
		virtual bool										CreateInterpolationForPoint(double wavelength, const HELIODETIC_POINT& point, SKTRAN_TableOpticalProperties_PointCache& interpolator) const;

};


class SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant : public SKTRAN_TableOpticalProperties_3D_UnitSphere
{
	protected:
		virtual bool										CalcAltIndices(const HELIODETIC_POINT& loc, double* weights, size_t* indices, size_t& numindex) const;

	public:
		virtual bool										ConfigureGeometry(const SKTRAN_Specifications_Base* specs);

};