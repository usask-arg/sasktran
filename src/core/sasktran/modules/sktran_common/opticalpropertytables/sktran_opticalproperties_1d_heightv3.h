//#include "sktran_common_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height_V3		2007-11-20*/
/** @ingroup optprop
 *	A class used to cache the optical properties of the atmosphere. This class
 *	implements optical properties which are only a function of altitude.
 *	NOte the classes SKTRAN_QuadratureScatteringMatrixxxxxx are also affected by
 *	by this.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_TableOpticalProperties_1D_Height_V3 : public SKTRAN_TableOpticalProperties_Base //public SKTRAN_TableOpticalProperties_V21
{

	protected:
		size_t											m_numheighttoindex;			//!< Number of elements in the height to index table (typically 10,000)
		size_t*											m_heighttoindextable;		//!< A table [m_numheighttoindex] that quickly converts a radius to an index with specified resolution
		double											m_heightindexresolution;	//!< The resolution o fthe radius to index rsolution, typically 10 meters
		double											m_minheight;
//		nx1dArray<double>								m_effectiveextinction;		//!< Array [numshells] of effectiveextinction .
		double											m_wavelen;
		double											m_mjd;
		const skBRDF*									m_albedo;
		const SKTRAN_GridDefOpticalPropertiesRadii_V21*	m_altitudegrid;				//!< The grid that defines the spatial grid for specifying optical properties
		const SKTRAN_GridDefScatterAngle_V21*			m_scatteranglegrid;			//!< The grid that defines the angular grid for specifying scattered rays.
		std::vector<double>**							m_scatextinction;			//!< JTW: Array [numshells] of scattering extinction .
		std::vector<double>**							m_extinction;				//!< Array [numshells] of extinction .
		//nx2dArray<SKTRAN_PhaseMatrixScalar>**			m_singleScatt;				//!< Array (numshells, numscatterangles) Indexed by shell and angle.
		//std::vector<double>								m_singleScatt;				//!< Array (numshells, numscatterangles) Indexed by shell and angle.
		//const SKTRAN_CoordinateTransform_V2*			m_coordinates;				//!< 
		double											M_GROUNDTOLERANCE;			//!< Treat points above this negative value as though they were at ground (i.e., forgive rounding error)

	private:
		bool										ConfigureAltitudeToIndexTable		();
		double                                      GetExtinctionPerCM					( double r, const std::vector<double>&	extinction  ) const;
		bool                                        GetUniquePointWeights               ( const HELIODETIC_POINT& point, double cosAngle, SKTRAN_GridIndex gridindices[4], double gridweights[4], size_t& numNonZero ) const;

		// Interface for rays using minimum container storage
		/* DEPRECATED FUNCTION */ double			TotalExtinctionPerCM				( const SKTRAN_RayStorage_Base* r, size_t index ) const;
		/* DEPRECATED FUNCTION */ bool				GetAlbedo							( const SKTRAN_RayStorage_Base* r, size_t index, double* albedo ) const;
		/* DEPRECATED FUNCTION */ bool				GetEffectiveExtinctionPerCMWithHeight1( const SKTRAN_RayStorage_Base* r, size_t startPtIndex, double* sigma0, double* sigma1 ) const;
		/* DEPRECATED FUNCTION */ bool				GetScatteringCoefficientCM2			( const SKTRAN_RayStorage_Base* r, size_t index, double cosangle,  SKTRAN_PhaseMatrixScalar* phasematrix ) const;
		/* DEPRECATED FUNCTION */ double			ScatteringExtinctionPerCM			( const SKTRAN_RayStorage_Base* r, size_t index ) const ;

	protected:
		virtual bool								Allocate							( size_t numcells, size_t numangles );
		        void								ReleaseResources					();
		virtual void								ReleaseObjects						();
		//double										GetExtinctionPerCM				( const HELIODETIC_POINT& location, const nx1dArray<double>& extinction  ) const;
		virtual bool								IndexOfPointBelowOrEqual			( double h0, size_t* lowindex ) const;
		bool										IndexOfPointEqualOrAbove			( double h0, size_t* hihindex ) const;
		virtual void								nxdebugFunction();
		inline SKTRAN_GridIndex                     TableSubToInd                       ( SKTRAN_GridIndex altidx, SKTRAN_GridIndex angidx ) const;

	private:

	public:
													SKTRAN_TableOpticalProperties_1D_Height_V3		();
		virtual 								   ~SKTRAN_TableOpticalProperties_1D_Height_V3		();
		virtual bool								ConfigureGeometry								( const SKTRAN_GridDefScatterAngle_V21& scatteranglegrid, const SKTRAN_GridDefOpticalPropertiesRadii_V21& altitudegrid);
		virtual bool								ConfigureGeometry								( const SKTRAN_Specifications_Base* specs ) ;
		virtual bool								GetEffectiveExtinctionPerCMWithHeightShell1		( const HELIODETIC_POINT& point, double h0, double r0, double h1, double r1, double* sigma0, double* sigma1 ) const;
		

		// Interface for general point-by-point access
		virtual bool								ConfigureOptical								( double wavelen, SKTRAN_AtmosphericOpticalState_V21& opticalstate ) override;
		virtual bool								IsOptionTrue									( SKTRAN_TableOpticalProperties_Base::OPTIONSENUM options) const override;
		virtual double								TotalExtinctionPerCM							( const HELIODETIC_POINT& point  ) const override;				//, double* shell_sigma_, double* shell_r, bool getshellbelow_radius ) const;
		virtual bool								GetBRDF											( const HELIODETIC_POINT& point, double mu_in, double mu_out, double cosdphi, double* brdf ) const override;
		virtual bool								GetBRDFGeodetic									( const GEODETIC_INSTANT& point, double mu_in, double mu_out, double cosdphi, double* brdf ) const override;

		virtual bool								GetEffectiveExtinctionPerCMWithHeight1			( const HELIODETIC_POINT& startpoint, const HELIODETIC_POINT& endpoint, double* sigma0, double* sigma1 ) const override;
		virtual bool								GetScatteringCoefficientCM2						( const HELIODETIC_POINT& point, double cosangle,  SKTRAN_PhaseMatrixScalar* phasematrix ) const override;
		virtual double								ScatteringExtinctionPerCM						( const HELIODETIC_POINT& point  ) const override;

		// Currently only implemented for point-by-point access -- would be very useful for minimum container interface
		virtual bool								GetLinearExtinctionPerCMVector				( const std::vector< HELIODETIC_POINT>& quadpoints, std::vector<double>& sigmak, std::vector<double>& sigmaf, size_t numpoints ) const ;


//		virtual size_t								GetNumPerturbations							() const { return 0;}
		size_t										NumShells									() const { return (m_altitudegrid     == NULL)? 0 : m_altitudegrid->NumAltitudes();}
		size_t										NumScatterAngles							() const { return (m_scatteranglegrid == NULL)? 0 : m_scatteranglegrid->NumAngles();}
		
		virtual bool                                GetScatteringMatrixCM2						( const HELIODETIC_POINT& point, double cosAngle, SKTRAN_ScatMat_MIMSNC&  pmatrix   ) const override;
		virtual bool                                GetResultOfUnpolarizedScatterCM2            ( const HELIODETIC_POINT& point, double cosAngle, SKTRAN_Stokes_NC& stokesvec ) const override;

        virtual bool                                CreateInterpolationForPoint                 ( const HELIODETIC_POINT& point, SKTRAN_TableOpticalProperties_PointCache& interpolator ) const override;

	public:
		virtual double								TotalExtinctionPerCM(double wavelength, const HELIODETIC_POINT& point) const override { nxLog::Record(NXLOG_ERROR, "SKTRAN_TableOpticalProperties_1D_Height_V3::TotalExtinctionPerCM, the wavelength-dependent function has not been implemented."); return std::numeric_limits<double>::infinity(); }
		virtual bool								GetBRDF(double wavelength, const HELIODETIC_POINT& point, double mu_in, double mu_out, double cosdphi, double* brdf) const override { return false; }
		virtual bool								GetBRDFGeodetic(double wavelength, const GEODETIC_INSTANT& point, double mu_in, double mu_out, double cosdphi, double* brdf) const override { return false; }
		virtual bool								GetEffectiveExtinctionPerCMWithHeight1(double wavelength, const HELIODETIC_POINT& startpoint, const HELIODETIC_POINT& endpoint, double* sigma0, double* sigma1) const override { return false; }
		virtual bool								GetEffectiveExtinctionPerCMWithHeight1(double wavelength, const SKTRAN_RayStorage_Base* r, size_t startPtIndex, double* sigma0, double* sigma1) const override { return false; }
		virtual bool								GetScatteringCoefficientCM2(double wavelength, const HELIODETIC_POINT& point, double cosangle, SKTRAN_PhaseMatrixScalar* scatcoeff) const override { return false; }
		virtual double								ScatteringExtinctionPerCM(double wavelength, const HELIODETIC_POINT& point) const { nxLog::Record(NXLOG_ERROR, "SKTRAN_TableOpticalProperties_1D_Height_V3::ScatteringExtinctionPerCM, the wavelength-dependent function has not been implemented."); return std::numeric_limits<double>::infinity(); }
		virtual bool								GetScatteringMatrixCM2(double wavelength, const HELIODETIC_POINT& point, double cosAngle, SKTRAN_ScatMat_MIMSNC& pmatrix) const override { return false; }
		virtual bool								GetResultOfUnpolarizedScatterCM2(double wavelength, const HELIODETIC_POINT& point, double cosAngle, SKTRAN_Stokes_NC& stokesvec) const override { return false; }
		virtual bool								CreateInterpolationForPoint(double wavelength, const HELIODETIC_POINT& point, SKTRAN_TableOpticalProperties_PointCache& interpolator) const override { return false; }
};

