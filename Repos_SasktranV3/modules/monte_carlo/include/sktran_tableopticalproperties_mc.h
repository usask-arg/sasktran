#pragma once

#include <omp.h>
#include "sktran_montecarlo_internals.h"

class SKTRAN_TableOpticalProperties_MCBase
{
	protected:
		std::vector<SKTRAN_PhaseMatrixScalar>					m_scatterCdfMatrices;  // Contains the normalized CDFs for each altitude
		std::vector<SKTRAN_PhaseMatrixScalar>                   m_cdfLookupSpace;      // Memory space used during cdf lookup -- needed to make class threadsafe without tons of memory allocation/deallocation
		std::vector< std::vector< SKTRAN_PhaseMatrixScalar >::iterator >                  m_lookupSpaceIndices;  // Used to access the lookup space. These should really be in some internal class so I don't screw up the access. 
		size_t                                      m_numAngles;
		size_t										m_numalt;
		size_t										m_numloc;
		std::shared_ptr<SKTRAN_TableOpticalProperties_Inelastic_Base> m_inelasticOptProp;

	private:
		inline bool     makeScatterCdfs_trapz                   ( std::unique_ptr< SKTRAN_PolarizationProperties_Base >& scatProps, size_t numAngles );

	protected:
		virtual void	releaseTable							();
		virtual bool	allocateCdf								( const SKTRAN_Specifications_Base *specs );
		virtual bool	allocateCdf								( const SKTRAN_GridDefScatterAngle_V21& scatteranglegrid, const SKTRAN_GridDefOpticalPropertiesRadii_V21& altitudegrid );
		virtual bool	allocateCdf								( const SKTRAN_GridDefScatterAngle_V21& scatteranglegrid, const SKTRAN_GridDefOpticalPropertiesRadii_V21& altitudegrid, const SKTRAN_UnitSphere_V2& unitsphere );
		virtual bool	allocateCdf								( const SKTRAN_GridDefScatterAngle_V21& scatteranglegrid, const SKTRAN_GridDefOpticalPropertiesRadii_V21& altitudegrid, const SKTRAN_UnitSphere_V2& unitsphere, const SKTRAN_GridDefWavelength& wavelengthgrid );
		//virtual bool	makeScatterCdfs							( nx2dArray<SKTRAN_PhaseMatrixScalar>** const singleScatt, const size_t& numAlts, const size_t& numAngles);
		//virtual bool	makeScatterCdfs							( const nx3dArray<SKTRAN_PhaseMatrixScalar>& singlescatt, size_t numalt, size_t numangle, size_t numloc );
		virtual bool    makeScatterCdfs                         ( std::unique_ptr< SKTRAN_PolarizationProperties_Base >& scatProps, size_t numAngles );
		virtual bool	GetCosScatteringAngleWeights			( const double& randNum, SKTRAN_GridIndex& loindex, double& loweight, SKTRAN_GridIndex& hiindex, double& hiweight ) const;
		virtual bool	GetCosScatteringAngleWeights			( double randNum, size_t* altindex, double* altweight, size_t numalt, size_t* locindex, double* locweight, size_t numloc, size_t* scatindex, double* scatweight ) const;
		virtual bool	GetCosScatteringAngleWeights			( double randNum, size_t* wavindex, double*wavweight, size_t numwav, size_t* altindex, double* altweight, size_t numalt, size_t* locindex, double* locweight, size_t numloc, size_t* scatindex, double* scatweight) const;
		virtual bool	GetCosScattAngleCdf						( const SKTRAN_GridIndex &lowercell, const double &lowerw, const SKTRAN_GridIndex &uppercell, const double &upperw, std::vector<SKTRAN_PhaseMatrixScalar>::iterator cdf ) const; // Get the scattering cdf. Normally this is normalized, but may not be in odd conditions (e.g. if tracing rays above where optical properties are defined, cdf=0)
		virtual bool	GetCosScattAngleCdf						( size_t* altindex, double* altweight, size_t numalt, size_t* locindex, double* locweight, size_t numloc, std::vector<SKTRAN_PhaseMatrixScalar>::iterator cdf ) const; // Get the scattering cdf. Normally this is normalized, but may not be in odd conditions (e.g. if tracing rays above where optical properties are defined, cdf=0)
		virtual bool	GetCosScattAngleCdf						( size_t* wavindex, double* wavweight, size_t numwav, size_t* altindex, double* altweight, size_t numalt, size_t* locindex, double* locweight, size_t numloc, std::vector<SKTRAN_PhaseMatrixScalar>::iterator cdf) const; 

		virtual bool	FindScatterAngleBoundingIndices			( std::vector< SKTRAN_PhaseMatrixScalar >::const_iterator cdf, double target, size_t& lowIndex, size_t& uppIndex, double& lowerw, double& upperw) const ;
	
	public:
		virtual bool												SetInelasticProperties( std::shared_ptr<SKTRAN_TableOpticalProperties_Inelastic_Base>& inelasticOptProp ) { m_inelasticOptProp = inelasticOptProp; return m_inelasticOptProp != nullptr; }
		virtual SKTRAN_TableOpticalProperties_Inelastic_Base const* InelasticProperties( ) const { return m_inelasticOptProp.get(); }
	public:
						SKTRAN_TableOpticalProperties_MCBase	();
		virtual			~SKTRAN_TableOpticalProperties_MCBase	();
		virtual bool	GetCosScatteringAngle					( const HELIODETIC_POINT& point, const double& randNum, double &cosScatAngle, SKTRAN_ScatMat_MIMSNC* pmatrix = nullptr ) const = 0; // pmatrix is the phase matrix for scatter through -cosScatAngle
		virtual bool	GetCosScatteringAngle					( const double& wavelength, const HELIODETIC_POINT& point, const double& randNum, double &cosScatAngle, SKTRAN_ScatMat_MIMSNC* pmatrix = nullptr ) const = 0; // pmatrix is the phase matrix for scatter through -cosScatAngle
		virtual bool	MakeThreadsafeFor                       ( size_t numThreads ) = 0;

};

class SKTRAN_TableOpticalProperties_1D_Height_MC : public SKTRAN_TableOpticalProperties_1D_Height_V3, public SKTRAN_TableOpticalProperties_MCBase
{

	protected:
		virtual bool    AllocateCdfLookupSpace                      ( size_t numThreads );
		void            ReleaseResources                            ( );
	//	void            RemoveScattExtinctionFromPhaseInfo          ( );

	public:
						SKTRAN_TableOpticalProperties_1D_Height_MC	( );
		virtual		   ~SKTRAN_TableOpticalProperties_1D_Height_MC	( );
		//virtual bool	Allocate( size_t numcells, size_t numangles );
		virtual bool	ConfigureGeometry							( const SKTRAN_Specifications_Base* specs) override;
		virtual bool	ConfigureOptical							( double wavelen, SKTRAN_AtmosphericOpticalState_V21& opticalstate ) override;
		virtual bool	GetCosScatteringAngle						( const HELIODETIC_POINT& point, const double& randNum, double &cosScatAngle, SKTRAN_ScatMat_MIMSNC* pmatrix = nullptr ) const override;
		virtual bool	GetCosScatteringAngle						( const double& wavelength, const HELIODETIC_POINT& point, const double& randNum, double &cosScatAngle, SKTRAN_ScatMat_MIMSNC* pmatrix = nullptr ) const override { return false; }
		virtual bool    MakeThreadsafeFor                           ( size_t numThreads ) override;

};

class SKTRAN_TableOpticalProperties_3D_UnitSphere_MC : public SKTRAN_TableOpticalProperties_3D_UnitSphere, public SKTRAN_TableOpticalProperties_MCBase
{
	friend class SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere;

	protected:
		virtual bool    AllocateCdfLookupSpace                      ( size_t numThreads );
		void            ReleaseResources                            ( );

	public:
						SKTRAN_TableOpticalProperties_3D_UnitSphere_MC	( );
		virtual		   ~SKTRAN_TableOpticalProperties_3D_UnitSphere_MC	( );
		virtual bool	ConfigureGeometry								( const SKTRAN_Specifications_Base* specs) override;
		virtual bool	ConfigureOptical								( double wavelen, SKTRAN_AtmosphericOpticalState_V21& opticalstate ) override;
		virtual bool	GetCosScatteringAngle							( const HELIODETIC_POINT& point, const double& randNum, double &cosScatAngle, SKTRAN_ScatMat_MIMSNC* pmatrix = nullptr ) const override;
		virtual bool	GetCosScatteringAngle							( const double& wavelength, const HELIODETIC_POINT& point, const double& randNum, double &cosScatAngle, SKTRAN_ScatMat_MIMSNC* pmatrix = nullptr ) const override;
		virtual bool    MakeThreadsafeFor								( size_t numThreads ) override;		
		

};

class SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant_MC : public SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant, public SKTRAN_TableOpticalProperties_MCBase
{
	// TODO find a better way to do this (a lot of copy/pasted code)
	friend class SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere;

	protected:
		virtual bool    AllocateCdfLookupSpace                      ( size_t numThreads );
		void            ReleaseResources                            ( );

	public:
						SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant_MC	( );
		virtual		   ~SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant_MC	( );
		virtual bool	ConfigureGeometry								( const SKTRAN_Specifications_Base* specs) override;
		virtual bool	ConfigureOptical								( double wavelen, SKTRAN_AtmosphericOpticalState_V21& opticalstate ) override;
		virtual bool	GetCosScatteringAngle							( const HELIODETIC_POINT& point, const double& randNum, double &cosScatAngle, SKTRAN_ScatMat_MIMSNC* pmatrix = nullptr ) const override;
		virtual bool	GetCosScatteringAngle							( const double& wavelength, const HELIODETIC_POINT& point, const double& randNum, double &cosScatAngle, SKTRAN_ScatMat_MIMSNC* pmatrix = nullptr ) const override;
		virtual bool    MakeThreadsafeFor								( size_t numThreads ) override;		
		
};

//
//class SKTRAN_TableOpticalProperties_1D_Height_Polarized_MC : public SKTRAN_TableOpticalProperties_1D_Height_Polarized_V3, public SKTRAN_TableOpticalProperties_MCBase
//{
//
//	protected:
//		virtual bool    AllocateCdfLookupSpace                      ( size_t numThreads );
//		void            ReleaseResources                            ( );
//		//void            RemoveScattExtinctionFromPhaseInfo          ( );
//
//	public:
//						SKTRAN_TableOpticalProperties_1D_Height_Polarized_MC	( );
//		virtual		   ~SKTRAN_TableOpticalProperties_1D_Height_Polarized_MC	( );
//		//virtual bool	Allocate( size_t numcells, size_t numangles );
//		virtual bool	ConfigureGeometry							( const SKTRAN_Specifications_Base* specs);
//		virtual bool	ConfigureOptical							( double wavelen, SKTRAN_AtmosphericOpticalState_V21& opticalstate );
//		virtual bool	GetCosScatteringAngle						( const HELIODETIC_POINT& point, const double& randNum, double &cosScatAngle, skRTPhaseMatrix* pmatrix = NULL) const;
//		virtual bool    MakeThreadsafeFor                           ( size_t numThreads );
//
//				
//};
//
//