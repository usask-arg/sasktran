#pragma once

class SKTRAN_MCPhoton_Base;

class SKTRAN_EmissionTable_Base : public SKTRAN_Source_Term 
{
    public:

    virtual      ~SKTRAN_EmissionTable_Base					( );
    virtual bool  ConfigureOptical							( double wavelen, SKTRAN_AtmosphericEmission* emissionObject, const SKTRAN_TableOpticalProperties_Base * opticaltable  ) = 0;

    virtual bool  SourceTermAtPoint							( const SKTRAN_SourceTermQueryObject_Base& qobj, double*           source   ) const override = 0;
    virtual bool  SourceTermAtPoint							( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC* source   ) const override;

    virtual bool  GroundSourceAtPoint						( const SKTRAN_SourceTermQueryObject_Base& qobj, double*           source   ) const override = 0;
    virtual bool  GroundSourceAtPoint						( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC* source   ) const override;

    virtual bool  MonteCarlo_SingleScatteredRadianceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, double&           radiance ) const override = 0;
    virtual bool  MonteCarlo_SingleScatteredRadianceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC& radiance ) const override;

    virtual bool  MonteCarlo_GroundScatteredRadianceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, double&           radiance ) const override = 0;
    virtual bool  MonteCarlo_GroundScatteredRadianceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC& radiance ) const override; 

    virtual bool  InitializeGeometry                        ( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords,
                                                              const SKTRAN_GridDefRayTracingShells_V21&             altgrid,
                                                              const SKTRAN_GridDefCosSZA_V21&                       cosszagrid,
                                                              const SKTRAN_GridDefSLON_V21&                         slongrid    ) = 0;
};


class SKTRAN_EmissionTable_DoNothing : public SKTRAN_EmissionTable_Base
{
    protected:
        virtual bool  FillTable_ClassSpecific					( ) { return true;}

    public:
                      SKTRAN_EmissionTable_DoNothing			 ( ) { ;}
                     ~SKTRAN_EmissionTable_DoNothing			 ( ) { ;}
        virtual bool  MakeThreadSafeFor							 ( size_t numThreads ){ return true;}
        virtual bool  ConfigureOptical							 ( double wavelen, SKTRAN_AtmosphericEmission* emissionObject, const SKTRAN_TableOpticalProperties_Base * opticaltable  ) override { return true;}
        virtual bool  SourceTermAtPoint							 ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_StokesScalar* source   ) const override { *source       = 0.0;  return true;}
		virtual bool  GroundSourceAtPoint						 ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_StokesScalar* source   ) const override { *source       = 0.0;  return true;}
		virtual bool  MonteCarlo_SingleScatteredRadianceAtPoint  ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_StokesScalar& radiance ) const override { radiance      = 0.0;  return true;}
		virtual bool  MonteCarlo_GroundScatteredRadianceAtPoint  ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_StokesScalar& radiance ) const override { radiance      = 0.0;  return true;}
        virtual bool  SourceTermAtPoint							 ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC*    source   ) const override { source  ->SetTo(0.0); return true;} 
        virtual bool  GroundSourceAtPoint						 ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC*    source   ) const override { source  ->SetTo(0.0); return true;} 
        virtual bool  MonteCarlo_SingleScatteredRadianceAtPoint  ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC&    radiance ) const override { radiance .SetTo(0.0); return true;}
        virtual bool  MonteCarlo_GroundScatteredRadianceAtPoint  ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC&    radiance ) const override { radiance .SetTo(0.0); return true;}

		virtual bool  SourceTermAtPoint							 ( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_StokesScalar* source   ) const override { *source       = 0.0;  return true;}
		virtual bool  GroundSourceAtPoint						 ( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_StokesScalar* source   ) const override { *source       = 0.0;  return true;}
		virtual bool  MonteCarlo_SingleScatteredRadianceAtPoint  ( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_StokesScalar& radiance ) const override { radiance      = 0.0;  return true;}
		virtual bool  MonteCarlo_GroundScatteredRadianceAtPoint  ( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_StokesScalar& radiance ) const override { radiance      = 0.0;  return true;}
        virtual bool  SourceTermAtPoint							 ( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC*    source   ) const override { source  ->SetTo(0.0); return true;}
        virtual bool  GroundSourceAtPoint						 ( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC*    source   ) const override { source  ->SetTo(0.0); return true;}
        virtual bool  MonteCarlo_SingleScatteredRadianceAtPoint  ( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC&    radiance ) const override { radiance .SetTo(0.0); return true;}
        virtual bool  MonteCarlo_GroundScatteredRadianceAtPoint  ( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC&    radiance ) const override { radiance .SetTo(0.0); return true;}

		virtual bool  MonteCarlo_SingleScatteredRadianceAtPoint  ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_MCPhoton_Base* photon) const override { return true; }
		virtual bool  MonteCarlo_GroundScatteredRadianceAtPoint  ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_MCPhoton_Base* photon) const override { return true; }

        virtual bool  InitializeGeometry                         ( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords,
                                                                   const SKTRAN_GridDefRayTracingShells_V21&             altgrid,
                                                                   const SKTRAN_GridDefCosSZA_V21&                       cosszagrid,
                                                                   const SKTRAN_GridDefSLON_V21&                         slongrid    ) override { return true;}
};


class SKTRAN_EmissionTable_NoTable : public SKTRAN_EmissionTable_Base 
{
    private:
        mutable SKTRAN_AtmosphericEmission      m_table;
        //skClimatology_MSIS90*                   m_neutralatmosphere;
        skSolarSpectrum_SAO2010*                m_solarSpectrum;
        skEmission_Tabulated_HeightWavelength*  m_emission;
        const SKTRAN_CoordinateTransform_V2*    m_coords;
        double                                  m_mjd;
        double                                  m_wavenumber;

    std::string m_filename;

    public:
                  SKTRAN_EmissionTable_NoTable ( );
    virtual      ~SKTRAN_EmissionTable_NoTable ( );
                  
    virtual bool  ConfigureOptical							( double wavelen, SKTRAN_AtmosphericEmission* emissionObject, const SKTRAN_TableOpticalProperties_Base * opticaltable  ) override;
	virtual bool  SourceTermAtPoint							( const SKTRAN_SourceTermQueryObject_Base& qobj, double* source   ) const override;
    virtual bool  GroundSourceAtPoint						( const SKTRAN_SourceTermQueryObject_Base& qobj, double* source   ) const override;
    virtual bool  MonteCarlo_SingleScatteredRadianceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, double& radiance ) const override;
    virtual bool  MonteCarlo_GroundScatteredRadianceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, double& radiance ) const override;

    void SetTableFile       ( std::string& filename );
    bool Initialize         ( );
    void InitializeGeometry ( const SKTRAN_CoordinateTransform_V2* coords );

};


class SKTRAN_EmissionTable_1D : public SKTRAN_EmissionTable_Base 
{
    private:
        std::vector<double>                 m_emission;
        std::vector<double>                 m_emissionSample;
        SKTRAN_GridDefRayTracingShells_V21  m_altgrid;
        double                              m_groundemission;
        GEODETIC_INSTANT                    m_point;


    public:
                  SKTRAN_EmissionTable_1D        ( );
    virtual      ~SKTRAN_EmissionTable_1D        ( );

    virtual bool  ConfigureOptical							( double wavelen, SKTRAN_AtmosphericEmission* emissionObject, const SKTRAN_TableOpticalProperties_Base * opticaltable ) override;

    virtual bool  SourceTermAtPoint							( const SKTRAN_SourceTermQueryObject_Base& qobj, double* source   ) const override;
    virtual bool  GroundSourceAtPoint						( const SKTRAN_SourceTermQueryObject_Base& qobj, double* source   ) const override;
    virtual bool  MonteCarlo_SingleScatteredRadianceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, double& radiance ) const override;
    virtual bool  MonteCarlo_GroundScatteredRadianceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, double& radiance ) const override;

    virtual bool  InitializeGeometry                         ( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords,
                                                               const SKTRAN_GridDefRayTracingShells_V21&             altgrid,
                                                               const SKTRAN_GridDefCosSZA_V21&                       cosszagrid,
                                                               const SKTRAN_GridDefSLON_V21&                         slongrid   ) override;
};


