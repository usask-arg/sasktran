#include "sktran_montecarlo_internals.h"

class SKTRAN_GridDefAirMassFactorShells : public SKTRAN_GridDefRayTracingShells_V21
{
private:
	bool						m_extendedToGround;
	bool						m_extendedToTOA;

public:
	SKTRAN_GridDefAirMassFactorShells();
	~SKTRAN_GridDefAirMassFactorShells() {}

public:
	bool						ConfigureHeights(const double* shellAlts, size_t numshells);
	bool						ConfigureHeights(const std::vector<double>& shellAlts);
	bool						ConfigureHeights(const double* shellAlts, size_t numshells, double surfaceHeight, double toaHeight);
	bool						ConfigureHeights(const std::vector<double>& shellAlts, double surfaceHeight, double toaHeight);

public:
	bool						ExtendedToGround() const { return m_extendedToGround; }
	bool						ExtendedToTOA() const { return m_extendedToTOA; }
};


class SKTRAN_MCAirMassFactorCalculator_Base
{
	protected:
		std::shared_ptr<const SKTRAN_CoordinateTransform_V2> m_coords;
		std::shared_ptr<SKTRAN_GridDefAirMassFactorShells> m_shellGrid_amf;
		SKTRAN_TableOpticalProperties_Base* m_amfopticalpropertiestable;
		SKTRAN_OpticalPropertiesIntegrator_Base* m_amfopticalpropsintegrator;

		std::shared_ptr<SKTRAN_RayFactory_Base>		m_amfRayFactory_los;
		std::shared_ptr<SKTRAN_RayFactory_Base>		m_amfRayFactory_secondary;
		std::shared_ptr<SKTRAN_RayFactory_Base>		m_amfRayFactory_solar;

		bool m_loscurved;
		bool m_solarcurved;
		bool m_secondarycurved;

		std::unique_ptr<const SKTRAN_RayOptical_Base>			m_losray;
		std::vector<std::unique_ptr<SKTRAN_RayOptical_Base>>	m_solarray;
		std::vector<std::unique_ptr<SKTRAN_RayOptical_Base>>	m_secondaryray;

		CLIMATOLOGY_HANDLE m_amfspecieshandle;

	protected:
		const SKTRAN_CoordinateTransform_V2*	Coordinates() const { return m_coords.get(); }
		bool									TraceRay( const HELIODETIC_VECTOR& observer, const HELIODETIC_UNITVECTOR& look, bool curved, bool optical, SKTRAN_RayOptical_Base* ray ) const;
		bool									FindScatterPoint(const SKTRAN_RayOptical_Base* ray, const HELIODETIC_POINT& scatterPoint, size_t& numCells, double& finalCellWeight) const;
		bool									IndexRayCells( const SKTRAN_RayOptical_Base* ray, size_t numcells, std::vector<size_t>& amfCellIndices ) const;
	
	public:
								SKTRAN_MCAirMassFactorCalculator_Base();
							   ~SKTRAN_MCAirMassFactorCalculator_Base();
		bool					ReleaseResources();

		bool					SetCoords						( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords )				{ m_coords = coords; return true; }
		bool					SetSpecies						( CLIMATOLOGY_HANDLE speciesHandle )											{ m_amfspecieshandle = speciesHandle; return true; }
		bool					SetRayTracingShells				( std::shared_ptr<SKTRAN_GridDefAirMassFactorShells> amfshells )				{ m_shellGrid_amf = amfshells; return true; }
		bool					SetOpticalPropertiesTable		( SKTRAN_TableOpticalProperties_Base* opticalpropertiestable )					{ m_amfopticalpropertiestable = opticalpropertiestable; if (opticalpropertiestable != nullptr) opticalpropertiestable->AddRef(); return true; }
		bool					SetOpticalPropertiesIntegrator	( SKTRAN_OpticalPropertiesIntegrator_Base* opticalpropsintegrator )				{ m_amfopticalpropsintegrator = opticalpropsintegrator; if (opticalpropsintegrator != nullptr) opticalpropsintegrator->AddRef(); return true; }
		bool					SetRayFactory_LOS				( std::shared_ptr<SKTRAN_RayFactory_Base>	rayFactory_los, bool curved )		{ m_amfRayFactory_los = rayFactory_los; m_loscurved = curved; return true; }
		bool					SetRayFactory_SOLAR				( std::shared_ptr<SKTRAN_RayFactory_Base>	rayFactory_solar, bool curved )		{ m_amfRayFactory_solar = rayFactory_solar; m_solarcurved = curved; return true; }
		bool					SetRayFactory_SECONDARY			( std::shared_ptr<SKTRAN_RayFactory_Base>	rayFactory_secondary, bool curved )	{ m_amfRayFactory_secondary = rayFactory_secondary; m_secondarycurved = curved; return true; }

		virtual size_t				NumAMFCells() const;
		virtual std::vector<double>	AMFShellHeights() const;
		
		virtual bool			CalculateOpticalPropertiesTable	( double wavelen, SKTRAN_AtmosphericOpticalState_V21* opticalstate, bool userupdateclimatology ) = 0;
		virtual bool			InitializeLogger				( SKTRAN_MCAirMassFactorLogger* logger ) const = 0;
		virtual bool			AllocatePhotons					( std::vector<std::unique_ptr<SKTRAN_MCPhoton_Base>>& mcphotons ) const = 0;
		virtual bool			ClearPhoton						( SKTRAN_MCPhoton_Base* mcphoton ) const = 0;
		virtual bool			AllocateRayOptical				( size_t numthreads ) = 0;
		virtual bool			TraceMotherRay					( SKTRAN_RayOptical_Base const* motherRay ) = 0;

		virtual bool			CalculateSlantContribution		( size_t order, SKTRAN_RayOptical_Base const* scatterRay, const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* mcphoton, size_t threadid ) = 0;


};

class SKTRAN_MCAirMassFactorCalculator_Length : public SKTRAN_MCAirMassFactorCalculator_Base
{
	public:
								SKTRAN_MCAirMassFactorCalculator_Length();
							   ~SKTRAN_MCAirMassFactorCalculator_Length() {}

		virtual bool			CalculateOpticalPropertiesTable ( double wavelen, SKTRAN_AtmosphericOpticalState_V21* opticalstate, bool userupdateclimatology ) override { return true; }
		virtual bool			InitializeLogger				( SKTRAN_MCAirMassFactorLogger* logger ) const override;
		virtual bool			AllocatePhotons					( std::vector<std::unique_ptr<SKTRAN_MCPhoton_Base>>& mcphotons ) const override;
		virtual bool			ClearPhoton						( SKTRAN_MCPhoton_Base* mcphoton ) const override;
		virtual bool			AllocateRayOptical				( size_t numthreads) override;
		virtual bool			TraceMotherRay					( SKTRAN_RayOptical_Base const* motherRay) override;

		virtual bool			CalculateSlantContribution		( size_t order, SKTRAN_RayOptical_Base const* scatterRay, const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* mcphoton, size_t threadid ) override;
};

class SKTRAN_MCAirMassFactorCalculator_OpticalDepth : public SKTRAN_MCAirMassFactorCalculator_Base
{
private:
	bool					PartialOpticalDepth( SKTRAN_RayOptical_Base const* ray, size_t cellindex, double fraction, double& opticaldepth ) const;
	bool					VerticalOpticalDepth( std::vector<double>& vod ) const;

public:
							SKTRAN_MCAirMassFactorCalculator_OpticalDepth();
						   ~SKTRAN_MCAirMassFactorCalculator_OpticalDepth() {}

	virtual bool			CalculateOpticalPropertiesTable ( double wavelen, SKTRAN_AtmosphericOpticalState_V21* opticalstate, bool userupdateclimatology ) override;
	virtual bool			InitializeLogger				( SKTRAN_MCAirMassFactorLogger* logger ) const override;
	virtual bool			AllocatePhotons					( std::vector<std::unique_ptr<SKTRAN_MCPhoton_Base>>& mcphotons ) const override;
	virtual bool			ClearPhoton						( SKTRAN_MCPhoton_Base* mcphoton ) const override;
	virtual bool			AllocateRayOptical				( size_t numthreads ) override;
	virtual bool			TraceMotherRay					( SKTRAN_RayOptical_Base const* motherRay ) override;

	virtual bool			CalculateSlantContribution		( size_t order, SKTRAN_RayOptical_Base const* scatterRay, const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* mcphoton, size_t threadid ) override;
};

class SKTRAN_MCAirMassFactorCalculator_DoNothing : public SKTRAN_MCAirMassFactorCalculator_Base
{
public:
							SKTRAN_MCAirMassFactorCalculator_DoNothing() {}
						   ~SKTRAN_MCAirMassFactorCalculator_DoNothing() {}

	virtual size_t				NumAMFCells() const override { return 0; }
	virtual std::vector<double>	AMFShellHeights() const override { std::vector<double> v; return v; }

	virtual bool			CalculateOpticalPropertiesTable	( double wavelen, SKTRAN_AtmosphericOpticalState_V21* opticalstate, bool userupdateclimatology ) override { return true; }
	virtual bool			InitializeLogger				( SKTRAN_MCAirMassFactorLogger* logger ) const override;
	virtual bool			AllocatePhotons					( std::vector<std::unique_ptr<SKTRAN_MCPhoton_Base>>& mcphotons ) const override { return true; }
	virtual bool			ClearPhoton						( SKTRAN_MCPhoton_Base* mcphoton ) const override { return true; }
	virtual bool			AllocateRayOptical				( size_t numthreads ) override { return true; }
	virtual bool			TraceMotherRay					( SKTRAN_RayOptical_Base const* motherRay ) override { return true; }

	virtual bool			CalculateSlantContribution		( size_t order, SKTRAN_RayOptical_Base const* scatterRay, const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* mcphoton, size_t threadid) override { return true; }

};
