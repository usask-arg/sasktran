


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_3D		2014-2-7*/
/** @ingroup solartransmission
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_SolarTransmission_3D : public SKTRAN_SolarTransmission_2D
{
	private:
		SKTRAN_GridDefSLON_V21							m_slongrid;
		mutable std::vector<SKTRAN_StokesScalar>        m_transmission;
		std::vector<double> 							m_deflectionangle;
		size_t 	m_numcosangles;
		size_t  m_numslon;
		size_t  m_numalts;

		bool m_prefilltable;
		bool m_refractionenabled;
		SKTRAN_GridDefCosSZA_V21						m_prefillgrid;

	private:
		void											ReleaseResources();
		bool											SlonWeights							( double slon, double* slonweights, size_t* slonindex, size_t& numindex ) const;
		bool											FillTableAtIndex					( size_t szaidx, size_t slonidx, size_t altidx ) const;
		inline size_t                                          sub2ind                             ( size_t szaidx, size_t slonidx, size_t altidx ) const;
		bool                                            Transmission_Interpolate            ( double cossza, double alt, double slon, double& transmission ) const;
		bool                                            Deflection_Interpolate              ( double cossza, double alt, double slon, double& deflection) const;

		std::unique_ptr<SKTRAN_RayOptical_Base>			CreateRayAndCalcFullTransmission    ( const HELIODETIC_VECTOR& loc ) const;
		bool 											PrefillTable();

		bool 											InterpolateAndFillSLON				( int slonidx, std::vector<double>& altitude, std::vector<double>& cossza, std::vector<double>& opticaldepth, std::vector<double>& deflection);


	protected:
		virtual bool									FillTable_ClassSpecific				( ) override;
		double                                          CosAngleToSource                    ( const HELIODETIC_UNITVECTOR& look, const HELIODETIC_POINT* location) const override;

		using SKTRAN_SolarTransmission_2D::SetGeometry;	// We Don't want to use this function externally 

	public:
		SKTRAN_SolarTransmission_3D(bool prefill = false, bool accountforrefraction = false);
		virtual ~SKTRAN_SolarTransmission_3D();

		virtual bool			TransmissionAtPoint		(	const HELIODETIC_POINT&					point,
															double&									transmission ) const;

		virtual bool			TransmissionAtVector	(	const HELIODETIC_VECTOR&				vec,
															double&									transmission) const;
		
		virtual bool			DeflectionAtPoint	   (	const HELIODETIC_POINT&				    vec,
															double&									transmission) const;


		virtual bool			InitializeGeometry		( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>&	coords,
														  const SKTRAN_GridDefRayTracingShells_V21&			altgrid,
														  const SKTRAN_GridDefCosSZA_V21&					cosszagrid,
														  const SKTRAN_GridDefSLON_V21&						slongrid );

		void 					SetPrefillGrid			( const SKTRAN_GridDefCosSZA_V21& cosszagrid ) { m_prefillgrid = cosszagrid; }

		void 					DumpTable();
};
