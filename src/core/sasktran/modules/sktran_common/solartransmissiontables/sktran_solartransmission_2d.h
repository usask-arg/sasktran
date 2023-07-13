


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_2D		2014-2-7*/
/** @ingroup solartransmission
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_SolarTransmission_2D : public SKTRAN_SolarTransmission_Base
{
	private:
		nx2dArray<SKTRAN_StokesScalar>							m_transmission;
		std::shared_ptr<const SKTRAN_CoordinateTransform_V2>	m_coords;

	private:
		std::vector<SKTRAN_Distance>							m_heights;
		SKTRAN_GridDefRayTracingShells_V21						m_heightgrid;
		SKTRAN_GridDefCosSZA_V21								m_cosszagrid;
		size_t													m_numThreadsOnFill;

	protected:
		const SKTRAN_CoordinateTransform_V2*					CoordinatesPtr							() const { return m_coords.get();}
		void													ReleaseResources						();
		const SKTRAN_GridDefCosSZA_V21&							CosSZAGrid								() const { return m_cosszagrid;}
		const SKTRAN_GridDefRayTracingShells_V21&				HeightGrid								() const { return m_heightgrid;}

		virtual bool											AltWeightsForProfile					( double alt, double* altweights, size_t* altindex, size_t& numindex ) const;
		virtual bool											CosSzaWeights							( double cossza, double* szaweights, size_t* szaindex, size_t& numindex ) const;
		virtual bool											CreateRayAndCalcTransmission			( const HELIODETIC_VECTOR& loc, double& transmission ) const;
		virtual bool											FillTable_ClassSpecific					();

	public:
								SKTRAN_SolarTransmission_2D	  ();
		virtual				   ~SKTRAN_SolarTransmission_2D	  ();
		virtual bool			SetGeometry                   ( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>&		coords,
																const SKTRAN_GridDefRayTracingShells_V21&			heightgrid,
                                                                const SKTRAN_GridDefCosSZA_V21&						cosszagrid );

		virtual bool			MakeThreadSafeFor             ( size_t numthreads ) override;

		virtual bool			TransmissionAtPoint           ( const HELIODETIC_POINT&              point,
															    double&                              transmission ) const;

		virtual bool			TransmissionAtVector          ( const HELIODETIC_VECTOR&             vec,
																double&                              transmission ) const;

		virtual bool			SourceTermAtPoint		      ( const SKTRAN_SourceTermQueryObject_Base& qobj,
															    double*                              source )    const override;

//		virtual bool			GroundSourceAtPoint	          ( const SKTRAN_SourceTermQueryObject_Base& qobj, 
//															    double*                              source )    const override;

		virtual bool            MonteCarlo_SingleScatteredRadianceAtPoint( const SKTRAN_SourceTermQueryObject_Base& qobj, 
															    double&                               radiance ) const override;

		virtual bool            MonteCarlo_GroundScatteredRadianceAtPoint( const SKTRAN_SourceTermQueryObject_Base& qobj,
                                                                double&                               radiance ) const override;

};
