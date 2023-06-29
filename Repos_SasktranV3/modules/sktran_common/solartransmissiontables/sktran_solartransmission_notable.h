

/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_NoTable		2013-06-12*/
/** @ingroup solartransmission
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_SolarTransmission_NoTable : public SKTRAN_SolarTransmission_Base
{

	private:
		mutable int m_numCalls;

    protected:
		virtual bool                FillTable_ClassSpecific         () override;

	public:
									SKTRAN_SolarTransmission_NoTable ();
		virtual					   ~SKTRAN_SolarTransmission_NoTable ();

		void						ReleaseResources               ( );

		virtual bool				MakeThreadSafeFor		       ( size_t                                   numThreads ) override;

		virtual bool				TransmissionAtPoint            ( const HELIODETIC_POINT&                  point,
                                                                     double&                                  transmission) const override;

		virtual bool				TransmissionAtVector           ( const HELIODETIC_VECTOR&                 vec, 
                                                                     double&                                  transmission ) const override;

		virtual bool				SourceTermAtPoint		       ( const SKTRAN_SourceTermQueryObject_Base& qobj, 
                                                                     double*                                  source )    const override;

//		virtual bool				GroundSourceAtPoint            ( const SKTRAN_SourceTermQueryObject_Base& qobj, 
//                                                                     double*                                  source )    const override;

        virtual bool                MonteCarlo_SingleScatteredRadianceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, 
																	 double&								  radiance ) const override;

        virtual bool                MonteCarlo_GroundScatteredRadianceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj,
																	 double&								  radiance)  const override;


		virtual bool				TransmissionAtPoint            ( const double&							  wavelength,
																	 const HELIODETIC_POINT&                  point,
                                                                     double&                                  transmission) const override;

		virtual bool				TransmissionAtVector           ( const double&							  wavelength,
																	 const HELIODETIC_VECTOR&                 vec,
                                                                     double&                                  transmission ) const override;

		virtual bool				SourceTermAtPoint		       ( const double&							  wavelength,
																	 const SKTRAN_SourceTermQueryObject_Base& qobj,
                                                                     double*                                  source )    const override;

//		virtual bool				GroundSourceAtPoint            ( const SKTRAN_SourceTermQueryObject_Base& qobj, 
//                                                                     double*                                  source )    const override;

        virtual bool                MonteCarlo_SingleScatteredRadianceAtPoint ( const double&							 wavelength,
																			    const SKTRAN_SourceTermQueryObject_Base& qobj,
																				double&								     radiance ) const override;

        virtual bool                MonteCarlo_GroundScatteredRadianceAtPoint ( const double&							 wavelength,
																				const SKTRAN_SourceTermQueryObject_Base& qobj,
																				double&								     radiance)  const override;
};




/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_NoTable_reuseRays		2014-2-7*/
/** @ingroup solartransmission
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_SolarTransmission_NoTable_reuseRays : public SKTRAN_SolarTransmission_NoTable
{

	protected:
//		HELIODETIC_UNITVECTOR												m_sununit;
		nxThreadStorageMap<std::unique_ptr<SKTRAN_RayOptical_Base> >		m_rayopt1;
		size_t																m_numThreadsSafeFor;

	protected:
		bool						Allocate(size_t numThreads);
		void                        DeallocateInternalRays( );

	public:
		void						ReleaseResources					();
		bool						InitializeThreadSafeStorageEntry	( std::unique_ptr<SKTRAN_RayOptical_Base>* ptr ) const;

	public:
									SKTRAN_SolarTransmission_NoTable_reuseRays();
								   ~SKTRAN_SolarTransmission_NoTable_reuseRays();
		
		virtual bool				MakeThreadSafeFor		(   size_t										numThreads)			override;

		virtual bool				TransmissionAtPoint		(	const HELIODETIC_POINT&						point, 
																double&										transmission )		const override;

		virtual bool				TransmissionAtVector	(	const HELIODETIC_VECTOR&					vec, 
																double&										transmission ) const override;

		virtual bool				TransmissionAtPoint		(	const double&								wavelength,
																const HELIODETIC_POINT&						point,
																double&										transmission)		const override;

		virtual bool				TransmissionAtVector	(	const double&								wavelength,
																const HELIODETIC_VECTOR&					vec,
																double&										transmission) const override;

};

class SKTRAN_SolarTransmission_NoTable_NoLOSSource : public SKTRAN_SolarTransmission_NoTable
{
	virtual bool				SourceTermAtPoint(const SKTRAN_SourceTermQueryObject_Base& qobj,
		double*                                  source)    const override;
};