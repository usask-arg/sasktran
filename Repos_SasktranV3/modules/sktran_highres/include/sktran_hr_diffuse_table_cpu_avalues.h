

/*-----------------------------------------------------------------------------
 *				 class Aval_IteratorManager_Base		 2016- 12- 2*/
/** These classes provide SKTRAN_HR_Diffuse_Table_CPU access to the phase function and scattering matrices. 
 *  A class to allow RadStore classes to access Avalues via an iterator (for efficiency). 
 *  This is needed because lifetime management depends on what is stored in the Aval class. 
 **/
/*---------------------------------------------------------------------------*/

template< class AvalType >
class Aval_IteratorManager_Base
{
	public:
									    Aval_IteratorManager_Base			( ) { ;}										// Iterators can only be constructed by Aval classes 
									    Aval_IteratorManager_Base			( const Aval_IteratorManager_Base& ) { ;}
	public: 
		virtual					       ~Aval_IteratorManager_Base			( ) { ;}
		virtual void                    Advance                             ( size_t numToAdvance ) = 0;             // Each Aval class must define its own iterator managers
        virtual SKTRAN_HR_WEIGHT_TYPE   ApplyAndAdvance                     ( SKTRAN_HR_WEIGHT_TYPE radiance ) = 0;
};

 using Aval_ScalarIteratorManager_Base  = Aval_IteratorManager_Base< SKTRAN_HR_WEIGHT_TYPE >;
 using Aval_ScatMatIteratorManager_Base = Aval_IteratorManager_Base< SKTRAN_ScatMat_MIMSNC >;

 class EtaCalculator_Base;

/*-----------------------------------------------------------------------------
 *					class Avals_Base		 2016- 12- 1*/
/** Base class defining internal storage area for the scattering terms used in
 *	diffuse calculations. Two of the derived classes Avals_Scalar and 
 *	Avals_MatrixStore use a global storage area to hold all of the scattering terms
 *	for all rays on all diffuse points.  These two classes are the primary techniques
 *	used in the HR engine and are suitable for standard radiation field calculations
 *
 *	A third class, Avals_MatrixTable, developed by Seth Dueck uses a strategy where 
 *	it allows optical properties to be cached for each diffuse point, so that lookups 
 *	into all but the scatter angle grid can be avoided. The purpose of this third class 
 *	is to deploy a strategy that is suitable for more extreme scenarios where we may
 *	have large numbers of incoming and outgoing rays coupled with polarization calculations.
 *
 *  Store "A-values": the scalar scattering function, scatter matrices, and phase matrices 
 *  Access is threadsafe
 *  By convention, scalar Avals are accessed in "incomingRayIndex-major" order. 
 *	If this changes, it needs to be changed in all the iterator manager clases 
 *	and in the scalar radstore classes
**/
/*---------------------------------------------------------------------------*/

class Avals_Base
{
	protected:
		const SKTRAN_TableOpticalProperties_Base*					m_opticaltable;

	private:
		void ReleaseResources ( );


	public:
																	Avals_Base                              ( );
		virtual													   ~Avals_Base                              ( );
		        void												SetOpticalTable                         ( const SKTRAN_TableOpticalProperties_Base* table );
	public:
		virtual bool												AllocateStorage                         ( size_t incomingcounter, size_t outgoingcounter, size_t scattvalcounter, size_t pointscounter ) = 0;
		virtual bool												ComputeMultipliersAndAdjustScatterArray ( const SKTRAN_HR_Diffuse_Point& point ) = 0;
		virtual bool												CalcScatteringMatrixPoint				( const SKTRAN_HR_Diffuse_Point& point ) = 0;
        virtual void												ApplyValue								( const SKTRAN_HR_Diffuse_Point& point, size_t inidx, size_t outidx, SKTRAN_Stokes_NC* radiance ) const = 0;
		virtual void												ApplyValue								( const SKTRAN_HR_Diffuse_Point& point, size_t inidx, size_t outidx, SKTRAN_HR_WEIGHT_TYPE* radiance) const = 0;
        virtual void												CreateEtaCalculator						( std::unique_ptr< EtaCalculator_Base >& etaCalculator ) const = 0; 
		virtual std::unique_ptr< Aval_ScalarIteratorManager_Base >  ScalarAvalIteratorManager				( const SKTRAN_HR_Diffuse_Point& point ) const = 0;
		virtual std::unique_ptr< Aval_ScatMatIteratorManager_Base > MatrixAvalIteratorManager				( const SKTRAN_HR_Diffuse_Point& point ) const = 0;
};



/*-----------------------------------------------------------------------------
 *					Avals_ScalarStore		 2016- 12- 1*/
/** Internal storage area for the scattering terms used in diffuse calculations.
 *	This matrix-like class stores wavelength dependent, scalar, scattering 
 *	matrix values for all the incoming-outgoing rays for each point in the
 *	entire diffuse field. The class also manages surface scattering BRDF terms 
 *	for ground points.
 *
 *	The scattering values for atmospheric points are nominally organized as a
 *	3-D array of [outidx,inidx,pointindex] although this strategy is modified for
 *	ground points. The function,  point.ScatterPropertyIdx(outrayidx, inrayidx) in
 *	class SKTRAN_HR_Diffuse_Point is used to index entries in the table 
 **/
/*---------------------------------------------------------------------------*/

class Avals_ScalarStore : public Avals_Base
{
    private: 
		class Aval_AScalarIteratorManager : public Aval_ScalarIteratorManager_Base
		{
			private: 
				const SKTRAN_HR_WEIGHT_TYPE*			m_currentAval;
                                                        Aval_AScalarIteratorManager ( ); // Don't want to create an iterator that points to nowhere 

			public:
														Aval_AScalarIteratorManager (const SKTRAN_HR_WEIGHT_TYPE* parent) { m_currentAval = parent;} // Only Avals_ScalarStore may construct this Aval_ScalarIteratorManager
													   ~Aval_AScalarIteratorManager ( ) { ;}
				
                virtual void                            Advance                     ( size_t numToAdvance ) override;
                virtual SKTRAN_HR_WEIGHT_TYPE           ApplyAndAdvance             ( SKTRAN_HR_WEIGHT_TYPE radiance ) override;

                // friend class Avals_ScalarStore; // Only Avals_ScalarStore should be able to construct this Aval_ScalarIteratorManager (constructor should be private), but some people don't like using friend... 
		};

	private:
		std::vector<SKTRAN_HR_WEIGHT_TYPE>								m_Avals;

	private:
		bool															CalcScatteringMatrixPoint_Participating ( const SKTRAN_HR_Diffuse_Point& point );
		bool															CalcScatteringMatrixPoint_Boundary      ( const SKTRAN_HR_Diffuse_Point& point );
		void															GetValue								( const SKTRAN_HR_Diffuse_Point& point, size_t inidx, size_t outidx, SKTRAN_ScatMat_MIMSNC* aval_matrix ) const;

	public:
		virtual bool													AllocateStorage                         ( size_t incomingcounter, size_t outgoingcounter, size_t scattvalcounter, size_t pointscounter ) override;
		virtual bool													ComputeMultipliersAndAdjustScatterArray ( const SKTRAN_HR_Diffuse_Point& point ) override;
		virtual bool													CalcScatteringMatrixPoint               ( const SKTRAN_HR_Diffuse_Point& point ) override;
        virtual void													ApplyValue								( const SKTRAN_HR_Diffuse_Point& point, size_t inidx, size_t outidx, SKTRAN_Stokes_NC* radiance ) const override;
		virtual void													ApplyValue								( const SKTRAN_HR_Diffuse_Point& point, size_t inidx, size_t outidx, SKTRAN_HR_WEIGHT_TYPE* radiance ) const override;

        virtual void													CreateEtaCalculator						( std::unique_ptr< EtaCalculator_Base >& etaCalculator ) const override { etaCalculator.reset( new EtaCalculator_DoRotation ); } 
		virtual std::unique_ptr< Aval_ScalarIteratorManager_Base >		ScalarAvalIteratorManager				( const SKTRAN_HR_Diffuse_Point& point ) const override;
		virtual std::unique_ptr< Aval_ScatMatIteratorManager_Base >		MatrixAvalIteratorManager				( const SKTRAN_HR_Diffuse_Point& point ) const override { return std::unique_ptr<Aval_ScatMatIteratorManager_Base>(nullptr);} 
};


/*-----------------------------------------------------------------------------
 *					Avals_MatrixStore		 2016- 12- 1*/
/** **/
/*---------------------------------------------------------------------------*/

template< class MatrixType >
class Avals_MatrixStore : public Avals_Base
{
    public: 
    class Aval_BScalarIteratorManager : public Aval_ScalarIteratorManager_Base
    {
        private: 
            std::vector< SKTRAN_HR_WEIGHT_TYPE > m_vals;
            std::vector< SKTRAN_HR_WEIGHT_TYPE >::iterator m_valIterator;

            Aval_BScalarIteratorManager ( ) { ;} // Only Avals_MatrixStore < MatrixType > may construct this Aval_ScalarIteratorManager
            Aval_BScalarIteratorManager ( const Aval_BScalarIteratorManager& other ) { ;};

        public: 
            ~Aval_BScalarIteratorManager ( ) { ;}
            virtual void Advance ( size_t numToAdvance ) override { m_valIterator += numToAdvance; }
            virtual SKTRAN_HR_WEIGHT_TYPE ApplyAndAdvance ( SKTRAN_HR_WEIGHT_TYPE radiance ) override { return radiance*(*m_valIterator++);}

        friend class Avals_MatrixStore < MatrixType >; // Only Avals_MatrixStore < MatrixType > may construct this Aval_ScalarIteratorManager
    };

    private:
		std::vector< MatrixType >								m_Avals;

    private:
		bool													CalcScatteringMatrixPoint_Participating ( const SKTRAN_HR_Diffuse_Point& point );
		bool													CalcScatteringMatrixPoint_Participating_Impl( const SKTRAN_HR_Diffuse_Point& point, const nxVector& indir, const nxVector& outdir, double cubatureWeight, MatrixType& p, std::unique_ptr< EtaCalculator_Base >& etaCalculator ); 
		bool													CalcScatteringMatrixPoint_Boundary      ( const SKTRAN_HR_Diffuse_Point& point );
		const MatrixType&										GetValueRef								( const SKTRAN_HR_Diffuse_Point& point, size_t inidx, size_t outidx ) const; // {return m_Avals[ point.ScatterPropertyIdx(outidx, inidx) ]; }
		virtual void											CreateEtaCalculatorForStoragePhase		( std::unique_ptr< EtaCalculator_Base >& etaCalculator ) const;

    public:
																Avals_MatrixStore                       ( ) {;}
    virtual													   ~Avals_MatrixStore                       ( ) {;}
    virtual bool												AllocateStorage                         ( size_t incomingcounter, size_t outgoingcounter, size_t scattvalcounter, size_t pointscounter ) override;
    virtual bool												ComputeMultipliersAndAdjustScatterArray ( const SKTRAN_HR_Diffuse_Point& point ) override;
    virtual bool												CalcScatteringMatrixPoint               ( const SKTRAN_HR_Diffuse_Point& point ) override;
    virtual void												ApplyValue								( const SKTRAN_HR_Diffuse_Point& point, size_t inidx, size_t outidx, SKTRAN_Stokes_NC* radiance      ) const override;
	virtual void												ApplyValue								( const SKTRAN_HR_Diffuse_Point& point, size_t inidx, size_t outidx, SKTRAN_HR_WEIGHT_TYPE* radiance ) const override;
	virtual void												CreateEtaCalculator						( std::unique_ptr< EtaCalculator_Base >& etaCalculator ) const override;
    virtual std::unique_ptr< Aval_ScalarIteratorManager_Base >  ScalarAvalIteratorManager				( const SKTRAN_HR_Diffuse_Point& point ) const override;
	virtual std::unique_ptr< Aval_ScatMatIteratorManager_Base > MatrixAvalIteratorManager				( const SKTRAN_HR_Diffuse_Point& point ) const override { return std::unique_ptr<Aval_ScatMatIteratorManager_Base>(nullptr);} 
};



/*-----------------------------------------------------------------------------
 *					Avals_MatrixTable		 2016- 12- 1*/
/** **/
/*---------------------------------------------------------------------------*/

class Avals_MatrixTable : public Avals_Base
{

    private: 
		class Aval_CScalarIteratorManager : public Aval_ScalarIteratorManager_Base
		{
			private: 
				std::vector< SKTRAN_HR_WEIGHT_TYPE >				m_vals;
                
																	Aval_CScalarIteratorManager			(size_t numincoming, size_t numoutgoing) {NXTRACE_ONCEONLY(firsta,("**** 2017-07-25 ***** Aval_CScalarIteratorManager, The implementation of this thing is *quite* inefficient -- I haven't actually investigated how bad it is, though. \n"));  m_vals.resize( numincoming*numoutgoing) ;}
																   

			public:
                                                                   ~Aval_CScalarIteratorManager			( ){ ;}
				//std::vector< SKTRAN_HR_WEIGHT_TYPE >::iterator		begin								( ) {return m_vals.begin();}
                virtual void                                        Advance                             ( size_t numToAdvance ) override { ;}
                virtual SKTRAN_HR_WEIGHT_TYPE                       ApplyAndAdvance                     ( SKTRAN_HR_WEIGHT_TYPE radiance ) override { return 0.0;} 

                friend class Avals_MatrixTable; // Only Avals_MatrixTable may construct this Aval_ScalarIteratorManager 
    };

    private:
		std::vector<double>										m_incomingWeightMultipliers;
		std::vector<SKTRAN_TableOpticalProperties_PointCache>	m_interpFunctions;

    private:
		bool													CalcScatteringMatrixPoint_Participating ( const SKTRAN_HR_Diffuse_Point& point );
		bool													CalcScatteringMatrixPoint_Boundary      ( const SKTRAN_HR_Diffuse_Point& point );
		void													GetValue								( const SKTRAN_HR_Diffuse_Point& point, size_t inidx, size_t outidx, SKTRAN_ScatMat_MIMSNC* Aval ) const;

	public:
																Avals_MatrixTable                       ( ) {;}
    virtual													   ~Avals_MatrixTable                       ( ) {;}
    virtual bool												AllocateStorage                         ( size_t incomingcounter, size_t outgoingcounter, size_t scattvalcounter, size_t pointscounter ) override;
    virtual bool												ComputeMultipliersAndAdjustScatterArray ( const SKTRAN_HR_Diffuse_Point& point ) override;
    virtual bool												CalcScatteringMatrixPoint               ( const SKTRAN_HR_Diffuse_Point& point ) override;
	virtual void												ApplyValue								( const SKTRAN_HR_Diffuse_Point& point, size_t inidx, size_t outidx, SKTRAN_Stokes_NC*      radiance ) const override;
	virtual void												ApplyValue								( const SKTRAN_HR_Diffuse_Point& point, size_t inidx, size_t outidx, SKTRAN_HR_WEIGHT_TYPE* radiance ) const override;
    virtual void												CreateEtaCalculator						( std::unique_ptr< EtaCalculator_Base >& etaCalculator ) const { etaCalculator= std::unique_ptr< EtaCalculator_Base >( new EtaCalculator_DoRotation ); }
    virtual std::unique_ptr< Aval_ScalarIteratorManager_Base >  ScalarAvalIteratorManager				( const SKTRAN_HR_Diffuse_Point& point ) const override;
	virtual std::unique_ptr< Aval_ScatMatIteratorManager_Base > MatrixAvalIteratorManager				( const SKTRAN_HR_Diffuse_Point& point ) const override { return std::unique_ptr<Aval_ScatMatIteratorManager_Base>(nullptr);} 
};


