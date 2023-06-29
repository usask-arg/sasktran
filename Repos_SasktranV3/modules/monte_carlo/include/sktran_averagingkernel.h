#include "sktran_montecarlo_internals.h"

class SKTRAN_PhotonLog_Base
{
	public:
		virtual       ~SKTRAN_PhotonLog_Base       ( ) {;}
		virtual bool   ConfigureObserverGeometry   ( const SKTRAN_RayOptical_Base& observerRay, double earthRadius ) = 0;
	    virtual bool   ConfigureKernel             ( double minAlt, size_t numAlts, double deltaAlt, double minCosAngle, size_t numCosAngles, double deltaCosAngle, size_t numThreads ) = 0;
		virtual void   WipeKernel                  ( ) = 0;
		virtual bool   AddToKernel                 ( const SKTRAN_MCPhoton_Base* photon, size_t order, size_t threadid) = 0;
		virtual bool   PrintKernel                 ( std::string filenameNoExtension ) = 0;

};

class SKTRAN_PhotonLog_Null : public SKTRAN_PhotonLog_Base
{
	public:
		virtual       ~SKTRAN_PhotonLog_Null       ( ) {;}

        virtual bool   ConfigureObserverGeometry   ( const SKTRAN_RayOptical_Base& observerRay, double earthRadius ) override { return true;}
		virtual bool   ConfigureKernel             ( double minAlt, size_t numAlts, double deltaAlt, double minCosAngle, size_t numCosAngles, double deltaCosAngle, size_t numThreads ) override { return true;}
		virtual void   WipeKernel                  ( ) override { return;}
		virtual bool   AddToKernel                 ( const SKTRAN_MCPhoton_Base* photon, size_t order, size_t threadid) override { return true;}
		virtual bool   PrintKernel                 ( std::string filenameNoExtension ) override { return true;}

};

class SKTRAN_PhotonLog_AveKernel : public SKTRAN_PhotonLog_Base
{

	private:
		std::vector< std::vector<double> > m_grid;
		double m_earthRadius;
		double m_minAlt, m_deltaAlt, m_minCosAngle, m_deltaCosAngle;
		size_t m_numAlts, m_numCosAngles;
		double m_altTolerance;	// Altitudes beyond min/max round to min/max if within this range
		HELIODETIC_VECTOR m_lookUnit, m_upUnit, m_rightUnit;

	private:
		bool           FindGridWeights            ( HELIODETIC_VECTOR vec, bool groundScatter, size_t* indices, double* weights );

	public:
                       SKTRAN_PhotonLog_AveKernel ( );
        virtual       ~SKTRAN_PhotonLog_AveKernel ( );
		virtual bool   ConfigureObserverGeometry  ( const SKTRAN_RayOptical_Base& observerRay, double earthRadius ) override;
	    virtual bool   ConfigureKernel            ( double minAlt, size_t numAlts, double deltaAlt, double minCosAngle, size_t numCosAngles, double deltaCosAngle, size_t numThreads ) override;
		virtual void   WipeKernel                 ( ) override;
		virtual bool   AddToKernel                 ( const SKTRAN_MCPhoton_Base* photon, size_t order, size_t threadid) override;
		virtual bool   PrintKernel                ( std::string filenameNoExtension ) override;

};

class SKTRAN_PhotonLog_RadianceOnLos : public SKTRAN_PhotonLog_Base
{
	protected:
		HELIODETIC_VECTOR m_observer;
		std::vector<double> m_quadPtDistances;
		std::vector< std::vector< double  > >          m_weights;
		std::vector< std::vector< SKTRAN_Stokes_NC > > m_vals;
		std::vector< size_t >                          m_next_lowerIndex;
		std::vector< double >                          m_next_lowerWeight;
		std::vector< double >                          m_next_upperWeight;

	public:
		virtual       ~SKTRAN_PhotonLog_RadianceOnLos       ( ) {;}
		virtual bool   ConfigureObserverGeometry   ( const SKTRAN_RayOptical_Base& observerRay, double earthRadius ) override;
	    virtual bool   ConfigureKernel             ( double minAlt, size_t numAlts, double deltaAlt, double minCosAngle, size_t numCosAngles, double deltaCosAngle, size_t numThreads ) override;
		virtual void   WipeKernel                  ( ) override;
		virtual bool   AddToKernel                 ( const SKTRAN_MCPhoton_Base* photon, size_t order, size_t threadid) override;
		virtual bool   PrintKernel                 ( std::string filenameNoExtension ) override;

};


class SKTRAN_PhotonLog_PhotonsOnLos : public SKTRAN_PhotonLog_Base
{
	protected:
		HELIODETIC_VECTOR m_observer;
		std::vector<double> m_quadPtDistances;
		std::vector< std::vector< std::vector< double > > >          m_weights;
		std::vector< std::vector< std::vector< SKTRAN_Stokes_NC > > > m_vals; // [threadidx] [qptidx] [photonidx]
		std::vector< size_t >                          m_next_lowerIndex;
		std::vector< double >                          m_next_lowerWeight;
		std::vector< double >                          m_next_upperWeight;

	public:
		virtual       ~SKTRAN_PhotonLog_PhotonsOnLos       ( ) {;}
		virtual bool   ConfigureObserverGeometry   ( const SKTRAN_RayOptical_Base& observerRay, double earthRadius ) override;
	    virtual bool   ConfigureKernel             ( double minAlt, size_t numAlts, double deltaAlt, double minCosAngle, size_t numCosAngles, double deltaCosAngle, size_t numThreads ) override;
		virtual void   WipeKernel                  ( ) override;
		virtual bool   AddToKernel                 ( const SKTRAN_MCPhoton_Base* photon, size_t order, size_t threadid) override;
		virtual bool   PrintKernel                 ( std::string filenameNoExtension ) override;

};



class SKTRAN_PhotonLog_ScatterPtOnLos : public SKTRAN_PhotonLog_Base
{
	protected:
		std::vector< std::vector < HELIODETIC_UNITVECTOR > >  m_units;
		std::vector< std::vector < SKTRAN_Stokes_NC > >       m_incoming;
		std::vector< std::vector < SKTRAN_Stokes_NC > >       m_scattered;


	public:
		virtual       ~SKTRAN_PhotonLog_ScatterPtOnLos       ( ) {;}
		virtual bool   ConfigureObserverGeometry   ( const SKTRAN_RayOptical_Base& observerRay, double earthRadius ) override;
	    virtual bool   ConfigureKernel             ( double minAlt, size_t numAlts, double deltaAlt, double minCosAngle, size_t numCosAngles, double deltaCosAngle, size_t numThreads ) override;
		virtual void   WipeKernel                  ( ) override;
		virtual bool   AddToKernel                 ( const SKTRAN_MCPhoton_Base* photon, size_t order, size_t threadid) override;
		virtual bool   PrintKernel                 ( std::string filenameNoExtension ) override;

};


class SKTRAN_PhotonLog_StDev : public SKTRAN_PhotonLog_Base
{

	private:
		std::vector<double> m_sumRad_1fo;
		std::vector<double> m_sumRad_1so;
		std::vector<double> m_sumRad_1ho;
		std::vector<double> m_sumRad_2fo;
		std::vector<double> m_sumRad_2so;
		std::vector<double> m_sumRad_2ho; // Sum of radiance measured at order >=3, squared
		std::vector<double> m_sumRad_2ho_intermediate; // Hold intermediate state so we get just variance in ho scatter
		std::vector<int>    m_numPhotons;
		size_t m_historyInterval;
		size_t m_numIntervals;
		std::vector<size_t>                 m_histidx;
		std::vector< std::vector <double> > m_history_fo;
		std::vector< std::vector <double> > m_history_so;
		std::vector< std::vector <double> > m_history_ho;

	private:
		void CheckForHistoryPush( size_t threadid );

	public:
                       SKTRAN_PhotonLog_StDev    ( );
        virtual       ~SKTRAN_PhotonLog_StDev    ( );
		virtual bool   ConfigureObserverGeometry ( const SKTRAN_RayOptical_Base& observerRay, double earthRadius ) override;
	    virtual bool   ConfigureKernel           ( double minAlt, size_t numAlts, double deltaAlt, double minCosAngle, size_t numCosAngles, double deltaCosAngle, size_t numThreads ) override;
		virtual void   WipeKernel                ( ) override;
		virtual bool   AddToKernel                 ( const SKTRAN_MCPhoton_Base* photon, size_t order, size_t threadid) override;
		virtual bool   PrintKernel               ( std::string filenameNoExtension ) override;

		bool           ConfigureHistoryIntervals ( size_t interval, size_t numIntervals );
};

class SKTRAN_MCVarianceLogger
{
private:
	std::vector<double> m_variance;
	std::vector<size_t> m_numSamples;

public:
	void   SetNumThreads(size_t numThreads) { m_variance.resize(numThreads); m_numSamples.resize(numThreads); }
	void   UpdateThreadVariance(double variance, size_t numSamples, size_t threadid) { m_variance[threadid] = variance; m_numSamples[threadid] = numSamples; }

	double GetTotalVariance() const {
		double var = 0.0;
		double totalNum = 0.0;
		for (int tidx = 0; tidx < (int)m_variance.size(); ++tidx) {
			//var += 0<m_numSamples[tidx] ? m_variance[tidx]/((double)m_numSamples[tidx]) : 0.0; 
			var += m_variance[tidx] * ((double)m_numSamples[tidx])*((double)m_numSamples[tidx]);
			totalNum += (double)m_numSamples[tidx];
		}
		var = 1 < totalNum ? var / ((totalNum - 1)*(totalNum - 1)) : 1.0;
		return var;
	}
};


#define MC_NUMDISTINCTORDERS ((size_t)8)

class SKTRAN_OptimalScatterSequenceManager_Base;

class SKTRAN_MCThreadRadianceLogger
{
	private:
		class orderInfo
		{
			public:
				SKTRAN_Stokes_NC hist1;
				double hist2;
				double variance;
				SKTRAN_Stokes_NC current;
				size_t numSamples;

				void SetToZero ( );
				void UpdateVar ( );
				void AddToMe   ( const orderInfo& other );
		};

	public:
		class runningSums
		{
			public:
				std::vector<SKTRAN_Stokes_NC> m_temprad;
				std::vector<double>           m_cov;
				std::vector<orderInfo>        m_oinfo;

				// new
				std::vector<SKTRAN_Stokes_NC> radBuffer;
				std::vector<SKTRAN_Stokes_NC> radSum;
				std::vector<double> rad2SumVar;
				std::vector<double> rad2SumCov;
				std::vector<size_t> numSamples;
				std::vector<double> varEstimate;
				std::vector<double> covEstimate;
				std::vector<double> varDerivative;
				std::vector<double> covDerivative;
				std::vector<double> varContribution; // contribution of each scatter sequence (including its subsequences) to the variance - probability of the sequence being chosen is proportional to this
				bool minSamplesComplete; // flag indicating when minimum samples are collected and optimal sample selection can begin
				bool primary; // flag indicating if the wavelength associated with this runningSums class is the primary wavelength
				double wavelength;
		};

	private:
		//static const size_t m_numdistinctorders=20;
		size_t           m_oinfoMaxIndex;
		size_t           m_minNumSamplesHigherOrder;
		size_t           m_hardMaxOrder;
		std::vector<double>           m_minFracHO;
		SKTRAN_MCThreadRadianceLogger* m_parent;
		size_t m_primaryWavelengthIndex;
		std::vector<runningSums> m_runningSums; // one per wavelength
		const SKTRAN_OptimalScatterSequenceManager_Base* m_scatterManager;

		size_t m_chunkSize;
		size_t m_chunkCounter;

	private:
		void Merge               ( const SKTRAN_MCThreadRadianceLogger& other ); // Merge #other's data into mine
		size_t FindIdx           ( size_t order ) const;

	public:
		SKTRAN_MCThreadRadianceLogger ( );
		SKTRAN_MCThreadRadianceLogger ( SKTRAN_MCThreadRadianceLogger& other );
		~SKTRAN_MCThreadRadianceLogger( );

		bool SecondaryMeasurementExists() const; 

		bool ConfigureWavelengths( const std::vector<double>& wavelengths, double currentWavelength, size_t primaryWavelengthIndex );
		bool ConfigureOptimalScatterSequence(const SKTRAN_OptimalScatterSequenceManager_Base* scatterManager);

		bool SetChunkSize(size_t chunkSize) { m_chunkSize = chunkSize; m_chunkCounter = chunkSize; return true; }
		bool SetMinFractionHigherOrder(double minfrac);
		bool SetMinFractionHigherOrder(const std::vector<double>& minfrac);
		void ResetLog          ( );
		bool Submit            ( size_t ossIdx, size_t order, const SKTRAN_MCPhoton_Base* photon );
		bool DeclareRayDone    ( size_t scatterSequenceIndex, SKTRAN_MCVarianceLogger& varianceLogger, size_t threadid);
		bool DiscardRay		   ( );
		void SetMaxOrder       ( size_t hardmax );
		//void SetMaxOrder       ( const std::vector<size_t>& hardmax );
		size_t GetNumDistinctOrders() const { return MC_NUMDISTINCTORDERS; }

		bool OptimalScatterSequenceIndex(double randNum, size_t& ossIdx, size_t& maxOrder);
		SKTRAN_Stokes_NC TotalMeasurement() const;
		SKTRAN_Stokes_NC TotalMeasurement(size_t wlidx) const;
		double TargetMeasurement() const;
		double TargetMeasurement(size_t wlidx) const;
		SKTRAN_Stokes_NC MeasurementAtScatterIndex(size_t scidx) const;
		SKTRAN_Stokes_NC MeasurementAtScatterIndex(size_t scidx, size_t wlidx) const;

		bool ExportStatistics(size_t losIdx, double wl);

		double Variance() const;
		double Variance(size_t wlidx) const;

		double SecondaryMeasurement() const;
		double SecondaryMeasurement(size_t wlidx) const;
		double SecondaryVariance() const;
		double SecondaryVariance(size_t wlidx) const;

		const runningSums& RunningSums() const { return m_runningSums[m_primaryWavelengthIndex]; }
		const runningSums& RunningSums(size_t wlidx) const { return m_runningSums[wlidx]; }

		// These should only be used if you kind of know how the class works
		double           rad(size_t order) const { size_t idx = FindIdx(order); return m_runningSums[m_primaryWavelengthIndex].radSum[idx].I() / m_runningSums[m_primaryWavelengthIndex].numSamples[idx]; }
		double           var(size_t order) const { return m_runningSums[m_primaryWavelengthIndex].varEstimate[FindIdx(order)]; }
		double           M1(size_t order) const { return M1_vec(order).I(); }
		SKTRAN_Stokes_NC M1_vec(size_t order) const { return m_runningSums[m_primaryWavelengthIndex].radSum[FindIdx(order)]; }
		double           M2(size_t order) const { return m_runningSums[m_primaryWavelengthIndex].rad2SumVar[FindIdx(order)]; }
		size_t           n(size_t order) const { return m_runningSums[m_primaryWavelengthIndex].numSamples[FindIdx(order)]; }

};




class SKTRAN_MCAirMassFactorLogger
{
	private:
		SKTRAN_MCAirMassFactorLogger* m_parent;

		// AMF parameters
		size_t				m_numamfcells;
		SKTRAN_Stokes_NC	m_temprad;
		std::vector<double> m_tempwsc; // wsc = weighted slant column = (radiance) * (slant column) = (radiance) * (slant length) * (average number density)
		std::vector<double> m_amfCellColumns;
		// sums needed to compute the mean and the variance:
		size_t				m_amfSumsNumSamples;
		double				m_amfSumsI; // radiance (all orders combined)
		double              m_amfSumsII; // radiance squared
		std::vector<double> m_amfSumsW; // weighted slant lengths = slant lengths times radiance
		std::vector<double> m_amfSumsWI; // weighted slant lengths times radiance
		std::vector<double> m_amfSumsWW; // weighted slant lengths squared

		void Merge( const SKTRAN_MCAirMassFactorLogger& other );

	public:
		SKTRAN_MCAirMassFactorLogger();
		SKTRAN_MCAirMassFactorLogger( SKTRAN_MCAirMassFactorLogger& other );
		~SKTRAN_MCAirMassFactorLogger();

		void Submit( int order, const SKTRAN_MCPhoton_Base* photon);
		void DeclareRayDone();
		void DiscardRay();
		
		void Initialize( size_t numcells, const std::vector<double>& verticalColumns );
		void Initialize( SKTRAN_MCAirMassFactorLogger& other );

		double AirMassFactor( size_t amfidx) const;
		double AirMassFactorVariance( size_t amfidx) const;

};


class SKTRAN_MCConvergenceDetector_Base
{	
	private:
		double m_userDesiredStdev;
		int    m_minNumRaysPerThread;

	public:
		virtual ~SKTRAN_MCConvergenceDetector_Base ( ) {;}

	virtual bool SetNumThreads          ( size_t numthreads ) = 0;
	virtual void UpdateThreadData       ( const SKTRAN_MCThreadRadianceLogger& radlog ) = 0;
	virtual bool IsConvergenceReached   ( ) = 0;
	
	void         SetUserDesiredStdev    ( double stdev ) { m_userDesiredStdev = stdev;}
	double       GetUserDesiredStdev    ( ) const        { return m_userDesiredStdev;}
	void         SetMinNumRaysPerThread ( int minNumRaysPerThread ) { m_minNumRaysPerThread = minNumRaysPerThread;}
	int          GetMinNumRaysPerThread ( ) const        { return m_minNumRaysPerThread;}

};


class SKTRAN_MCConvergenceDetector_ThreadIsolated : public SKTRAN_MCConvergenceDetector_Base
{
	private:
		double m_variance;
		double m_radiance;
		size_t m_numSamples;

	public:
		SKTRAN_MCConvergenceDetector_ThreadIsolated();

		virtual bool SetNumThreads    ( size_t numthreads ) override;
		virtual void UpdateThreadData ( const SKTRAN_MCThreadRadianceLogger& radlog ) override;
		virtual bool IsConvergenceReached ( ) override;

};

//
//class SKTRAN_MCConvergenceDetector_ThreadCrosstalk : public SKTRAN_MCConvergenceDetector_Base
//{
//	private:
//		std::function<void( const SKTRAN_MCThreadRadianceLogger& radlog)> m_trueUpdateThreadData;
//
//		
//	public:
//		SKTRAN_MCConvergenceDetector_ThreadCrosstalk();
//		SKTRAN_MCConvergenceDetector_ThreadCrosstalk( SKTRAN_MCConvergenceDetector_ThreadCrosstalk& parent );
//
//		virtual bool SetNumThreads    ( size_t numthreads ) override;
//		virtual void UpdateThreadData ( const SKTRAN_MCThreadRadianceLogger& radlog ) override;
//		virtual bool IsConvergenceReached ( ) override;
//
//};
