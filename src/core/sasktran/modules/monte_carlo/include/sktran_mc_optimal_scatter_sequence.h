#include "sktran_montecarlo_internals.h"

typedef std::vector<size_t>::iterator size_it;
typedef std::vector<uint64_t>::iterator uint64_it;
typedef std::vector<bool>::iterator bool_it;

class SKTRAN_OptimalScatterSequenceManager_Base
{
	protected:
		size_t m_hardMax;
		size_t m_numDistinctOrders;
		std::string m_filename;
		std::vector<size_t> m_covToLowerVarIdx;
		std::vector<size_t> m_covToUpperVarIdx;

	protected:

		std::string PrintSequence(uint64_t sequenceBits, size_t order, bool elasticRaman) const;
		virtual std::string PrintSequence(size_t varIdx) const = 0;

		virtual bool CalculateVariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& variance) const = 0;
		virtual bool CalculateCovariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& covariance) const = 0;
		virtual bool CalculateDerivative(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& varDerivative, std::vector<double>& covDerivative) const = 0;
		virtual bool CalculateVarianceContribution(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& contribution) const = 0;

	public:
		SKTRAN_OptimalScatterSequenceManager_Base(); 
		virtual ~SKTRAN_OptimalScatterSequenceManager_Base();

		virtual bool SecondaryMeasurement() const { return false; }

		virtual bool SetNumDistinctOrders(size_t numDistinctOrders) { m_numDistinctOrders = numDistinctOrders; return true; }
		virtual bool SetMaxOrder(size_t hardmax);

		virtual bool OptimalScatterSequenceIndex( const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double randNum, size_t& ossIdx, bool& minSamplesComplete ) const = 0; // select an optimal scatter sequence, return the corresponding index
		virtual bool SortSamples( const size_t& ossIdx, SKTRAN_MCThreadRadianceLogger::runningSums& rSums ) const = 0;
		virtual bool SubmitSample(const size_t ossIdx, const size_t order, const SKTRAN_MCPhoton_Base* mcphoton, std::vector<SKTRAN_MCThreadRadianceLogger::runningSums>& rSums) const = 0;

		// ossIdx = optimal scatter sequence index
		// order = the reverse position in this sequence (ie order 3 in sequence CRCC is a Raman scatter)
		virtual bool ElasticScatter( size_t ossIdx, size_t order, bool& elastic ) const = 0;				// if ossIdx -> CRCC, order = 3: returns true because order 3 (counting backwards) is Raman
		virtual bool Order( size_t ossIdx, size_t& order ) const = 0;										// if ossIdx -> CRCC: returns 4

		//virtual bool VarianceIndex( size_t ossIdx, size_t order, size_t& index ) const = 0;						// if ossIdx -> CRCC, order = 3: variance index of RCC
		//virtual bool CovarianceIndex( size_t ossIdx, size_t order1, size_t order2, size_t& index ) const = 0;	// if ossIdx -> CRCC, order1 = 3, order2 = 2: covariance index of RCC and CC

		virtual bool ConfigureRunningSums(SKTRAN_MCThreadRadianceLogger::runningSums& rSums) const = 0;
		virtual bool ProcessVariance(SKTRAN_MCThreadRadianceLogger::runningSums& rSums) const;
		virtual bool CalculateMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, SKTRAN_Stokes_NC& measurement) const = 0;
		virtual bool CalculateTotalVariance(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& variance) const = 0;
		virtual bool CalculateSecondaryMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& measurement) const = 0;
		virtual bool CalculateSecondaryVariance(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& measurement) const = 0;
		virtual bool CalculateTargetMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& variance) const = 0; // measurement used to determine target precision
		virtual bool CalculateTargetVariance(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& variance, size_t& numSamples) const = 0; // variance used to compare with target precision

		virtual size_t NumVarianceTerms( ) const = 0;
		virtual size_t NumCovarianceTerms( ) const = 0;

		virtual bool ConfigureStatisticsExport(const std::string& filename) { m_filename = filename; return true; }
		virtual bool ExportStatistics(const std::vector<SKTRAN_MCThreadRadianceLogger::runningSums>& rSums, const std::string& suffix) const;

};

class SKTRAN_OptimalScatterSequenceManager_Uniform : public SKTRAN_OptimalScatterSequenceManager_Base
{
	protected:
		size_t m_minNumSamplesHigherOrder;
		double m_minFracHO;

	protected:
		virtual bool CalculateVariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& variance) const override;
		virtual bool CalculateCovariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& covariance) const override;
		virtual bool CalculateDerivative(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& varDerivative, std::vector<double>& covDerivative) const override { return true; }
		virtual bool CalculateVarianceContribution(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& contribution) const override;

		virtual std::string PrintSequence(size_t varIdx) const override;

	public:
		SKTRAN_OptimalScatterSequenceManager_Uniform();
		virtual ~SKTRAN_OptimalScatterSequenceManager_Uniform();

		bool SetMinFractionHigherOrder(double minfrac) { m_minFracHO = minfrac; return minfrac >= 0.0 && minfrac <= 1.0; }
		bool SetMinNumSamplesHigherOrder(size_t numSamples) { m_minNumSamplesHigherOrder = numSamples; return true; }
		virtual bool SetMaxOrder(size_t hardmax) override;


		virtual bool OptimalScatterSequenceIndex( const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double randNum, size_t& ossIdx, bool& minSamplesComplete ) const override;
		virtual bool SortSamples( const size_t& ossIdx, SKTRAN_MCThreadRadianceLogger::runningSums& rSums ) const override;
		virtual bool SubmitSample(const size_t ossIdx, const size_t order, const SKTRAN_MCPhoton_Base* mcphoton, std::vector<SKTRAN_MCThreadRadianceLogger::runningSums>& rSums) const override;

		virtual bool ElasticScatter(size_t ossIdx, size_t order, bool& elastic) const override;
		virtual bool Order(size_t ossIdx, size_t& order) const override;

		virtual bool ConfigureRunningSums(SKTRAN_MCThreadRadianceLogger::runningSums& rSums) const override;
		virtual bool CalculateMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, SKTRAN_Stokes_NC& measurement) const override;
		virtual bool CalculateTotalVariance(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& variance) const override;
		virtual bool CalculateSecondaryMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& measurement) const override { return true; };
		virtual bool CalculateSecondaryVariance(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& variance) const override { return true; };
		virtual bool CalculateTargetMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& measurement) const override; // measurement used to determine target precision
		virtual bool CalculateTargetVariance(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& variance, size_t& numSamples) const; // variance used to compare with target precision

		virtual size_t NumVarianceTerms() const override { return m_numDistinctOrders; }
		virtual size_t NumCovarianceTerms() const override { return m_numDistinctOrders * (m_numDistinctOrders - 1) / 2; }


};

class SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic : public SKTRAN_OptimalScatterSequenceManager_Base
{
	protected:
		size_t m_numSequences;
		size_t m_numCovarianceTerms;
		size_t m_numMaxSequences;

		size_t m_minSamplesHigherOrder;
		std::vector<double> m_minFractionHigherOrder;
		double m_minFractionHigherOrderSum;

		std::vector<size_t> m_maxTotalOrders;
		std::vector<size_t> m_maxRamanOrders;

		uint64_t  m_bits[63]; // integers for picking out 1 bit with the & operation
		std::vector<uint64_t> m_seq; // maps seqIdx -> integer with bits representing the sequence (0 <-> elastic, 1 <-> inelastic, MSB <-> 1st order scatter / final physical scatter)
		std::vector<size_t> m_seqToTotalOrder; // maps seqIdx -> number of total scatters
		std::vector<size_t> m_seqToRamanOrder; // maps seqIdx -> number of Raman scatters
		std::vector<size_t> m_seqToSeqPlusCabannes; // maps seqIdx -> seqIdx of original sequence with an additional Cabannes scatter, if it exists
		std::vector<size_t> m_seqToSeqPlusRaman; // maps seqIdx -> seqIdx of original sequence with an additional Raman scatter, if it exists

		std::vector<size_t> m_covToLowerSeqIdx; // maps covIdx -> lower seqIdx
		std::vector<size_t> m_covToUpperSeqIdx; // maps covIdx -> upper seqIdx
		std::vector<size_t> m_covGroupByUpperIdx; // maps covIdx (grouped by upper) -> covIdx (grouped by lower)

		std::vector<size_t> m_seqToCovLowerIdx; // maps seqIdx -> start of group in m_covToLowerSeqIdx (contains all covariance pairs with seqIdx as the lower index)
		std::vector<size_t> m_seqToCovLowerSize; // maps seqIdx -> size of group in m_covToLowerSeqIdx (contains all covariance pairs with seqIdx as the lower index)
		std::vector<size_t> m_seqToCovUpperIdx; // maps seqIdx -> start of group in m_covGroupByUpperIdx (contains all covariance pairs with seqIdx as the upper index)
		std::vector<size_t> m_seqToCovUpperSize; // maps seqIdx -> size of group in m_covGroupByUpperIdx (contains all covariance pairs with seqIdx as the upper index)

		std::vector<size_t> m_maxSeqIdx; // list of "maximum sequence" indices: sequences with maximum total scatters for each order of Raman scatters
		std::vector<size_t> m_maxSeqToMinSamples; // minimum samples for each of the maximum sequences
		std::vector<double> m_maxSeqToMinFraction; // minimum fraction for each of the maximum sequences

	private:
		bool CalculateMaxOrders(const std::vector<size_t>& maxTotalOrders, std::vector<size_t>& maxRamanOrders) const;
		bool CalculateNumberOfTerms(size_t& numSequences, size_t& numCovarianceTerms) const;
		bool DefineSequences(std::vector<uint64_t>& sequences, std::vector<size_t>& seqToTotalOrder, std::vector<size_t>& seqToRamanOrder) const;
		bool DefineCovarianceIndicesLower(std::vector<size_t>& covToLowerSeqIdx, std::vector<size_t>& covToUpperSeqIdx, std::vector<size_t>& seqToCovLowerIdx, std::vector<size_t>& seqToCovLowerSize) const;
		bool DefineCovarianceIndicesUpper(std::vector<size_t>& covGroupByUpperIdx, std::vector<size_t>& seqToCovUpperIdx, std::vector<size_t>& seqToCovUpperSize) const;
		bool DefineSequenceIncrements(std::vector<size_t>& seqToSeqPlusCabannes, std::vector<size_t>& seqToSeqPlusRaman) const;
		bool CalculateMinimumFractions(size_t& numMaxSeq, std::vector<size_t>& maxSeq, std::vector<double>& maxSeqToMinFraction, std::vector<size_t>& maxSeqToMinSamples) const;

		bool SequenceIterator(size_t T, size_t R, size_it& totalOrder, size_it& ramanOrder, uint64_it& sequence) const;
		bool SequenceIterator(size_t T, size_t R, size_t t, size_t r, uint64_t n, size_it& totalOrder, size_it& ramanOrder, uint64_it& sequence) const;

	protected:
		bool Combination(size_t n, size_t r, size_t& c) const;

		virtual bool CalculateVariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& variance) const override;
		virtual bool CalculateCovariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& covariance) const override;
		virtual bool CalculateDerivative(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& varDerivative, std::vector<double>& covDerivative) const override { return true; }
		virtual bool CalculateVarianceContribution(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& contribution) const override;

		virtual std::string PrintSequence(size_t varidx) const override;

	public:
		SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic();
		virtual ~SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic() {}

		virtual bool   SetMaxOrder( size_t hardmax ) override;

		bool		   SetMinFractionHigherOrder(const std::vector<double>& minFrac);
		bool		   SetMaxRamanOrders(const std::vector<size_t>& maxOrders);

		virtual bool OptimalScatterSequenceIndex( const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double randNum, size_t& ossIdx, bool& minSamplesComplete ) const override;
		virtual bool SortSamples( const size_t& ossIdx, SKTRAN_MCThreadRadianceLogger::runningSums& rSums ) const override;
		virtual bool SubmitSample( const size_t ossIdx, const size_t order, const SKTRAN_MCPhoton_Base* mcphoton, std::vector<SKTRAN_MCThreadRadianceLogger::runningSums>& rSums ) const override;

		//bool SetMinNumSamplesHigherOrder(size_t numSamples) { m_minNumSamplesHigherOrder = numSamples; return true; }

		virtual bool ElasticScatter(size_t ossIdx, size_t order, bool& elastic) const override;
		virtual bool Order(size_t ossIdx, size_t& order) const override;

		virtual bool ConfigureRunningSums(SKTRAN_MCThreadRadianceLogger::runningSums& rSums) const;
		virtual bool CalculateMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, SKTRAN_Stokes_NC& measurement) const override;
		virtual bool CalculateTotalVariance(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& variance) const override;
		virtual bool CalculateSecondaryMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& measurement) const override { return true; }
		virtual bool CalculateSecondaryVariance(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& variance) const override { return true; }
		virtual bool CalculateTargetMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& measurement) const override; // measurement used to determine target precision
		virtual bool CalculateTargetVariance(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& variance, size_t& numSamples) const; // variance used to compare with target precision

		virtual size_t NumVarianceTerms() const override { return m_numSequences; }
		virtual size_t NumCovarianceTerms() const override { return m_numCovarianceTerms; }
};

class SKTRAN_OptimalScatterSequenceManager_UniformSecondary : public SKTRAN_OptimalScatterSequenceManager_Uniform
{
	protected:

		virtual bool CalculateVariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& variance) const override;
		virtual bool CalculateCovariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& covariance) const override;

		virtual std::string PrintSequence(size_t varIdx) const override;

	public:
		SKTRAN_OptimalScatterSequenceManager_UniformSecondary();
		virtual ~SKTRAN_OptimalScatterSequenceManager_UniformSecondary() {}

		virtual bool SecondaryMeasurement() const override { return true; }

		bool SetMinFractionHigherOrder(double minfrac) { m_minFracHO = minfrac; return minfrac >= 0.0 && minfrac <= 1.0; }
		bool SetMinNumSamplesHigherOrder(size_t numSamples) { m_minNumSamplesHigherOrder = numSamples; return true; }
		virtual bool SetMaxOrder(size_t hardmax) override;

		virtual bool SortSamples(const size_t& ossIdx, SKTRAN_MCThreadRadianceLogger::runningSums& rSums) const override;
		virtual bool SubmitSample(const size_t ossIdx, const size_t order, const SKTRAN_MCPhoton_Base* mcphoton, std::vector<SKTRAN_MCThreadRadianceLogger::runningSums>& rSums) const override;

		virtual bool ConfigureRunningSums(SKTRAN_MCThreadRadianceLogger::runningSums& rSums) const override;
		virtual bool CalculateSecondaryMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& measurement) const override = 0;
		virtual bool CalculateSecondaryVariance(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& variance) const override = 0;
		virtual bool CalculateTargetMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& measurement) const override; // measurement used to determine target precision
		virtual bool CalculateTargetVariance(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& variance, size_t& numSamples) const; // variance used to compare with target precision

		virtual size_t NumVarianceTerms() const override { return 2 * m_numDistinctOrders; }
		virtual size_t NumCovarianceTerms() const override { return m_numDistinctOrders * (2 * m_numDistinctOrders - 1); }
};

class SKTRAN_OptimalScatterSequenceManager_UniformRing: public SKTRAN_OptimalScatterSequenceManager_UniformSecondary
{
	protected:
		virtual bool CalculateSecondaryMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& measurement) const override;
		virtual bool CalculateSecondaryVariance(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& variance) const override;
};

class SKTRAN_OptimalScatterSequenceManager_UniformFillingIn: public SKTRAN_OptimalScatterSequenceManager_UniformSecondary
{
	protected:
		virtual bool CalculateSecondaryMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& measurement) const override;
		virtual bool CalculateSecondaryVariance(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& variance) const override;
};

class SKTRAN_OptimalScatterSequenceManager_UniformElastic : public SKTRAN_OptimalScatterSequenceManager_UniformSecondary
{
	public:
		virtual bool CalculateSecondaryMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& measurement) const override;
		virtual bool CalculateSecondaryVariance(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& variance) const override;
};

class SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary : public SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic
{
	protected:
		size_t m_numVarianceTerms; // size of m_var, etc.
		size_t m_numSubSequences; // size of m_subSeq, m_subSeqToSubVarIdx
		size_t m_numSubVarianceTerms; // size of m_subVar, m_subVarToBufferIdx

		std::vector<uint64_t> m_var; // similar to m_seq, but has two entries for each scatter sequence with at least one Raman scatter (one with wavelength shift, one without)
		std::vector<size_t> m_varToTotalOrder; // maps varIdx -> number of total scatters
		std::vector<size_t> m_varToRamanOrder; // maps varIdx -> number of Raman scatters
		std::vector<bool> m_varToElasticRaman; // maps varIdx -> true for elastic Raman (neglect wavelength shift), false for inelastic Raman

		std::vector<size_t> m_varToSeqIdx; // maps varIdx -> seqIdx corresponding to the same scatter sequence
		std::vector<size_t> m_seqToVarIdx; // maps seqIdx -> the first varIdx corresponding to the same scatter sequence
		std::vector<size_t> m_seqToVarSize; // maps seqIdx -> the size of the group in varIdx containing the same scatter sequence (should only be 1 or 2)
		
		//std::vector<size_t> m_covToLowerVarIdx; // maps covIdx -> lower varIdx
		//std::vector<size_t> m_covToUpperVarIdx; // maps covIdx -> upper varIdx

		std::vector<size_t> m_varToCovLowerIdx; // maps varIdx -> beginning of group in m_covToLowerVarIdx (contains all covariance pairs with varIdx as the lower index)
		std::vector<size_t> m_varToCovLowerSize; // maps varIdx -> size of group in m_covToLowerVarIdx (contains all covariance pairs with varIdx as the lower index)
		std::vector<size_t> m_varToCovUpperIdx; // maps varIdx -> beginning of group in m_covGroupByUpperIdx (contains all covariance pairs with varIdx as the upper index)
		std::vector<size_t> m_varToCovUpperSize; // maps varIdx -> size of group in m_covGroupByUpperIdx (contains all covariance pairs with varIdx as the upper index)

		std::vector<size_t> m_subSeq; // contains seqIdx for the sub-sequence leading up to every sequence (ex: {C}{R}{C CC}{C CR}{R RC}{R RR}{C CC CCC} etc.)
		std::vector<size_t> m_seqToSubSeqIdx; // maps seqIdx -> start of the sub-sequence group in m_subSeq

		std::vector<size_t> m_subVar; // contains varIdx for the sub-sequence leading up to every sequence (ex: {C}{R R'}{C CC}{C CR CR'}{R R' RC RC'} etc.)
		std::vector<size_t> m_seqToSubVarIdx; // maps seqIdx -> the beginning of the corresponding group in m_subVar
		std::vector<size_t> m_seqToSubVarSize; // maps seqIdx -> the size of the corresponding group in m_subVar
		std::vector<size_t> m_subSeqToSubVarIdx; // maps subSeqIdx -> subVarIdx of the matching sequence in m_subVar
		std::vector<size_t> m_subVarToBufferIdx; // maps subVarIdx -> the index where the samples will be stored in radBuffer

	private:
		bool CalculateNumberOfTerms(size_t& numVarianceTerms, size_t& numCovarianceTerms, size_t& numSubSequences, size_t& numSubVarianceTerms) const;
		bool DefineSequences(std::vector<uint64_t>& sequences, std::vector<size_t>& seqToTotalOrder, std::vector<size_t>& seqToRamanOrder, std::vector<bool>& seqToElasticRaman) const;
		bool MapSequences(std::vector<size_t>& varToSeqIdx, std::vector<size_t>& seqToVarIdx, std::vector<size_t>& seqToVarSize);
		bool DefineCovarianceIndicesLower(std::vector<size_t>& covToLowerSeqIdx, std::vector<size_t>& covToUpperSeqIdx, std::vector<size_t>& seqToCovLowerIdx, std::vector<size_t>& seqToCovLowerSize, std::vector<size_t>& covToLowerVarIdx, std::vector<size_t>& covToUpperVarIdx, std::vector<size_t>& varToCovLowerIdx, std::vector<size_t>& varToCovLowerSize) const;
		bool DefineCovarianceIndicesUpper(std::vector<size_t>& covGroupByUpperIdx, std::vector<size_t>& varToCovUpperIdx, std::vector<size_t>& varToCovUpperSize, std::vector<size_t>& seqToCovUpperIdx, std::vector<size_t>& seqToCovUpperSize) const;
		bool DefineSubSequenceIndices(std::vector<size_t>& subSeq, std::vector<size_t>& seqToSubSeqIdx) const;
		bool DefineSubVarianceIndices(std::vector<size_t>& subVar, std::vector<size_t>& seqToSubVarIdx, std::vector<size_t>& seqToSubVarSize, std::vector<size_t>& subSeqToSubVarIdx, std::vector<size_t>& subVarToBufferIdx) const;
		
		bool SequenceIterator(size_t T, size_t R, size_it& totalOrder, size_it& ramanOrder, bool_it& elasticRaman, uint64_it& sequence) const;
		bool SequenceIterator(size_t T, size_t R, size_t t, size_t r, uint64_t n, size_it& totalOrder, size_it& ramanOrder, bool_it& elasticRaman, uint64_it& sequence) const;

	protected:
		virtual bool CalculateRadiances(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& inelasticRadiance, double& elasticRadiance) const;
		virtual std::string PrintSequence(size_t varIdx) const override; // seq=true for idx=seqIdx, seq=false for idx=varIdx

		virtual bool CalculateVariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& variance) const override;
		virtual bool CalculateCovariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& covariance) const override;
		virtual bool CalculateVarianceContribution(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& contribution) const override;

	public:
		SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary() {}
		virtual ~SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary() {}

		virtual bool SecondaryMeasurement() const override { return true; }

		virtual bool SetMaxOrder(size_t hardmax) override;

		virtual bool OptimalScatterSequenceIndex(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double randNum, size_t& ossIdx, bool& minSamplesComplete) const override;

		virtual bool SortSamples(const size_t& ossIdx, SKTRAN_MCThreadRadianceLogger::runningSums& rSums) const override;
		virtual bool SubmitSample(const size_t ossIdx, const size_t order, const SKTRAN_MCPhoton_Base* mcphoton, std::vector<SKTRAN_MCThreadRadianceLogger::runningSums>& rSums) const override;

		virtual bool ConfigureRunningSums(SKTRAN_MCThreadRadianceLogger::runningSums& rSums) const override;
		virtual bool CalculateMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, SKTRAN_Stokes_NC& measurement) const override;
		virtual bool CalculateTotalVariance(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& variance) const override;
		virtual bool CalculateSecondaryMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& measurement) const override = 0;
		virtual bool CalculateSecondaryVariance(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& variance) const override = 0;
		virtual bool CalculateTargetMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& measurement) const override; // measurement used to determine target precision
		virtual bool CalculateTargetVariance(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& variance, size_t& numSamples) const override; // variance used to compare with target precision

		virtual size_t NumVarianceTerms() const override { return m_numVarianceTerms; }
		virtual size_t NumCovarianceTerms() const override { return m_numCovarianceTerms; }
};

class SKTRAN_OptimalScatterSequenceManager_OptimizedRing : public SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary
{
	protected:
		virtual bool CalculateDerivative(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& varDerivative, std::vector<double>& covDerivative) const override;
		virtual bool CalculateSecondaryMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& measurement) const override;
		virtual bool CalculateSecondaryVariance(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& variance) const override;
};

class SKTRAN_OptimalScatterSequenceManager_OptimizedFillingIn : public SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary
{
	protected:
		virtual bool CalculateDerivative(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& varDerivative, std::vector<double>& covDerivative) const override;
		virtual bool CalculateSecondaryMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& measurement) const override;
		virtual bool CalculateSecondaryVariance(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& variance) const override;
};

class SKTRAN_OptimalScatterSequenceManager_OptimizedElastic: public SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary
{
	protected:
		virtual bool CalculateDerivative(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& varDerivative, std::vector<double>& covDerivative) const override;
		virtual bool CalculateSecondaryMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& measurement) const override;
		virtual bool CalculateSecondaryVariance(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double& variance) const override;
};
