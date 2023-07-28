#include "include/sktran_montecarlo_internals.h"


std::string SKTRAN_OptimalScatterSequenceManager_Base::PrintSequence(uint64_t sequenceBits, size_t order, bool elasticRaman) const
{
	bool ok = true;

	std::string sequenceStr;
	uint64_t selector;
	size_t numC = 0;

	if (elasticRaman) sequenceStr = "*";

	for (size_t i = 0; i < order; i++)
	{
		selector = pow(2, order - 1 - i);
		if (sequenceBits & selector)
		{
			if (numC == 1) sequenceStr += "C";
			if (numC == 2) sequenceStr += "CC";
			if (numC > 2) sequenceStr += std::to_string(numC) + "C";
			sequenceStr += "R";
			numC = 0;
		}
		else
		{
			numC += 1;
		}
	}
	if (numC == 1) sequenceStr += "C";
	if (numC == 2) sequenceStr += "CC";
	if (numC > 2) sequenceStr += std::to_string(numC) + "C";

	return sequenceStr;
}

SKTRAN_OptimalScatterSequenceManager_Base::SKTRAN_OptimalScatterSequenceManager_Base()
{
	m_hardMax = 50;
	m_numDistinctOrders = 50;
}

SKTRAN_OptimalScatterSequenceManager_Base::~SKTRAN_OptimalScatterSequenceManager_Base()
{
}

bool SKTRAN_OptimalScatterSequenceManager_Base::SetMaxOrder(size_t hardmax)
{
	m_hardMax = hardmax; 
	m_numDistinctOrders = std::min(hardmax, m_numDistinctOrders);
	return true;
}

bool SKTRAN_OptimalScatterSequenceManager_Base::ProcessVariance(SKTRAN_MCThreadRadianceLogger::runningSums & rSums) const
{
	bool ok = true;
	ok = ok && CalculateVariance(rSums, rSums.varEstimate);
	ok = ok && CalculateCovariance(rSums, rSums.covEstimate);
	ok = ok && CalculateDerivative(rSums, rSums.varDerivative, rSums.covDerivative);
	if (rSums.primary) ok = ok && CalculateVarianceContribution(rSums, rSums.varContribution);
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_Base::ExportStatistics(const std::vector<SKTRAN_MCThreadRadianceLogger::runningSums>& rSums, const std::string& suffix) const
{
	bool ok = true;
	if (!m_filename.empty())
	{
		FILE* pFile = fopen((m_filename + suffix).c_str(), "w");

		fprintf(pFile, "# variance terms: %6zd | covariance terms: %6zd | wavelengths: %6zd\n", NumVarianceTerms(), NumCovarianceTerms(), rSums.size());

		fprintf(pFile, "\n####### wavelengths\n");
		for (auto&& r : rSums) fprintf(pFile, "%19.12f ", r.wavelength);

		fprintf(pFile, "\n\n# i ###### sequence ########### samples\n");
		for (size_t i = 0; i < NumVarianceTerms(); i++)
		{
			fprintf(pFile, "%3zd %15s ", i, PrintSequence(i).c_str());
			fprintf(pFile, "%19zd ", rSums[0].numSamples[i]);
			fprintf(pFile, "\n");
		}

		fprintf(pFile, "\n# i ###### sequence ########## mean(wl)\n");
		for (size_t i = 0; i < NumVarianceTerms(); i++)
		{
			fprintf(pFile, "%3zd %15s ", i, PrintSequence(i).c_str());
			for (auto&& r : rSums)
			{
				fprintf(pFile, "%19.12e ", r.radSum[i].I() / r.numSamples[i]);
			}
			fprintf(pFile, "\n");
		}
		fprintf(pFile, "\n# i ###### sequence ###### variance(wl)\n");
		for (size_t i = 0; i < NumVarianceTerms(); i++)
		{
			fprintf(pFile, "%3zd %15s ", i, PrintSequence(i).c_str());
			for (auto&& r : rSums)
			{
				fprintf(pFile, "%19.12e ", r.varEstimate[i]);
			}
			fprintf(pFile, "\n");
		}

		fprintf(pFile, "\n# c # l # u ## sequence l ## sequence u #### covariance(wl)\n");
		for (size_t cidx = 0; cidx < NumCovarianceTerms(); cidx++)
		{
			size_t lidx = m_covToLowerVarIdx[cidx];
			size_t uidx = m_covToUpperVarIdx[cidx];
			fprintf(pFile, "%3zd %3zd %3zd %13s %13s ", cidx, lidx, uidx, PrintSequence(lidx).c_str(), PrintSequence(uidx).c_str());
			for (auto&& r : rSums)
			{
				fprintf(pFile, "%19.12e ", r.covEstimate[cidx]);
			}
			fprintf(pFile, "\n");
		}

		fclose(pFile);
	}
	return ok;
}



SKTRAN_OptimalScatterSequenceManager_Uniform::SKTRAN_OptimalScatterSequenceManager_Uniform()
{
	m_minNumSamplesHigherOrder = 100;
}

SKTRAN_OptimalScatterSequenceManager_Uniform::~SKTRAN_OptimalScatterSequenceManager_Uniform()
{
}

bool SKTRAN_OptimalScatterSequenceManager_Uniform::SetMaxOrder(size_t hardmax)
{
	bool ok = true;
	ok = ok && SKTRAN_OptimalScatterSequenceManager_Base::SetMaxOrder(hardmax);

	size_t numCov = m_numDistinctOrders * (m_numDistinctOrders - 1) / 2;
	m_covToLowerVarIdx.resize(numCov);
	m_covToUpperVarIdx.resize(numCov);
	size_t cidx = 0;
	for (size_t lidx = 0; lidx < m_numDistinctOrders; lidx++)
	{
		for (size_t uidx = lidx + 1; uidx < m_numDistinctOrders; uidx++)
		{
			m_covToLowerVarIdx[cidx] = lidx;
			m_covToUpperVarIdx[cidx++] = uidx;
		}
	}
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_Uniform::OptimalScatterSequenceIndex(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double r, size_t& ossIdx, bool& minSamplesComplete ) const
{
	bool ok = true;
	size_t optorder;

	double target;

	if (!minSamplesComplete) // If we haven't samples enough of the ho rays yet
	{
		size_t n = rSums.numSamples.back();
		if (n < m_minNumSamplesHigherOrder)
		{
			optorder = m_hardMax - 1;
		}
		if (++n >= m_minNumSamplesHigherOrder) minSamplesComplete = true;
		ossIdx = optorder;
		return ok;
	}

	if (r < m_minFracHO || m_hardMax == 1)  // If we want to make sure we sample the ho space by at least some min amount
	{
		optorder = m_hardMax - 1; // Scatter to higher order
	}
	else {
		r = (r - m_minFracHO) / (1.0 - m_minFracHO); // Rescale uni(0,1) to account for using it twice
		double numhscattersave = 20.0; // Assuming 20 orders of scatter for h


		{
			const std::vector<double>& pro = rSums.varContribution;
			target = r * *(pro.end() - 1);
			auto it = std::upper_bound(pro.begin(), pro.end(), target);
			if (it < pro.end() - 1) {
				optorder = std::distance(pro.begin(), it) + 1;
			}
			else {
				optorder = m_hardMax - 1;
			}
			//printf("%i\n",optorder);
		}
	}

	ossIdx = optorder;
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_Uniform::SortSamples(const size_t & ossIdx, SKTRAN_MCThreadRadianceLogger::runningSums & rSums) const
{
	bool ok = true;

	size_t maxOrder = ossIdx < m_numDistinctOrders ? ossIdx + 1 : m_numDistinctOrders;
	size_t cidx = 0;
	for (size_t i = 0; i < maxOrder; i++)
	{
		rSums.radSum[i] += rSums.radBuffer[i];
		rSums.rad2SumVar[i] += std::pow(rSums.radBuffer[i].I(), 2.0);
		rSums.numSamples[i]++;
		for (size_t j = i + 1; j < maxOrder; j++)
		{
			rSums.rad2SumCov[cidx++] += rSums.radBuffer[i].I() * rSums.radBuffer[j].I();
		}
		cidx += m_numDistinctOrders - std::max(maxOrder, i + 1);
	}

	for (size_t i = 0; i < m_numDistinctOrders; i++) rSums.radBuffer[i].SetTo(0.0);
	   
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_Uniform::SubmitSample(const size_t ossIdx, const size_t order, const SKTRAN_MCPhoton_Base * mcphoton, std::vector<SKTRAN_MCThreadRadianceLogger::runningSums>& rSums) const
{
	bool ok = true;
	size_t idx = order < m_numDistinctOrders ? order - 1 : m_numDistinctOrders - 1;

	auto pr = mcphoton->photonRadiances().cbegin();
	for (auto rs = rSums.begin(); rs != rSums.end(); rs++, pr++)
	{
		rs->radBuffer[idx] += pr->GetRecentContribVec();
	}
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_Uniform::Order(size_t ossIdx, size_t& order) const
{
	bool ok = true;
#if defined(NXDEBUG)
	ok = ok && ossIdx < m_hardMax;
#endif
	order = ossIdx + 1;
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_Uniform::CalculateVariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& variance) const
{
	bool ok = true;

	for (size_t idx = 0; idx < m_numDistinctOrders; idx++)
	{
		//rSums.varEstimate[idx] = 0 < rSums.numSamples[idx] ? (rSums.rad2SumVar[idx] - (pow(rSums.radSum[idx].I(), 2.0) / ((double)rSums.numSamples[idx]))) * pow((double)rSums.numSamples[idx], -1.0) : 0.0;
		variance[idx] = 0 < rSums.numSamples[idx] ? (rSums.rad2SumVar[idx] - (pow(rSums.radSum[idx].I(), 2.0) / ((double)rSums.numSamples[idx]))) * pow((double)rSums.numSamples[idx], -2.0) : 0.0;
	}

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_Uniform::CalculateCovariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& covariance) const
{
	bool ok = true;

	size_t cidx = 0;
	for (size_t lidx = 0; lidx < m_numDistinctOrders; lidx++)
	{
		for (size_t uidx = lidx + 1; uidx < m_numDistinctOrders; uidx++)
		{
			//rSums.covEstimate[cidx++] = 0 < rSums.numSamples[uidx] ? (rSums.rad2SumCov[cidx] - ((rSums.radSum[uidx].I() * rSums.radSum[lidx].I()) / ((double)rSums.numSamples[lidx]))) * pow((double)rSums.numSamples[uidx], -1.0) : 0.0;
			covariance[cidx] = 0 < rSums.numSamples[uidx] ? (rSums.rad2SumCov[cidx] - ((rSums.radSum[uidx].I() * rSums.radSum[lidx].I()) / ((double)rSums.numSamples[lidx]))) * pow((double)rSums.numSamples[lidx], -1.0) * pow((double)rSums.numSamples[uidx], -1.0) : 0.0;
			++cidx;
		}
	}

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_Uniform::CalculateVarianceContribution(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& contribution) const
{
	bool ok = true;
	if (m_numDistinctOrders > 1)
	{
		std::vector<double> var(m_numDistinctOrders - 1);
		double runningVar;
		double runningCov;
		size_t covidx = m_numDistinctOrders - 1;
		for (size_t idx = 1; idx < m_numDistinctOrders; ++idx) {
			runningVar = rSums.rad2SumVar[idx]; // oinfo[idx].variance; this is wrong, should take from orderInfo.UpdateVar()
			runningVar = 0 < rSums.numSamples[idx] ? (rSums.rad2SumVar[idx] - (pow(rSums.radSum[idx].I(), 2.0) / ((double)rSums.numSamples[idx]))) * pow(((double)rSums.numSamples[idx]), -2.0) : 0.0;
			runningCov = 0.0;
			for (size_t cidx = idx + 1; cidx < m_numDistinctOrders; ++cidx, ++covidx)
			{
				if (rSums.numSamples[cidx] > 0) runningCov += 2.0 * (rSums.rad2SumCov[covidx] - (rSums.radSum[idx].I()*(1.0 / rSums.numSamples[idx]) * rSums.radSum[cidx].I())) * pow((double)rSums.numSamples[cidx], -1.0) *  pow((double)rSums.numSamples[idx], -1.0);
			}
			var[idx - 1] = sqrt(std::max(runningVar + runningCov, 0.0)); // Assign to variance vector
		}
		//var.back() = sqrt(m_oinfo[MC_NUMDISTINCTORDERS-1].variance + (runningCov*m_oinfo[MC_NUMDISTINCTORDERS-1].variance/m_oinfo[MC_NUMDISTINCTORDERS-2].variance) );	// Estimate cov amongst higher orders as prop to var
		for (size_t idx = 0;idx < m_numDistinctOrders - 2; ++idx) {
			contribution[idx] = var[idx] - var[idx + 1];
		}
		contribution.back() = var[m_numDistinctOrders - 2];
		double sum = 0.0;
		for (size_t idx = 0; idx < contribution.size(); ++idx) {
			sum += pow(contribution[idx], 1.0);//* ((double)(idx+2));
			contribution[idx] = sum;
		}
		ok = ok && sum > 0;
	}
	return ok;
}

std::string SKTRAN_OptimalScatterSequenceManager_Uniform::PrintSequence(size_t varIdx) const
{
	return SKTRAN_OptimalScatterSequenceManager_Base::PrintSequence(0, varIdx + 1, false);
}

bool SKTRAN_OptimalScatterSequenceManager_Uniform::ConfigureRunningSums(SKTRAN_MCThreadRadianceLogger::runningSums & rSums) const
{
	size_t nvar = m_numDistinctOrders;
	size_t ncov = m_numDistinctOrders * (m_numDistinctOrders - 1) / 2;
	size_t nseq = m_hardMax > 1 ? m_numDistinctOrders - 1 : m_numDistinctOrders;
	rSums.radBuffer.resize(nvar);
	rSums.radSum.resize(nvar);
	rSums.rad2SumVar.resize(nvar);
	rSums.rad2SumCov.resize(ncov);
	rSums.numSamples.resize(nvar);
	rSums.varEstimate.resize(nvar);
	rSums.covEstimate.resize(ncov);
	rSums.varDerivative.resize(nvar);
	rSums.covDerivative.resize(ncov);
	rSums.varContribution.resize(nseq);
	rSums.minSamplesComplete = false;
	return true;
}

bool SKTRAN_OptimalScatterSequenceManager_Uniform::CalculateMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, SKTRAN_Stokes_NC& measurement) const
{
	bool ok = true;
	measurement.SetTo(0.0);
	for (size_t idx = 0; idx < m_numDistinctOrders; idx++)
	{
		if (0 < rSums.numSamples[idx])
		{
			measurement += rSums.radSum[idx] * (1.0 / rSums.numSamples[idx]);
		}
	}
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_Uniform::CalculateTotalVariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, double& variance) const
{
	bool ok = true;
	variance = 0.0;
	for (size_t idx = 0; idx < m_numDistinctOrders; idx++)
	{
		variance += rSums.varEstimate[idx];
	}
	for (size_t cidx = 0; cidx < m_numDistinctOrders * (m_numDistinctOrders - 1) / 2; cidx++)
	{
		variance += 2.0 * rSums.covEstimate[cidx];
	}
	return ok;
}


bool SKTRAN_OptimalScatterSequenceManager_Uniform::CalculateTargetMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, double & measurement) const
{
	bool ok = true;
	SKTRAN_Stokes_NC vec;
	ok = ok && CalculateMeasurement(rSums, vec);
	measurement = vec.I();
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_Uniform::CalculateTargetVariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, double & variance, size_t& numSamples) const
{
	bool ok = true;
	ok = ok && CalculateTotalVariance(rSums, variance);
	numSamples = rSums.numSamples[0]; // number of first order samples represents total number of samples
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_Uniform::ElasticScatter(size_t ossIdx, size_t order, bool& elastic) const
{
	elastic = true;
	return true;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::Combination(size_t n, size_t r, size_t& c) const
{
	bool ok = true;

	r = n - r < r ? n - r : r; // use the equivalent version with smaller r
	c = 1;

	size_t m = std::numeric_limits<size_t>::max() / n; // keeping c < m guarantees the unsigned int won't loop around

	size_t f = n - r + 1; // starting factor in the numerator
	size_t d = 1; // starting factor in the denominator
	while (ok && f <= n)
	{
		c *= f++;
		c /= d++; // guaranteed to be divisible since we have already multiplied d consecutive numbers in the numerator
		ok = c < m;
	}

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::SequenceIterator(size_t T, size_t R, size_it& totalOrder, size_it& ramanOrder, uint64_it& sequence) const
{
	return SequenceIterator(T, R, 0, 0, 0, totalOrder, ramanOrder, sequence);
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::SequenceIterator(size_t T, size_t R, size_t t, size_t r, uint64_t n, size_it& totalOrder, size_it& ramanOrder, uint64_it& sequence) const
{
	// this function builds sequences of elastic scatters (0s) and Raman scatters (1s) recursively
	// T: number of total scatters remaining
	// R: maximum number of Raman scatters remaining
	// t: total number of scatters decided so far
	// r: number of Raman scatters decided so far
	// n: binary number representing the sequence so far

	bool ok = true;

	if (R == 0)
	{
		*(totalOrder++) = t + T;
		*(ramanOrder++) = r;
		*(sequence++) = n;
	}
	else if (T == 1)
	{
		*(totalOrder++) = t + 1;
		*(ramanOrder++) = r;
		*(sequence++) = n;
		*(totalOrder++) = t + 1;
		*(ramanOrder++) = r + 1;
		*(sequence++) = n + 1;
	}
	else
	{
		ok = ok && SequenceIterator(T - 1, R, t + 1, r, n, totalOrder, ramanOrder, sequence);
		ok = ok && T < 64;
		ok = ok && SequenceIterator(T - 1, R - 1, t + 1, r + 1, n + m_bits[T - 1], totalOrder, ramanOrder, sequence);
	}

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::CalculateMaxOrders(const std::vector<size_t>& maxTotalOrders, std::vector<size_t>& maxRamanOrders) const
{
	bool ok = true;
	const std::vector<size_t>& S = maxTotalOrders;
	std::vector<size_t>& R = maxRamanOrders;

	// make sure S is decreasing
	for (size_t i = 0; i < S.size() - 1; i++)
	{
		ok = ok && S[i] >= S[i + 1];
	}
	if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::CalculateMaxOrders, Max total orders not decreasing.");

	// S0 is the maximum number of scatters
	size_t S0 = S.front();
	ok = ok && S0 < 65; // can't store the sequence in a 64-bit integer if it's too long

	// R = [R1, R2, ..., RS0], where Rs is the maximum number of Raman scatters for s total scatters
	R.resize(S0);
	size_t s = 0, r = S.size() - 1;
	for (size_t i = 0; ok && i < S0; i++)
	{
		s = i + 1;
		R[i] = s < r ? s : r;
		while (ok && r < S0 && s >= S[r]) r--;
	}

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::CalculateNumberOfTerms(size_t & numSequences, size_t & numCovarianceTerms) const
{
	bool ok = true;

	size_t c, s, r;
	const std::vector<size_t>& R = m_maxRamanOrders;

	numSequences = 0; // total number of sequences (also total number of variance terms)
	numCovarianceTerms = 0; // total number of covariance terms
	for (size_t i = 0; ok && i < R.size(); i++)
	{
		s = i + 1;
		for (r = 0; ok && r <= R[i]; r++)
		{
			ok = ok && Combination(s, r, c);
			numSequences += c;
			numCovarianceTerms += c * (s - 1);
		}
	}
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::DefineSequences(std::vector<uint64_t>& sequences, std::vector<size_t>& seqToTotalOrder, std::vector<size_t>& seqToRamanOrder) const
{
	bool ok = true;
	size_t s;
	const std::vector<size_t>& R = m_maxRamanOrders;
	
	// fill in the sequences (m_seq) and their total number of scatters (m_seqToMaxOrder)
	sequences.resize(m_numSequences);
	seqToTotalOrder.resize(m_numSequences);
	seqToRamanOrder.resize(m_numSequences);

	auto seqIt = sequences.begin();
	auto totalIt = seqToTotalOrder.begin();
	auto ramanIt = seqToRamanOrder.begin();

	for (size_t i = 0; ok && i < R.size(); i++)
	{
		s = i + 1;
		ok = ok && SequenceIterator(s, R[i], totalIt, ramanIt, seqIt);
	}
	ok = ok && seqIt == sequences.end() && totalIt == seqToTotalOrder.end() && ramanIt == seqToRamanOrder.end();

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::DefineCovarianceIndicesLower(std::vector<size_t>& covToLowerSeqIdx, std::vector<size_t>& covToUpperSeqIdx, std::vector<size_t>& seqToCovLowerIdx, std::vector<size_t>& seqToCovLowerSize) const
{
	bool ok = true;

	size_t lOrder, uOrder;
	uint64_t lSeq, uSeq;
	size_t k = 0;

	covToLowerSeqIdx.resize(m_numCovarianceTerms);
	covToUpperSeqIdx.resize(m_numCovarianceTerms);
	seqToCovLowerIdx.resize(m_numSequences);
	seqToCovLowerSize.resize(m_numSequences);

	std::vector<size_t>::iterator il = covToLowerSeqIdx.begin();
	std::vector<size_t>::iterator iu = covToUpperSeqIdx.begin();
	for (size_t i = 0; ok && i < m_numSequences; i++) // loop through all sequences
	{
		lOrder = m_seqToTotalOrder[i];
		lSeq = m_seq[i];
		seqToCovLowerIdx[i] = k;
		for (size_t j = i + 1; ok && j < m_numSequences; j++) // loop through all following sequences
		{
			uOrder = m_seqToTotalOrder[j];
			uSeq = m_seq[j];
			// skip if the lower sequence is not a subset of the upper sequence
			if (uSeq - lSeq * m_bits[uOrder - lOrder] < m_bits[uOrder - lOrder])
			{
				// multiplying the lower number by powers of 2 shifts it to the left, adding 0s on the right
				// example: 001 (lSeq) is a subsequence of 001110 (uSeq), so 001110 (uSeq) - 001000 (lSeq * 2^3) = 000110 < 2^3
				// example: 001 (lSeq) is not a subsequence of 10110 (uSeq), so 10110 (uSeq) - 01000 (lSeq * 2^2) = 00110 >= 2^2
				// example: 010 (lSeq) is not a subsequence of 0011101 (uSeq), so 0011101 (uSeq) - 0100000 (lSeq * 2^4) < 0 (the subtraction wraps around because these are unsigned)
				*(il++) = i;
				*(iu++) = j;
				k++;
			}
		}
		seqToCovLowerSize[i] = k - seqToCovLowerIdx[i];
	}
	ok = ok && il == covToLowerSeqIdx.end() && iu == covToUpperSeqIdx.end();

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::DefineCovarianceIndicesUpper(std::vector<size_t>& covGroupByUpperIdx, std::vector<size_t>& seqToCovUpperIdx, std::vector<size_t>& seqToCovUpperSize) const
{
	bool ok = true;

	// reorder the covariance terms to be grouped by the upper term
	if (m_numCovarianceTerms > 0)
	{
		covGroupByUpperIdx.resize(m_numCovarianceTerms);
		for (size_t i = 0; i < m_numCovarianceTerms; i++) covGroupByUpperIdx[i] = i; // fill with ordered indices to start
		std::vector<size_t> covUpperScore(m_numCovarianceTerms); // unique score for each covariance pair
		ok = ok && m_numSequences < std::numeric_limits<size_t>::max() / m_numCovarianceTerms; // make sure the score won't overflow
		for (size_t i = 0; i < m_numCovarianceTerms; i++) covUpperScore[i] = m_covToUpperSeqIdx[i] * m_numSequences + m_covToLowerSeqIdx[i]; // sort by upper idx first, then by lower idx
		std::sort(covGroupByUpperIdx.begin(), covGroupByUpperIdx.end(), [&](size_t i, size_t j) { return covUpperScore[i] < covUpperScore[j]; }); // sort the indices according to the score

		// record the size and position of each group
		seqToCovUpperIdx.resize(m_numSequences);
		seqToCovUpperSize.resize(m_numSequences);
		for (size_t i = 0; i < m_numSequences; i++) seqToCovUpperSize[i] = m_seqToTotalOrder[i] - 1;
		if (m_numSequences > 0) seqToCovUpperIdx[0] = 0;
		for (size_t i = 1; i < m_numSequences; i++) seqToCovUpperIdx[i] = seqToCovUpperIdx[i - 1] + seqToCovUpperSize[i - 1];
	}

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::DefineSequenceIncrements(std::vector<size_t>& seqToSeqPlusCabannes, std::vector<size_t>& seqToSeqPlusRaman) const
{
	bool ok = true;

	size_t lowIdx, highIdx, order, orderIdx, covIdx, covIdxGroupUpper;

	seqToSeqPlusCabannes.resize(m_numSequences);	// maps seqIdx to the index of the sequence with one additional Cabannes scatter
	seqToSeqPlusRaman.resize(m_numSequences);		// maps seqIdx to the index of the sequence with one additional Raman scatter

	// fill with invalid value to start (for cases where the sequence with one additional Cabannes/Raman scatter does not exist)
	std::fill(seqToSeqPlusCabannes.begin(), seqToSeqPlusCabannes.end(), m_numSequences);
	std::fill(seqToSeqPlusRaman.begin(), seqToSeqPlusRaman.end(), m_numSequences);

	for (highIdx = 0; highIdx < m_numSequences; highIdx++)
	{
		// the desired list is already mostly stored in the covariance indices, we're just extracting it for easier access
		order = m_seqToTotalOrder[highIdx];
		if (order > 1)
		{
			orderIdx = order - 2;
			covIdxGroupUpper = m_seqToCovUpperIdx[highIdx] + orderIdx;
			covIdx = m_covGroupByUpperIdx[covIdxGroupUpper];
			lowIdx = m_covToLowerSeqIdx[covIdx];
			if (m_seq[highIdx] & m_bits[0]) seqToSeqPlusRaman[lowIdx] = highIdx;	// if the last scatter is Raman
			else seqToSeqPlusCabannes[lowIdx] = highIdx;							// if the last scatter is Cababnes	
		}
	}

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::CalculateMinimumFractions(size_t& numMaxSeq, std::vector<size_t>& maxSeq, std::vector<double>& maxSeqToMinFraction, std::vector<size_t>& maxSeqToMinSamples) const
{
	bool ok = true;

	if (m_maxTotalOrders.size() > 0 && m_minFractionHigherOrder.size() > 0)
	{
		ok = ok && m_maxTotalOrders.size() == m_minFractionHigherOrder.size();
		numMaxSeq = 0;
		std::vector<size_t> numMaxSeqPerOrder(m_maxTotalOrders.size());
		for (size_t i = 0; ok && i < m_maxTotalOrders.size(); i++)
		{
			ok = ok && Combination(m_maxTotalOrders[i], i, numMaxSeqPerOrder[i]);
			numMaxSeq += numMaxSeqPerOrder[i];
		}

		maxSeq.resize(numMaxSeq);
		maxSeqToMinFraction.resize(numMaxSeq);
		maxSeqToMinSamples.resize(numMaxSeq);

		auto maxIt = maxSeq.begin();
		auto minFracIt = maxSeqToMinFraction.begin();
		auto minSampIt = maxSeqToMinSamples.begin();

		size_t r;
		for (size_t i = 0; i < m_seq.size(); i++)
		{ 
			r = m_seqToRamanOrder[i];
			ok = ok && r < m_maxTotalOrders.size();
			if (m_maxTotalOrders[r] == m_seqToTotalOrder[i])
			{
				*(maxIt++) = i;
				*(minFracIt++) = m_minFractionHigherOrder[r] / numMaxSeqPerOrder[r];
				*(minSampIt++) = m_minSamplesHigherOrder;
			}
		}
		ok = ok && maxIt == maxSeq.end() && minFracIt == maxSeqToMinFraction.end() && minSampIt == maxSeqToMinSamples.end();
	}

	return ok;
}

SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic()
{
	m_minSamplesHigherOrder = 20;
	m_numSequences = 0;
	m_numCovarianceTerms = 0;
	m_bits[0] = 1;
	for (size_t i = 0; i < 62; i++) m_bits[i + 1] = 2 * m_bits[i];
	m_hardMax = 0; // this must be set
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::SetMaxOrder(size_t hardmax)
{
	bool ok = true;

	if (hardmax != m_hardMax)
	{
		ok = ok && SKTRAN_OptimalScatterSequenceManager_Base::SetMaxOrder(hardmax);
		m_maxTotalOrders[0] = m_numDistinctOrders;

		if (hardmax < m_maxTotalOrders.size() - 1)
		{
			nxLog::Record(NXLOG_WARNING, "SKTRAN_OptimalScatterSequenceManager_Inelastic::SetMaxOrder, The number of Raman orders set by SKTRAN_Specifications_MC::SetMaxRamanOrders exceeded the hard maximum specified in the call to SKTRAN_Engine_MC_V21::CalculateRadiance - it (and other vectors arranged by Raman order) have been truncated to the hard maximum + 1.");
			m_maxTotalOrders.resize(hardmax + 1);
			m_minFractionHigherOrder.resize(hardmax + 1);
		}

		for (size_t i = 0; i < m_maxTotalOrders.size(); i++)
		{
			if (m_maxTotalOrders[i] > hardmax)
			{
				m_maxTotalOrders[i] = hardmax;
				nxLog::Record(NXLOG_WARNING, "SKTRAN_OptimalScatterSequenceManager_Inelastic::SetMaxOrder, A maximum order specified by SKTRAN_Specifications_MC::SetMaxRamanOrders exceeded the hard maximum specified in the call to SKTRAN_Engine_MC_V21::CalculateRadiance - it has been truncated to the hard maximum.");
			}
		}

		ok = ok && CalculateMaxOrders(m_maxTotalOrders, m_maxRamanOrders);
		ok = ok && CalculateNumberOfTerms(m_numSequences, m_numCovarianceTerms);
		ok = ok && DefineSequences(m_seq, m_seqToTotalOrder, m_seqToRamanOrder);
		ok = ok && DefineCovarianceIndicesLower(m_covToLowerSeqIdx, m_covToUpperSeqIdx, m_seqToCovLowerIdx, m_seqToCovLowerSize);
		ok = ok && DefineCovarianceIndicesUpper(m_covGroupByUpperIdx, m_seqToCovUpperIdx, m_seqToCovUpperSize);
		ok = ok && DefineSequenceIncrements(m_seqToSeqPlusCabannes, m_seqToSeqPlusRaman);
		//ok = ok && CalculateMinimumFractions(m_numMaxSequences, m_maxSeqIdx, m_maxSeqToMinFraction, m_maxSeqToMinSamples);
		ok = ok && SetMinFractionHigherOrder(m_minFractionHigherOrder); // this function resets m_minFractionHigherOrder if it was truncated AND calls CalculateMinimumFractions
		m_covToLowerVarIdx = m_covToLowerSeqIdx; // without elastic Raman samples, 
		m_covToUpperVarIdx = m_covToUpperSeqIdx; // sequences and variance terms are one to one
	}

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::SetMinFractionHigherOrder(const std::vector<double>& minFrac)
{
	bool ok = true;
	m_minFractionHigherOrder = minFrac;
	m_minFractionHigherOrderSum = 0.0;
	for (auto&& frac : minFrac) m_minFractionHigherOrderSum += frac;
	ok = ok && CalculateMinimumFractions(m_numMaxSequences, m_maxSeqIdx, m_maxSeqToMinFraction, m_maxSeqToMinSamples);
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::SetMaxRamanOrders(const std::vector<size_t>& maxOrders)
{
	bool ok = true;

	m_maxTotalOrders = maxOrders;												// maxOrders contains the max total order for each Raman order, starting from 1
	m_maxTotalOrders.insert(m_maxTotalOrders.begin(), m_numDistinctOrders);		// add the max total order for Raman order 0 (capped by the number of distinct orders to consider)

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::OptimalScatterSequenceIndex(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double randNum, size_t& ossIdx, bool& minSamplesComplete) const
{
	bool ok = true;

	size_t maxsize_t = -1;
	size_t seqIdx, n, i;
	size_t minSeqIdx = maxsize_t, minSamples = maxsize_t;

	if (!minSamplesComplete)
	{
		for (i = 0; i < m_numMaxSequences; i++)
		{
			seqIdx = m_maxSeqIdx[i];
			n = rSums.numSamples[seqIdx];
			if (n < m_maxSeqToMinSamples[i])
			{
				ossIdx = seqIdx;
				break;
			}
		}
		if (!(++n < m_maxSeqToMinSamples[i] || ++i < m_numMaxSequences)) minSamplesComplete = true;
		return ok;
	}

	if (randNum < m_minFractionHigherOrderSum)
	{
		bool done = false;
		for (size_t i = 0; !done && i < m_numMaxSequences; i++)
		{
			if (randNum < m_maxSeqToMinFraction[i])
			{
				ossIdx = m_maxSeqIdx[i];
				done = true;
			}
			else
			{
				randNum -= m_maxSeqToMinFraction[i];
			}
		}
		ok = ok && done;
	}
	else
	{	
		const std::vector<double>& x = rSums.varContribution;
		size_t offset = m_hardMax > 1 ? 2 : 0; // don't choose to only sample a first order scatter (unless we are only calculating single scatter)
		randNum = (randNum - m_minFractionHigherOrderSum) / (1.0 - m_minFractionHigherOrderSum); // renormalize randNum
		double target = randNum * x.back();
		auto it = std::upper_bound(x.begin(), x.end(), target);
		ossIdx = std::distance(x.begin(), it) + offset;
	}
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::SortSamples(const size_t & ossIdx, SKTRAN_MCThreadRadianceLogger::runningSums & rSums) const
{
	bool ok = true;

	size_t i, j, idx, cidx, uidx, maxOrder, order;

	idx = ossIdx;
	maxOrder = m_seqToTotalOrder[idx];
	 
	for (i = 0; i < maxOrder; i++)
	{
		order = m_seqToTotalOrder[idx];

		rSums.radSum[idx] += rSums.radBuffer[order - 1];
		rSums.rad2SumVar[idx] += std::pow(rSums.radBuffer[order - 1].I(), 2.0);
		rSums.numSamples[idx]++;

		if (order > 1)
		{
			uidx = m_seqToCovUpperIdx[idx];				// cov index (in the group-by-upper-term order) of the first covariance term with upper term idx
			for (j = 0; j < order - 1; j++)
			{
				cidx = m_covGroupByUpperIdx[uidx++];	// cov index (in the group-by-lower-term order)
				ok = ok && m_covToUpperSeqIdx[cidx] == idx;	// sanity check
				rSums.rad2SumCov[cidx] += rSums.radBuffer[order - 1].I() * rSums.radBuffer[j].I();
			}
			idx = m_covToLowerSeqIdx[cidx]; // start the next loop with the upper term reduced by one order
		}
	}

	for (i = 0; i < rSums.radBuffer.size(); i++) rSums.radBuffer[i].SetTo(0.0);

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::SubmitSample(const size_t ossIdx, const size_t order, const SKTRAN_MCPhoton_Base * mcphoton, std::vector<SKTRAN_MCThreadRadianceLogger::runningSums>& rSums) const
{
	bool ok = true;
	size_t idx = order < m_numDistinctOrders ? order - 1 : m_numDistinctOrders - 1; // order > m_numDistinctOrders should only occur for pure Cabannes sequences

	auto pr = mcphoton->photonRadiances().cbegin();
	for (auto rs = rSums.begin(); rs != rSums.end(); rs++, pr++)
	{
		rs->radBuffer[idx] += pr->GetRecentContribVec();
	}
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::Order(size_t ossIdx, size_t& order) const
{
	bool ok = true;
	ok = ok && ossIdx < m_numSequences;
	if (ok) order = m_seqToTotalOrder[ossIdx];
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::CalculateVariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& variance) const
{
	bool ok = true;

	for (size_t idx = 0; idx < m_numSequences; idx++)
	{
		//rSums.varEstimate[idx] = 0 < rSums.numSamples[idx] ? (rSums.rad2SumVar[idx] - (pow(rSums.radSum[idx].I(), 2.0) / ((double)rSums.numSamples[idx]))) * pow((double)rSums.numSamples[idx], -1.0) : 0.0;
		variance[idx] = 0 < rSums.numSamples[idx] ? (rSums.rad2SumVar[idx] - (pow(rSums.radSum[idx].I(), 2.0) / ((double)rSums.numSamples[idx]))) * pow((double)rSums.numSamples[idx], -2.0) : 0.0;
	}

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::CalculateCovariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& covariance) const
{
	bool ok = true;
	size_t lidx, uidx;

	for (size_t cidx = 0; cidx < m_numCovarianceTerms; cidx++)
	{
		lidx = m_covToLowerSeqIdx[cidx];
		uidx = m_covToUpperSeqIdx[cidx];
		//rSums.covEstimate[cidx] = 0 < rSums.numSamples[uidx] ? (rSums.rad2SumCov[cidx] - ((rSums.radSum[uidx].I() * rSums.radSum[lidx].I()) / ((double)rSums.numSamples[lidx]))) * pow((double)rSums.numSamples[uidx], -1.0) : 0.0;
		covariance[cidx] = 0 < rSums.numSamples[uidx] ? (rSums.rad2SumCov[cidx] - ((rSums.radSum[uidx].I() * rSums.radSum[lidx].I()) / ((double)rSums.numSamples[lidx]))) * pow((double)rSums.numSamples[lidx], -1.0) * pow((double)rSums.numSamples[uidx], -1.0) : 0.0;
	}

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::CalculateVarianceContribution(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& contribution) const
{
	size_t offset = m_hardMax > 1 ? 2 : 0; // don't choose to only sample a first order scatter (unless we are only calculating single scatter)
	std::vector<double> w(m_numSequences - offset); // sum of variance and covariance terms
	std::vector<double>& x = contribution; // 

	for (size_t lowIdx = offset; lowIdx < m_numSequences; lowIdx++)
	{
		size_t wIdx = lowIdx - offset;
		w[wIdx] = rSums.varEstimate[lowIdx] * rSums.numSamples[lowIdx];
		for (size_t covIdx = m_seqToCovLowerIdx[lowIdx]; covIdx < m_seqToCovLowerIdx[lowIdx] + m_seqToCovLowerSize[lowIdx]; covIdx++)
		{
			w[wIdx] += 2.0 * rSums.covEstimate[covIdx] * rSums.numSamples[lowIdx];
		}
		w[wIdx] = sqrt(std::max(w[wIdx], 0.0));
	}

	double sum = 0.0;
	for (size_t lowIdx = offset; lowIdx < m_numSequences; lowIdx++)
	{
		size_t xIdx = lowIdx - offset;
		size_t cIdx = m_seqToSeqPlusCabannes[lowIdx];
		size_t rIdx = m_seqToSeqPlusRaman[lowIdx];

		sum += w[xIdx];
		if (cIdx < m_numSequences) sum -= w[cIdx - offset];
		if (rIdx < m_numSequences) sum -= w[rIdx - offset];
		x[xIdx] = sum;
	}
	return sum > 0;
}

std::string SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::PrintSequence(size_t varidx) const
{
	return SKTRAN_OptimalScatterSequenceManager_Base::PrintSequence(m_seq[varidx], m_seqToTotalOrder[varidx], false);
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::ConfigureRunningSums(SKTRAN_MCThreadRadianceLogger::runningSums & rSums) const
{
	size_t nvar = m_numSequences;
	size_t ncov = m_numCovarianceTerms;
	size_t nseq = m_hardMax > 1 ? m_numSequences - 2 : m_numSequences;
	rSums.radBuffer.resize(nvar);
	rSums.radSum.resize(nvar);
	rSums.rad2SumVar.resize(nvar);
	rSums.rad2SumCov.resize(ncov);
	rSums.numSamples.resize(nvar);
	rSums.varEstimate.resize(nvar);
	rSums.covEstimate.resize(ncov);
	rSums.varDerivative.resize(nvar);
	rSums.covDerivative.resize(ncov);
	rSums.varContribution.resize(nseq);
	return true;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::CalculateMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, SKTRAN_Stokes_NC& measurement) const
{
	bool ok = true;
	skRTStokesVector::SetToZero(measurement);

	for (size_t idx = 0; idx < m_numSequences; idx++)
	{
		if (0 < rSums.numSamples[idx])
		{
			measurement += rSums.radSum[idx] * (1.0 / rSums.numSamples[idx]);
		}
	}
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::CalculateTotalVariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, double& variance) const
{
	bool ok = true;
	variance = 0.0;
	for (size_t idx = 0; idx < m_numSequences; idx++)
	{
		variance += rSums.varEstimate[idx];
	}
	for (size_t cidx = 0; cidx < m_numCovarianceTerms; cidx++)
	{
		variance += 2.0 * rSums.covEstimate[cidx];
	}
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::CalculateTargetMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, double & measurement) const
{
	bool ok = true;
	SKTRAN_Stokes_NC vec;
	ok = ok && CalculateMeasurement(rSums, vec);
	measurement = vec.I();
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::CalculateTargetVariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, double & variance, size_t& numSamples) const
{
	bool ok = true;

	ok = ok && CalculateTotalVariance(rSums, variance);
	numSamples = rSums.numSamples[0] + rSums.numSamples[1]; // number of first order samples (cabannes + raman) represents total number of samples

	return ok;
}

//bool SKTRAN_OptimalScatterSequenceManager_Inelastic::VarianceIndex(size_t ossIdx, size_t order, size_t& index) const
//{
//	return false;
//}
//
//bool SKTRAN_OptimalScatterSequenceManager_Inelastic::CovarianceIndex(size_t ossIdx, size_t order1, size_t order2, size_t& index) const
//{
//	return false;
//}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::ElasticScatter(size_t ossIdx, size_t order, bool& elastic) const
{
	bool ok = true;
	size_t maxOrder;
	ok = ok && Order(ossIdx, maxOrder);
	ok = ok && order > 0 && order <= maxOrder && ossIdx < m_numSequences;
	if (ok) elastic = !(m_seq[ossIdx] & m_bits[maxOrder - order]);
	return ok;
}

SKTRAN_OptimalScatterSequenceManager_UniformSecondary::SKTRAN_OptimalScatterSequenceManager_UniformSecondary()
{
}

bool SKTRAN_OptimalScatterSequenceManager_UniformSecondary::SetMaxOrder(size_t hardmax)
{
	bool ok = true;
	ok = ok && SKTRAN_OptimalScatterSequenceManager_Base::SetMaxOrder(hardmax);

	size_t numCov = m_numDistinctOrders * (2 * m_numDistinctOrders - 1);
	size_t numVar = 2 * m_numDistinctOrders;
	size_t cidx = 0;
	m_covToLowerVarIdx.resize(numCov);
	m_covToUpperVarIdx.resize(numCov);
	for (size_t lidx = 0; lidx < numVar; lidx++)
	{
		for (size_t uidx = lidx + 1; uidx < numVar; uidx++)
		{
			m_covToLowerVarIdx[cidx] = lidx;
			m_covToUpperVarIdx[cidx++] = uidx;
		}
	}

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_UniformSecondary::SortSamples(const size_t & ossIdx, SKTRAN_MCThreadRadianceLogger::runningSums & rSums) const
{
	bool ok = true;

	size_t maxOrder = ossIdx < m_numDistinctOrders ? ossIdx + 1 : m_numDistinctOrders;
	size_t j = m_numDistinctOrders;
	for (size_t i = 0; i < maxOrder; i++, j++)
	{
		rSums.numSamples[i]++;
		rSums.numSamples[j]++;
	}

	size_t numVar = 2 * m_numDistinctOrders;
	size_t cidx = 0;
	for (size_t i = 0; i < numVar; i++)
	{
		rSums.radSum[i] += rSums.radBuffer[i];
		rSums.rad2SumVar[i] += pow(rSums.radBuffer[i].I(), 2.0);
		for (size_t j = i + 1; j < numVar; j++)
		{
			rSums.rad2SumCov[cidx++] += rSums.radBuffer[i].I() * rSums.radBuffer[j].I();
		}
	}

	for (size_t i = 0; i < rSums.radBuffer.size(); i++) rSums.radBuffer[i].SetTo(0.0);

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_UniformSecondary::SubmitSample(const size_t ossIdx, const size_t order, const SKTRAN_MCPhoton_Base * mcphoton, std::vector<SKTRAN_MCThreadRadianceLogger::runningSums>& rSums) const
{
	bool ok = true;
	size_t idx = order < m_numDistinctOrders ? order - 1 : m_numDistinctOrders - 1;
	size_t eidx = idx + m_numDistinctOrders;

	auto pr = mcphoton->photonRadiances(false).cbegin();
	auto epr = mcphoton->photonRadiances(true).cbegin();
	for (auto rs = rSums.begin(); rs != rSums.end(); rs++, pr++, epr++)
	{
		rs->radBuffer[idx] += pr->GetRecentContribVec();
		rs->radBuffer[eidx] += epr->GetRecentContribVec();
	}
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_UniformSecondary::ConfigureRunningSums(SKTRAN_MCThreadRadianceLogger::runningSums & rSums) const
{

	size_t nvar = 2 * m_numDistinctOrders;
	size_t ncov = m_numDistinctOrders * (2 * m_numDistinctOrders - 1);
	size_t nseq = m_hardMax > 1 ? m_numDistinctOrders - 1 : m_numDistinctOrders;
	rSums.radBuffer.resize(nvar);
	rSums.radSum.resize(nvar);
	rSums.rad2SumVar.resize(nvar);
	rSums.rad2SumCov.resize(ncov);
	rSums.numSamples.resize(nvar);
	rSums.varEstimate.resize(nvar);
	rSums.covEstimate.resize(ncov);
	rSums.varDerivative.resize(nvar);
	rSums.covDerivative.resize(ncov);
	rSums.varContribution.resize(nseq);
	return true;
}

bool SKTRAN_OptimalScatterSequenceManager_UniformSecondary::CalculateTargetMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, double & measurement) const
{
	return CalculateSecondaryMeasurement(rSums, measurement);
}

bool SKTRAN_OptimalScatterSequenceManager_UniformSecondary::CalculateTargetVariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, double & variance, size_t& numSamples) const
{
	bool ok = true;

	ok = ok && CalculateSecondaryVariance(rSums, variance);
	numSamples = rSums.numSamples[0]; // number of first order samples represents total number of samples

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_UniformSecondary::CalculateVariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& variance) const
{
	bool ok = true;

	ok = ok && SKTRAN_OptimalScatterSequenceManager_Uniform::CalculateVariance(rSums, variance);

	for (size_t idx = m_numDistinctOrders; idx < 2 * m_numDistinctOrders; idx++)
	{
		//rSums.varEstimate[idx] = 0 < rSums.numSamples[idx] ? (rSums.rad2SumVar[idx] - (pow(rSums.radSum[idx].I(), 2.0) / ((double)rSums.numSamples[idx]))) * pow((double)rSums.numSamples[idx], -1.0) : 0.0;
		variance[idx] = 0 < rSums.numSamples[idx] ? (rSums.rad2SumVar[idx] - (pow(rSums.radSum[idx].I(), 2.0) / ((double)rSums.numSamples[idx]))) * pow((double)rSums.numSamples[idx], -2.0) : 0.0;
	}

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_UniformSecondary::CalculateCovariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& covariance) const
{
	bool ok = true;

	size_t numVar = 2 * m_numDistinctOrders;
	size_t cidx = 0;
	for (size_t i = 0; i < numVar; i++)
	{
		for (size_t j = i + 1; j < numVar; j++)
		{
			//rSums.covEstimate[cidx] = 0 < rSums.numSamples[j] ? (rSums.rad2SumCov[cidx] - ((rSums.radSum[j].I() * rSums.radSum[i].I()) / ((double)rSums.numSamples[i]))) * pow((double)rSums.numSamples[j], -1.0) : 0.0;
			covariance[cidx] = 0 < rSums.numSamples[j] ? (rSums.rad2SumCov[cidx] - ((rSums.radSum[j].I() * rSums.radSum[i].I()) / ((double)rSums.numSamples[i]))) * pow((double)rSums.numSamples[i], -1.0) * pow((double)rSums.numSamples[j], -1.0) : 0.0;
			cidx++;
		}
	}
	
	return ok;
}

std::string SKTRAN_OptimalScatterSequenceManager_UniformSecondary::PrintSequence(size_t varIdx) const
{
	bool elasticRaman = varIdx >= m_numDistinctOrders;
	size_t order = elasticRaman ? varIdx - m_numDistinctOrders + 1 : varIdx + 1;
	return SKTRAN_OptimalScatterSequenceManager_Base::PrintSequence(0, order, elasticRaman);
}

bool SKTRAN_OptimalScatterSequenceManager_UniformRing::CalculateSecondaryMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, double& measurement) const
{
	bool ok = true;

	double elasticResult = 0.0, inelasticResult = 0.0;

	for (size_t idx = 0; idx < m_numDistinctOrders; idx++)
	{
		if (0 < rSums.numSamples[idx])
		{
			inelasticResult += rSums.radSum[idx].I() * (1.0 / rSums.numSamples[idx]);
		}
	}

	for (size_t idx = m_numDistinctOrders; idx < 2 * m_numDistinctOrders; idx++)
	{
		if (0 < rSums.numSamples[idx])
		{
			elasticResult += rSums.radSum[idx].I() * (1.0 / rSums.numSamples[idx]);
		}
	}

	measurement = log(inelasticResult / elasticResult);

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_UniformRing::CalculateSecondaryVariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, double& variance) const
{
	double elasticResult = 0.0, inelasticResult = 0.0;
	size_t i, j, idx;

	for (i = 0; i < m_numDistinctOrders; i++)
	{
		idx = i;
		if (0 < rSums.numSamples[idx])
		{
			inelasticResult += rSums.radSum[idx].I() * (1.0 / rSums.numSamples[idx]);
		}
	}

	for (i = 0; i < m_numDistinctOrders; i++)
	{
		idx = i + m_numDistinctOrders;
		if (0 < rSums.numSamples[idx])
		{
			elasticResult += rSums.radSum[idx].I() * (1.0 / rSums.numSamples[idx]);
		}
	}

	double inelasticDerivative = 1.0 / inelasticResult;
	double elasticDerivative = -1.0 / elasticResult;

	size_t numVar = 2 * m_numDistinctOrders;
	size_t cidx = 0;
	double di, dj;
	variance = 0.0;
	for (i = 0; i < numVar; i++)
	{
		di = i < m_numDistinctOrders ? inelasticDerivative : elasticDerivative;
		variance += rSums.varEstimate[i] * di * di;
		for (j = i + 1; j < numVar; j++)
		{
			dj = j < m_numDistinctOrders ? inelasticDerivative : elasticDerivative;
			variance += rSums.covEstimate[cidx++] * di * dj;
		}
	}
	
	return true;
}

bool SKTRAN_OptimalScatterSequenceManager_UniformFillingIn::CalculateSecondaryMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, double & measurement) const
{
	bool ok = true;

	double elasticResult = 0.0, inelasticResult = 0.0;

	for (size_t idx = 0; idx < m_numDistinctOrders; idx++)
	{
		if (0 < rSums.numSamples[idx])
		{
			inelasticResult += rSums.radSum[idx].I() * (1.0 / rSums.numSamples[idx]);
		}
	}

	for (size_t idx = m_numDistinctOrders; idx < 2 * m_numDistinctOrders; idx++)
	{
		if (0 < rSums.numSamples[idx])
		{
			elasticResult += rSums.radSum[idx].I() * (1.0 / rSums.numSamples[idx]);
		}
	}

	measurement = inelasticResult / elasticResult - 1;

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_UniformFillingIn::CalculateSecondaryVariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, double & variance) const
{
	double elasticResult = 0.0, inelasticResult = 0.0;
	size_t i, j, idx;

	for (i = 0; i < m_numDistinctOrders; i++)
	{
		idx = i;
		if (0 < rSums.numSamples[idx])
		{
			inelasticResult += rSums.radSum[idx].I() * (1.0 / rSums.numSamples[idx]);
		}
	}

	for (i = 0; i < m_numDistinctOrders; i++)
	{
		idx = i + m_numDistinctOrders;
		if (0 < rSums.numSamples[idx])
		{
			elasticResult += rSums.radSum[idx].I() * (1.0 / rSums.numSamples[idx]);
		}
	}

	double inelasticDerivative = 1.0 / elasticResult;
	double elasticDerivative = -inelasticDerivative * inelasticDerivative * inelasticResult;

	size_t numVar = 2 * m_numDistinctOrders;
	size_t cidx = 0;
	double di, dj;
	variance = 0.0;
	for (i = 0; i < numVar; i++)
	{
		di = i < m_numDistinctOrders ? inelasticDerivative : elasticDerivative;
		variance += rSums.varEstimate[i] * di * di;
		for (j = i + 1; j < numVar; j++)
		{
			dj = j < m_numDistinctOrders ? inelasticDerivative : elasticDerivative;
			variance += rSums.covEstimate[cidx++] * di * dj;
		}
	}

	return true;
}

bool SKTRAN_OptimalScatterSequenceManager_UniformElastic::CalculateSecondaryMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, double & measurement) const
{
	bool ok = true;
	measurement = 0.0;
	for (size_t idx = m_numDistinctOrders; idx < 2 * m_numDistinctOrders; idx++)
	{
		if (0 < rSums.numSamples[idx])
		{
			measurement += rSums.radSum[idx].I() * (1.0 / rSums.numSamples[idx]);
		}
	}
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_UniformElastic::CalculateSecondaryVariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, double & variance) const
{
	bool ok = true;
	variance = 0.0;
	for (size_t idx = m_numDistinctOrders; idx < 2 * m_numDistinctOrders; idx++)
	{
		variance += rSums.varEstimate[idx];
	}
	size_t numCovTerms = m_numDistinctOrders * (m_numDistinctOrders - 1) / 2;
	for (size_t cidx = numCovTerms; cidx < 2 * numCovTerms; cidx++)
	{
		variance += 2.0 * rSums.covEstimate[cidx];
	}
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary::CalculateNumberOfTerms(size_t & numVarianceTerms, size_t & numCovarianceTerms, size_t& numSubSequences, size_t& numSubVarianceTerms) const
{
	bool ok = true;

	size_t c, s, r;
	const std::vector<size_t>& R = m_maxRamanOrders;

	numVarianceTerms = 0; // total number of sequences (also total number of variance terms)
	numCovarianceTerms = 0; // total number of covariance terms
	numSubSequences = 0;
	numSubVarianceTerms = 0;
	for (size_t i = 0; ok && i < R.size(); i++)
	{
		s = i + 1; // s: total number of scatter
		// no Raman scatters:
		numVarianceTerms += 1; // inelastic and elastic Raman samples are the same, since there are no Raman scatters
		numCovarianceTerms += s - 1;
		numSubSequences += s;
		numSubVarianceTerms += s;
		for (r = 1; ok && r <= R[i]; r++) // r: number of Raman scatters
		{
			ok = ok && Combination(s, r, c);
			numVarianceTerms += 2 * c; // one seqeuence for inelastic Raman, one for elastic 
			numSubSequences += c * s;
			numSubVarianceTerms += s + r; // this accounts for the sequence with all the Ramans stacked at then end (ex: s=5, r=3: CCRRR)
			numCovarianceTerms += 2 * (s + r - 2) + 1; // also stacked at the end: two covariance terms per subsequence where the lower term has no Raman 2(s - r); four per sequence where it does 4(r - 1); one extra for the inelastic/elastic top term
			for (size_t j = 0; j < s - r; j++) // this accounts for the sequences with the first Raman scatter at the (j+1)th position (ex: s=5, r=3, j=1: CRCRR, CRRCR, CRRRC)
			{
				ok = ok && Combination(s - j - 1, r - 1, c); // c: # of combinations with s total and r raman where the (j+1)th scatter is the first Raman scatter
				numSubVarianceTerms += c * (2 * s - j); // 2 subvariance terms (one inelastic, one elastic) for every subsequence containing a Raman scatter (including the main sequence); 1 term for every subsequence without a Raman scatter
				// numCovarianceTerms += c * (4 * s - 3 - j); // 4 covariance terms (upper/lower inelastic/elastic) for every subsequence containing a Raman scatter (not including the main sequence); 2 terms for every subsequence without a Raman scatter; one extra term for the inelastic/elastic top term
				numCovarianceTerms += c * (4 * s - 2 * j - 3); // 4 covariance terms (upper/lower inelastic/elastic) for every subsequence containing a Raman scatter (not including the main sequence); 2 terms for every subsequence without a Raman scatter; one extra term for the inelastic/elastic top term
			}
		}
	}

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary::DefineSequences(std::vector<uint64_t>& sequences, std::vector<size_t>& seqToTotalOrder, std::vector<size_t>& seqToRamanOrder, std::vector<bool>& seqToElasticRaman) const
{
	bool ok = true;
	size_t s;
	const std::vector<size_t>& R = m_maxRamanOrders;

	// fill in the sequences and their total number of scatters
	sequences.resize(m_numVarianceTerms);
	seqToTotalOrder.resize(m_numVarianceTerms);
	seqToRamanOrder.resize(m_numVarianceTerms);
	seqToElasticRaman.resize(m_numVarianceTerms);

	auto seqIt = sequences.begin();
	auto totalIt = seqToTotalOrder.begin();
	auto ramanIt = seqToRamanOrder.begin();
	auto elasticIt = seqToElasticRaman.begin();

	for (size_t i = 0; ok && i < R.size(); i++)
	{
		s = i + 1;
		ok = ok && SequenceIterator(s, R[i], totalIt, ramanIt, elasticIt, seqIt);
	}
 	ok = ok && seqIt == sequences.end() && totalIt == seqToTotalOrder.end() && ramanIt == seqToRamanOrder.end() && elasticIt == seqToElasticRaman.end();

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary::MapSequences(std::vector<size_t>& varToSeqIdx, std::vector<size_t>& seqToVarIdx, std::vector<size_t>& seqToVarSize)
{
	bool ok = true;

	varToSeqIdx.resize(m_numVarianceTerms);
	seqToVarIdx.resize(m_numSequences);
	seqToVarSize.resize(m_numSequences);

	size_t varIdx = 0;
	size_t numVarTerms;
	for (size_t seqIdx = 0; seqIdx < m_numSequences; seqIdx++)
	{
		ok = ok && m_seq[seqIdx] == m_var[varIdx] && m_seqToTotalOrder[seqIdx] == m_varToTotalOrder[varIdx];
		if (ok)
		{
			seqToVarIdx[seqIdx] = varIdx;

			numVarTerms = 0;
			while (ok && varIdx < m_numVarianceTerms && m_seq[seqIdx] == m_var[varIdx] && m_seqToTotalOrder[seqIdx] == m_varToTotalOrder[varIdx])
			{
				varToSeqIdx[varIdx] = seqIdx;
				varIdx++;
				numVarTerms++;
			}
			seqToVarSize[seqIdx] = numVarTerms;
		}
	}
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary::DefineCovarianceIndicesLower(std::vector<size_t>& covToLowerSeqIdx, std::vector<size_t>& covToUpperSeqIdx, std::vector<size_t>& seqToCovLowerIdx, std::vector<size_t>& seqToCovLowerSize, std::vector<size_t>& covToLowerVarIdx, std::vector<size_t>& covToUpperVarIdx, std::vector<size_t>& varToCovLowerIdx, std::vector<size_t>& varToCovLowerSize) const
{
	bool ok = true;

	size_t lOrder, uOrder;
	uint64_t lSeq, uSeq;
	size_t covIdx = 0;

	covToLowerSeqIdx.resize(m_numCovarianceTerms);
	covToUpperSeqIdx.resize(m_numCovarianceTerms);     
	seqToCovLowerIdx.resize(m_numSequences);
	seqToCovLowerSize.resize(m_numSequences);
	covToLowerVarIdx.resize(m_numCovarianceTerms);
	covToUpperVarIdx.resize(m_numCovarianceTerms);
	varToCovLowerIdx.resize(m_numVarianceTerms);
	varToCovLowerSize.resize(m_numVarianceTerms);

	size_t lSeqIdx = -1, uSeqIdx = 0;
	std::vector<size_t>::iterator lsit = covToLowerSeqIdx.begin();
	std::vector<size_t>::iterator usit = covToUpperSeqIdx.begin();
	std::vector<size_t>::iterator lvit = covToLowerVarIdx.begin();
	std::vector<size_t>::iterator uvit = covToUpperVarIdx.begin();
	for (size_t lVarIdx = 0; ok && lVarIdx < m_numVarianceTerms; lVarIdx++) // loop through all variance terms
	{
		// if the sequence index doesn't change, we are considering the elastic sequence immediately after the matching inelastic sequence, and these variables do not need to be set again
		if (lSeqIdx != m_varToSeqIdx[lVarIdx]) 
		{
			lSeqIdx = m_varToSeqIdx[lVarIdx];
			lOrder = m_varToTotalOrder[lVarIdx];
			lSeq = m_var[lVarIdx];
			seqToCovLowerIdx[lSeqIdx] = covIdx;
		}
		varToCovLowerIdx[lVarIdx] = covIdx;

		for (size_t uVarIdx = lVarIdx + 1; ok && uVarIdx < m_numVarianceTerms; uVarIdx++) // loop through all following sequences
		{
			uSeqIdx = m_varToSeqIdx[uVarIdx];
			uOrder = m_varToTotalOrder[uVarIdx];
			uSeq = m_var[uVarIdx];
			// skip if the lower sequence is not a subset of the upper sequence
			if (uSeq - lSeq * m_bits[uOrder - lOrder] < m_bits[uOrder - lOrder])
			{
				*(lsit++) = lSeqIdx;
				*(usit++) = uSeqIdx;
				*(lvit++) = lVarIdx;
				*(uvit++) = uVarIdx;
				covIdx++;
			}
		}
		seqToCovLowerSize[lSeqIdx] = covIdx - seqToCovLowerIdx[lSeqIdx];
		varToCovLowerSize[lVarIdx] = covIdx - varToCovLowerIdx[lVarIdx];
	}
	ok = ok && lsit == covToLowerSeqIdx.end() && usit == covToUpperSeqIdx.end() && lvit == covToLowerVarIdx.end() && uvit == covToUpperVarIdx.end();

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary::DefineCovarianceIndicesUpper(std::vector<size_t>& covGroupByUpperIdx, std::vector<size_t>& varToCovUpperIdx, std::vector<size_t>& varToCovUpperSize, std::vector<size_t>& seqToCovUpperIdx, std::vector<size_t>& seqToCovUpperSize) const
{
	bool ok = true;

	// reorder the covariance terms to be grouped by the upper term
	if (m_numCovarianceTerms > 0)
	{
		covGroupByUpperIdx.resize(m_numCovarianceTerms);
		for (size_t i = 0; i < m_numCovarianceTerms; i++) covGroupByUpperIdx[i] = i; // fill with ordered indices to start
		std::vector<size_t> covUpperScore(m_numCovarianceTerms); // unique score for each covariance pair
		ok = ok && m_numSequences < std::numeric_limits<size_t>::max() / m_numCovarianceTerms; // make sure the score won't overflow
		for (size_t i = 0; i < m_numCovarianceTerms; i++) covUpperScore[i] = m_covToUpperVarIdx[i] * m_numVarianceTerms + m_covToLowerVarIdx[i] ; // sort by upper idx first, then by lower idx
		std::sort(covGroupByUpperIdx.begin(), covGroupByUpperIdx.end(), [&](size_t i, size_t j) { return covUpperScore[i] < covUpperScore[j]; }); // sort the indices according to the score

		// record the size and position of each group
		varToCovUpperIdx.resize(m_numVarianceTerms);
		varToCovUpperSize.resize(m_numVarianceTerms);
		seqToCovUpperIdx.resize(m_numSequences);
		seqToCovUpperSize.resize(m_numSequences);
		size_t covIdx = 0, seqIdx = -1;
		for (size_t varIdx = 0; ok && varIdx < m_numVarianceTerms; varIdx++)
		{
			if (varIdx < m_covToUpperVarIdx[covGroupByUpperIdx[covIdx]])
			{
				varToCovUpperIdx[varIdx] = covIdx;
				varToCovUpperSize[varIdx] = 0;
			}
			else
			{
				while (ok && covIdx < m_numCovarianceTerms && varIdx > m_covToUpperVarIdx[covGroupByUpperIdx[covIdx]]) covIdx++;
				varToCovUpperIdx[varIdx] = covIdx;
				while (ok && covIdx < m_numCovarianceTerms && varIdx == m_covToUpperVarIdx[covGroupByUpperIdx[covIdx]]) covIdx++;
				varToCovUpperSize[varIdx] = covIdx - varToCovUpperIdx[varIdx];
			}

			if (seqIdx == m_varToSeqIdx[varIdx]) // if the sequence for this variance term is the same as the last variance term, leave the index and increase the size
			{
				seqToCovUpperSize[seqIdx] += varToCovUpperSize[varIdx];
			}
			else // if this is a new sequence, set the index and the size
			{
				seqIdx = m_varToSeqIdx[varIdx];
				ok = ok && seqIdx < m_numSequences;
				if (ok) seqToCovUpperIdx[seqIdx] = varToCovUpperIdx[varIdx];
				if (ok) seqToCovUpperSize[seqIdx] = varToCovUpperSize[varIdx];
			}
		}
	}

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary::DefineSubSequenceIndices(std::vector<size_t>& subSeq, std::vector<size_t>& seqToSubSeqIdx) const
{
	bool ok = true;

	subSeq.resize(m_numSubSequences); // for every sequence, list the indicies of the sequences leading up to it
	seqToSubSeqIdx.resize(m_numSequences); // used to access the information in subSeq; if you want to list the indicies of all sequences below and including a given sequence (with order s), start at the index given by this array, then iterate forward s-1 times

	// define the sub-sequence array
	auto subSeqIt = subSeq.begin();
	auto seqToSubSeqIt = seqToSubSeqIdx.begin();
	size_t lOrder, uOrder;
	uint64_t divisor;
	for (size_t uSeqIdx = 0; ok && uSeqIdx < m_numSequences; uSeqIdx++)
	{
		*(seqToSubSeqIt++) = std::distance(subSeq.begin(), subSeqIt);
		uOrder = m_seqToTotalOrder[uSeqIdx];
		lOrder = 1;
		divisor = pow(2, uOrder - lOrder); // used to check if the lower sequence is a subsequnce of the upper sequence
		for (size_t lSeqIdx = 0; ok && lSeqIdx < uSeqIdx; lSeqIdx++)
		{
			divisor = pow(2, uOrder - lOrder); // used to check if the lower sequence is a subsequnce of the upper sequence
			if (m_seqToTotalOrder[lSeqIdx] == lOrder && m_seq[lSeqIdx] == m_seq[uSeqIdx] / divisor)
			{
				*(subSeqIt++) = lSeqIdx;
				lOrder++;
			}
		}
		ok = ok && lOrder == uOrder; // this inidcates a sub-sequences of every preceeding order was found
		if (ok) *(subSeqIt++) = uSeqIdx;
	}
	ok = ok && subSeqIt == subSeq.end();
	ok = ok && seqToSubSeqIt == seqToSubSeqIdx.end();

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary::DefineSubVarianceIndices(std::vector<size_t>& subVar, std::vector<size_t>& seqToSubVarIdx, std::vector<size_t>& seqToSubVarSize, std::vector<size_t>& subSeqToSubVarIdx, std::vector<size_t>& subVarToBufferIdx) const
{
	bool ok = true;

	subVar.resize(m_numSubVarianceTerms); // like m_subSeq, but includes the inelastic and elastic version of any sequences containing at least one Raman scatter; used primarily to sort samples into the correct running sums
										  // for every possible sequence, all samples that will accompany a sample of that sequence are listed here
										  // example: when sampling the sequence CRR, you get samples of C, CR, CR', CRR, and CRR' (*)
										  // subVar contains the variance indicies corresponding to a list like this, defined for every sequence
	seqToSubVarIdx.resize(m_numSequences); // maps seqIdx to the start of its group, like the list above (*)
	seqToSubVarSize.resize(m_numSequences); // maps seqIdx to the size of its group, like the list above (*)
	subSeqToSubVarIdx.resize(m_numSubSequences); // maps a sub-sequence in subSeq to the equivalent term in subVar; if the subsequence has at least one Raman scatter, then the elastic version will be found immediately after the location given in this array
	subVarToBufferIdx.resize(m_numSubVarianceTerms); // used to organize the buffer which holds samples as they are being calculated; simply assigns sequential indices to each list as described above (*)

	// define the sub-variance-term array
	auto subVarIt = subVar.begin();
	auto subSeqToSubVarIt = subSeqToSubVarIdx.begin();
	size_t seqIdx = 0;
	for (size_t subSeqIdx = 0; subSeqIdx < m_numSubSequences; subSeqIdx++)
	{
		*(subSeqToSubVarIt++) = std::distance(subVar.begin(), subVarIt);
		seqIdx = m_subSeq[subSeqIdx];
		if (m_seqToRamanOrder[seqIdx] == 0)
		{
			ok = ok && m_seqToVarSize[seqIdx] == 1; // sanity check; should be redundant with m_seqToRamanOrder[seqIdx] == 0
			*(subVarIt++) = m_seqToVarIdx[seqIdx];
		}
		else
		{
			ok = ok && m_seqToVarSize[seqIdx] == 2; // sanity check; should be redundant with m_seqToRamanOrder[seqIdx] != 0
			*(subVarIt++) = m_seqToVarIdx[seqIdx];
			*(subVarIt++) = m_seqToVarIdx[seqIdx] + 1;
		}
	}
	ok = ok && subVarIt == subVar.end();
	ok = ok && subSeqToSubVarIt == subSeqToSubVarIdx.end();

	for (size_t seqIdx = 0; seqIdx < m_numSequences; seqIdx++)
	{
		size_t subSeqIdx = m_seqToSubSeqIdx[seqIdx];
		size_t subVarIdx = subSeqToSubVarIdx[subSeqIdx];
		seqToSubVarIdx[seqIdx] = subVarIdx;
	}

	for (size_t seqIdx = 0; seqIdx < m_numSequences - 1; seqIdx++)
	{
		seqToSubVarSize[seqIdx] = seqToSubVarIdx[seqIdx + 1] - seqToSubVarIdx[seqIdx];
	}
	seqToSubVarSize[m_numSequences - 1] = m_numSubVarianceTerms - seqToSubVarIdx[m_numSequences - 1];

	auto bufferIt = subVarToBufferIdx.begin();
	for (size_t seqIdx = 0; seqIdx < m_numSequences; seqIdx++)
	{
		for (size_t i = 0; i < seqToSubVarSize[seqIdx]; i++) *(bufferIt++) = i;
	}
	ok = ok && bufferIt == subVarToBufferIdx.end();

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary::SequenceIterator(size_t T, size_t R, size_it & totalOrder, size_it & ramanOrder, bool_it & elasticRaman, uint64_it & sequence) const
{
	return SequenceIterator(T, R, 0, 0, 0, totalOrder, ramanOrder, elasticRaman, sequence);
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary::SequenceIterator(size_t T, size_t R, size_t t, size_t r, uint64_t n, size_it & totalOrder, size_it & ramanOrder, bool_it& elasticRaman, uint64_it & sequence) const
{
	// this function builds sequences of elastic scatters (0s) and Raman scatters (1s) recursively
	// T: number of total scatters remaining
	// R: maximum number of Raman scatters remaining
	// t: total number of scatters decided so far
	// r: number of Raman scatters decided so far
	// n: binary number representing the sequence so far

	bool ok = true;

	if (R == 0)
	{		
		*(totalOrder++) = t + T;
		*(ramanOrder++) = r;
		*(elasticRaman++) = false;
		*(sequence++) = n;
		
		if (r > 0)
		{
			*(totalOrder++) = t + T;
			*(ramanOrder++) = r;
			*(elasticRaman++) = true;
			*(sequence++) = n;
		}
	}
	else if (T == 1)
	{
		*(totalOrder++) = t + 1;
		*(ramanOrder++) = r;
		*(elasticRaman++) = false;
		*(sequence++) = n;

		if (r > 0)
		{
			*(totalOrder++) = t + 1;
			*(ramanOrder++) = r;
			*(elasticRaman++) = true;
			*(sequence++) = n;
		}

		*(totalOrder++) = t + 1;
		*(ramanOrder++) = r + 1;
		*(elasticRaman++) = false;
		*(sequence++) = n + 1;

		*(totalOrder++) = t + 1;
		*(ramanOrder++) = r + 1;
		*(elasticRaman++) = true;
		*(sequence++) = n + 1;
	}
	else
	{
		ok = ok && SequenceIterator(T - 1, R, t + 1, r, n, totalOrder, ramanOrder, elasticRaman, sequence);
		ok = ok && T < 64;
		ok = ok && SequenceIterator(T - 1, R - 1, t + 1, r + 1, n + m_bits[T - 1], totalOrder, ramanOrder, elasticRaman, sequence);
	}

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary::CalculateRadiances(const SKTRAN_MCThreadRadianceLogger::runningSums& rSums, double & inelasticRadiance, double & elasticRadiance) const
{
	bool ok = true;

	elasticRadiance = 0.0, inelasticRadiance = 0.0;

	for (size_t varIdx = 0; varIdx < m_numVarianceTerms; varIdx++)
	{
		if (rSums.numSamples[varIdx] > 0)
		{
			if (m_varToRamanOrder[varIdx] > 0)
			{
				if (m_varToElasticRaman[varIdx])
				{
					elasticRadiance += rSums.radSum[varIdx].I() * (1.0 / rSums.numSamples[varIdx]);
				}
				else
				{
					inelasticRadiance += rSums.radSum[varIdx].I() * (1.0 / rSums.numSamples[varIdx]);
				}
			}
			else
			{
				elasticRadiance += rSums.radSum[varIdx].I() * (1.0 / rSums.numSamples[varIdx]);
				inelasticRadiance += rSums.radSum[varIdx].I() * (1.0 / rSums.numSamples[varIdx]);
			}
		}
	}
	return ok;
}

std::string SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary::PrintSequence(size_t varIdx) const
{
	return SKTRAN_OptimalScatterSequenceManager_Base::PrintSequence(m_var[varIdx], m_varToTotalOrder[varIdx], m_varToElasticRaman[varIdx]);
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary::SetMaxOrder(size_t hardmax)
{
	bool ok = true;

	if (hardmax != m_hardMax)
	{
		ok = ok && SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic::SetMaxOrder(hardmax);
		ok = ok && CalculateNumberOfTerms(m_numVarianceTerms, m_numCovarianceTerms, m_numSubSequences, m_numSubVarianceTerms); // re-calculate the number of covariance terms
		ok = ok && DefineSequences(m_var, m_varToTotalOrder, m_varToRamanOrder, m_varToElasticRaman);
		ok = ok && MapSequences(m_varToSeqIdx, m_seqToVarIdx, m_seqToVarSize);
		ok = ok && DefineCovarianceIndicesLower(m_covToLowerSeqIdx, m_covToUpperSeqIdx, m_seqToCovLowerIdx, m_seqToCovLowerSize, m_covToLowerVarIdx, m_covToUpperVarIdx, m_varToCovLowerIdx, m_varToCovLowerSize);
		ok = ok && DefineCovarianceIndicesUpper(m_covGroupByUpperIdx, m_varToCovUpperIdx, m_varToCovUpperSize, m_seqToCovUpperIdx, m_seqToCovUpperSize);
		ok = ok && DefineSubSequenceIndices(m_subSeq, m_seqToSubSeqIdx);
		ok = ok && DefineSubVarianceIndices(m_subVar, m_seqToSubVarIdx, m_seqToSubVarSize, m_subSeqToSubVarIdx, m_subVarToBufferIdx);
	}

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary::OptimalScatterSequenceIndex(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, double randNum, size_t & ossIdx, bool& minSamplesComplete) const
{
	bool ok = true;

	size_t maxsize_t = -1;
	size_t seqIdx, varIdx, i, n;
	size_t minSeqIdx = maxsize_t, minSamples = maxsize_t;

	if (!minSamplesComplete)
	{
		for (i = 0; i < m_numMaxSequences; i++)
		{
			seqIdx = m_maxSeqIdx[i];
			varIdx = m_seqToVarIdx[seqIdx];
			n = rSums.numSamples[varIdx];
			if (n < m_maxSeqToMinSamples[i])
			{
				ossIdx = seqIdx;
				break;
			}
		}
		if (!(++n < m_maxSeqToMinSamples[i] || ++i < m_numMaxSequences))
		{
			minSamplesComplete = true;
		}
		return ok;
	}

	if (randNum < m_minFractionHigherOrderSum)
	{
		bool done = false;
		for (size_t i = 0; !done && i < m_numMaxSequences; i++)
		{
			if (randNum < m_maxSeqToMinFraction[i])
			{
				ossIdx = m_maxSeqIdx[i];
				done = true;
			}
			else
			{
				randNum -= m_maxSeqToMinFraction[i];
			}
		}
		ok = ok && done;
	}
	else
	{
		size_t offset = m_hardMax > 1 ? 2 : 0;
		randNum = (randNum - m_minFractionHigherOrderSum) / (1.0 - m_minFractionHigherOrderSum); // renormalize randNum
		const std::vector<double>& x = rSums.varContribution;


		double target = randNum * x.back();
		auto it = std::upper_bound(x.begin(), x.end(), target);
		ossIdx = std::distance(x.begin(), it) + offset;
	}
	
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary::SortSamples(const size_t & ossIdx, SKTRAN_MCThreadRadianceLogger::runningSums & rSums) const
{
	// take the values stored in the intermediate buffer and add them to the appropriate running sums

	bool ok = true;

	size_t subVarStartIdx = m_seqToSubVarIdx[ossIdx]; // beginning of the list of all variance terms that were sampled during this iteration
	size_t i, j, uSubVarIdx, uVarIdx, lSubVarIdx, lVarIdx, covUpperIdx, covIdx;	

	//printf("\nVariance term list:\n");
	//for (i = 0; i < m_numVarianceTerms; i++)
	//{
	//	printf("%zd (%s)\n", i, PrintSequence(i, false).c_str());
	//}

	//printf("\nCovariance term list:\n");
	//for (i = 0; i < m_numCovarianceTerms; i++)
	//{
	//	printf("%zd (%s-%s)\n", i, PrintSequence(m_covToUpperVarIdx[i], false).c_str(), PrintSequence(m_covToLowerVarIdx[i], false).c_str());
	//}

	//printf("\nFinal sequence index: %zd (%s)\n", ossIdx, PrintSequence(ossIdx, true).c_str());
	//for (i = 0; i < m_seqToSubVarSize[ossIdx]; i++)
	//{
	//	uSubVarIdx = subVarStartIdx + i;
	//	uVarIdx = m_subVar[uSubVarIdx];
	//	printf("Sub-variance-term index: %zd (%s)\n", uSubVarIdx, PrintSequence(uVarIdx, false).c_str());
	//}

	for (i = 0; i < m_seqToSubVarSize[ossIdx]; i++)
	{
		uSubVarIdx = subVarStartIdx + i;
		uVarIdx = m_subVar[uSubVarIdx];
		rSums.radSum[uVarIdx] += rSums.radBuffer[i];
		rSums.rad2SumVar[uVarIdx] += std::pow(rSums.radBuffer[i].I(), 2.0);
		rSums.numSamples[uVarIdx]++;

		covUpperIdx = m_varToCovUpperIdx[uVarIdx];
		for (j = 0; j < i; j++)
		{
			lSubVarIdx = subVarStartIdx + j;
			lVarIdx = m_subVar[lSubVarIdx];
			covIdx = m_covGroupByUpperIdx[covUpperIdx++];
			rSums.rad2SumCov[covIdx] += rSums.radBuffer[i].I() * rSums.radBuffer[j].I();
			//printf("%zd (%zd-%zd) (%s-%s)\n", covIdx, uVarIdx, lVarIdx, PrintSequence(uVarIdx, false).c_str(), PrintSequence(lVarIdx, false).c_str());
		}
	}

	for (i = 0; i < rSums.radBuffer.size(); i++) rSums.radBuffer[i].SetTo(0.0);

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary::SubmitSample(const size_t ossIdx, const size_t order, const SKTRAN_MCPhoton_Base * mcphoton, std::vector<SKTRAN_MCThreadRadianceLogger::runningSums>& rSums) const
{
	// place a sample (and potentially its elastic counterpart) in the intermediate buffer

	bool ok = true;
	size_t subSeqIdx = m_seqToSubSeqIdx[ossIdx]; // index of the beginning of the group in m_subSeq containing all sub-sequence of the sequence with index ossIdx
	subSeqIdx += std::min(order, m_numDistinctOrders) - 1; // the group is sorted by order; order > m_numDistinctOrders should only occur when ossIdx has no Raman scatters
	size_t subVarIdx = m_subSeqToSubVarIdx[subSeqIdx]; 
	size_t varIdx = m_subVar[subVarIdx];
	size_t bufferIdx = m_subVarToBufferIdx[subVarIdx];

	if (m_varToRamanOrder[varIdx] == 0) // if there are no Raman scatters, we only have one sample to add
	{
		auto pr = mcphoton->photonRadiances(false).cbegin();
		for (auto rs = rSums.begin(); rs != rSums.end(); rs++, pr++)
		{
			rs->radBuffer[bufferIdx] += pr->GetRecentContribVec();
		}
	}
	else // if there is at least one Raman scatter, we have an inelastic and an elastic sample to add
	{
		auto pr = mcphoton->photonRadiances(false).cbegin();
		auto epr = mcphoton->photonRadiances(true).cbegin();
		for (auto rs = rSums.begin(); rs != rSums.end(); rs++, pr++, epr++)
		{
			rs->radBuffer[bufferIdx] += pr->GetRecentContribVec();
			rs->radBuffer[bufferIdx + 1] += epr->GetRecentContribVec();
		}
	}

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary::CalculateVariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& variance) const
{
	bool ok = true;

	for (size_t idx = 0; idx < m_numVarianceTerms; idx++)
	{
		//rSums.varEstimate[idx] = 0 < rSums.numSamples[idx] ? (rSums.rad2SumVar[idx] - (pow(rSums.radSum[idx].I(), 2.0) / ((double)rSums.numSamples[idx]))) * pow((double)rSums.numSamples[idx], -1.0) : 0.0;
		variance[idx] = 0 < rSums.numSamples[idx] ? (rSums.rad2SumVar[idx] - (pow(rSums.radSum[idx].I(), 2.0) / ((double)rSums.numSamples[idx]))) * pow((double)rSums.numSamples[idx], -2.0) : 0.0;
	}

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary::CalculateCovariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& covariance) const
{
	bool ok = true;
	size_t lidx, uidx;

	for (size_t cidx = 0; cidx < m_numCovarianceTerms; cidx++)
	{
		lidx = m_covToLowerVarIdx[cidx];
		uidx = m_covToUpperVarIdx[cidx];
		//rSums.covEstimate[cidx] = 0 < rSums.numSamples[uidx] ? (rSums.rad2SumCov[cidx] - ((rSums.radSum[uidx].I() * rSums.radSum[lidx].I()) / ((double)rSums.numSamples[lidx]))) * pow((double)rSums.numSamples[uidx], -1.0) : 0.0;
		covariance[cidx] = 0 < rSums.numSamples[uidx] ? (rSums.rad2SumCov[cidx] - ((rSums.radSum[uidx].I() * rSums.radSum[lidx].I()) / ((double)rSums.numSamples[lidx]))) * pow((double)rSums.numSamples[lidx], -1.0) * pow((double)rSums.numSamples[uidx], -1.0) : 0.0;
	}

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary::CalculateVarianceContribution(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& contribution) const
{
	size_t offset = m_hardMax > 1 ? 2 : 0;
	std::vector<double> w(m_numSequences - offset); // sum of variance and covariance terms
	std::vector<double>& x = contribution;

	for (size_t seqIdx = offset; seqIdx < m_numSequences; seqIdx++)
	{
		size_t wIdx = seqIdx - offset;
		size_t varIdx = m_seqToVarIdx[seqIdx];
		w[wIdx] = rSums.varEstimate[varIdx] * rSums.varDerivative[varIdx] * rSums.numSamples[varIdx];
		//printf("\nseqIdx %zd %s\n", seqIdx, PrintSequence(seqIdx, true).c_str());
		//printf("varIdx %zd %s\n", varIdx, PrintSequence(varIdx, false).c_str());
		if (m_varToRamanOrder[varIdx] > 0) {
			w[wIdx] += rSums.varEstimate[varIdx + 1] * rSums.varDerivative[varIdx + 1] * rSums.numSamples[varIdx + 1];
			//printf("varIdx %zd %s\n", varIdx + 1, PrintSequence(varIdx + 1, false).c_str());
		}
		for (size_t covIdx = m_seqToCovLowerIdx[seqIdx]; covIdx < m_seqToCovLowerIdx[seqIdx] + m_seqToCovLowerSize[seqIdx]; covIdx++)
		{
			//if (!m_varToElasticRaman[m_covToLowerVarIdx[covIdx]] && !m_varToElasticRaman[m_covToUpperVarIdx[covIdx]]) 
			//{
			w[wIdx] += 2.0 * rSums.covEstimate[covIdx] * rSums.covDerivative[covIdx] * rSums.numSamples[varIdx];
			//}
			//printf("covIdx %zd %s-%s\n", covIdx, PrintSequence(m_covToLowerVarIdx[covIdx], false).c_str(), PrintSequence(m_covToUpperVarIdx[covIdx], false).c_str());
		}
		w[wIdx] = sqrt(std::max(w[wIdx], 0.0));
	}

	double sum = 0.0;
	for (size_t lowIdx = offset; lowIdx < m_numSequences; lowIdx++)
	{
		size_t xIdx = lowIdx - offset;
		size_t cIdx = m_seqToSeqPlusCabannes[lowIdx];
		size_t rIdx = m_seqToSeqPlusRaman[lowIdx];

		sum += w[xIdx];
		if (cIdx < m_numSequences) sum -= w[cIdx - offset];
		if (rIdx < m_numSequences) sum -= w[rIdx - offset];
		x[xIdx] = sum;
	}

	return sum > 0;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary::ConfigureRunningSums(SKTRAN_MCThreadRadianceLogger::runningSums & rSums) const
{
	size_t nvar = m_numVarianceTerms;
	size_t ncov = m_numCovarianceTerms;
	size_t nseq = m_hardMax > 1 ? m_numSequences - 2 : m_numSequences;
	rSums.radBuffer.resize(nvar);
	rSums.radSum.resize(nvar);
	rSums.rad2SumVar.resize(nvar);
	rSums.rad2SumCov.resize(ncov);
	rSums.numSamples.resize(nvar);
	rSums.varEstimate.resize(nvar);
	rSums.covEstimate.resize(ncov);
	rSums.varDerivative.resize(nvar);
	rSums.covDerivative.resize(ncov);
	rSums.varContribution.resize(nseq);
	return true;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary::CalculateMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, SKTRAN_Stokes_NC & measurement) const
{
	bool ok = true;
	skRTStokesVector::SetToZero(measurement);

	for (size_t idx = 0; idx < m_numVarianceTerms; idx++)
	{
		if (0 < rSums.numSamples[idx] && !m_varToElasticRaman[idx])
		{
			measurement += rSums.radSum[idx] * (1.0 / rSums.numSamples[idx]);
		}
	}
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary::CalculateTotalVariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, double & variance) const
{
	bool ok = true;
	variance = 0.0;
	for (size_t idx = 0; idx < m_numVarianceTerms; idx++)
	{
		if (!m_varToElasticRaman[idx]) variance += rSums.varEstimate[idx];
	}
	for (size_t cidx = 0; cidx < m_numCovarianceTerms; cidx++)
	{
		size_t lidx = m_covToLowerVarIdx[cidx];
		size_t uidx = m_covToUpperVarIdx[cidx];
		if (!m_varToElasticRaman[lidx] && !m_varToElasticRaman[uidx]) variance += 2.0 * rSums.covEstimate[cidx];
	}
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary::CalculateTargetMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, double & measurement) const
{
	return CalculateSecondaryMeasurement(rSums, measurement);
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary::CalculateTargetVariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, double & variance, size_t& numSamples) const
{
	bool ok = true;

	ok = ok && CalculateSecondaryVariance(rSums, variance);
	numSamples = rSums.numSamples[0] + rSums.numSamples[1]; // number of first order samples (cabannes + raman) represents total number of samples

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedRing::CalculateDerivative(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& varDerivative, std::vector<double>& covDerivative) const
{
	bool ok = true;

	double inelasticRadiance, elasticRadiance;
	ok = ok && CalculateRadiances(rSums, inelasticRadiance, elasticRadiance);

	double inelasticDerivative = 1.0 / inelasticRadiance;
	double elasticDerivative = -1.0 / elasticRadiance;
	double cabannesDerivative = inelasticDerivative + elasticDerivative;
	double lDeriv, uDeriv;

	size_t uVarIdx = 0, covIdx = 0;
	for (size_t lVarIdx = 0; lVarIdx < m_numVarianceTerms; lVarIdx++)
	{
		lDeriv = m_varToElasticRaman[lVarIdx] ? elasticDerivative : m_varToRamanOrder[lVarIdx] > 0 ? inelasticDerivative : cabannesDerivative;
		varDerivative[lVarIdx] = lDeriv * lDeriv;
		for (size_t j = 0; j < m_varToCovLowerSize[lVarIdx]; j++)
		{
			covIdx = m_varToCovLowerIdx[lVarIdx] + j;
			uVarIdx = m_covToUpperVarIdx[covIdx];
			uDeriv = m_varToElasticRaman[uVarIdx] ? elasticDerivative : m_varToRamanOrder[uVarIdx] > 0 ? inelasticDerivative : cabannesDerivative;
			covDerivative[covIdx] = lDeriv * uDeriv;
		}
	}

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedRing::CalculateSecondaryMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, double & measurement) const
{
	bool ok = true;
	double inelasticRadiance = 0.0, elasticRadiance = 0.0;
	ok = ok && CalculateRadiances(rSums, inelasticRadiance, elasticRadiance);
	measurement = log(inelasticRadiance/elasticRadiance);
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedRing::CalculateSecondaryVariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, double & variance) const
{
	bool ok = true;
	double inelasticRadiance = 0.0, elasticRadiance = 0.0;
	ok = ok && CalculateRadiances(rSums, inelasticRadiance, elasticRadiance);

	double inelasticDerivative = 1.0 / inelasticRadiance;
	double elasticDerivative = -1.0 / elasticRadiance;
	double cabannesDerivative = inelasticDerivative + elasticDerivative;
	   
	size_t cidx = 0, lidx = 0, uidx = 0;
	double dl, du, dv;
	variance = 0.0;

	for (size_t vidx = 0; vidx < m_numVarianceTerms; vidx++)
	{
		dv = m_varToElasticRaman[vidx] ? elasticDerivative : m_varToRamanOrder[vidx] > 0 ? inelasticDerivative : cabannesDerivative;
		variance += rSums.varEstimate[vidx] * dv * dv;
	}

	for (size_t cidx = 0; cidx < m_numCovarianceTerms; cidx++)
	{
		lidx = m_covToLowerVarIdx[cidx];
		uidx = m_covToUpperVarIdx[cidx];
		dl = m_varToElasticRaman[lidx] ? elasticDerivative : m_varToRamanOrder[lidx] > 0 ? inelasticDerivative : cabannesDerivative;
		du = m_varToElasticRaman[uidx] ? elasticDerivative : m_varToRamanOrder[uidx] > 0 ? inelasticDerivative : cabannesDerivative;
		variance += rSums.covEstimate[cidx] * dl * du;
	}

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedElastic::CalculateDerivative(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& varDerivative, std::vector<double>& covDerivative) const
{
	bool ok = true;
	size_t uVarIdx;
	std::fill(varDerivative.begin(), varDerivative.end(), 0.0);
	std::fill(covDerivative.begin(), covDerivative.end(), 0.0);
	for (size_t lVarIdx = 0; lVarIdx < m_numVarianceTerms; lVarIdx++)
	{
		if (m_varToElasticRaman[lVarIdx] || m_varToRamanOrder[lVarIdx] == 0)
		{
			varDerivative[lVarIdx] = 1.0;
			for (size_t covIdx = m_varToCovLowerIdx[lVarIdx]; covIdx < m_varToCovLowerIdx[lVarIdx] + m_varToCovLowerSize[lVarIdx]; covIdx++)
			{
				uVarIdx = m_covToUpperVarIdx[covIdx];
				if (m_varToElasticRaman[uVarIdx] || m_varToRamanOrder[uVarIdx] == 0)
				{
					covDerivative[covIdx] = 1.0;
				}
			}
		}
	}
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedElastic::CalculateSecondaryMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, double & measurement) const
{
	bool ok = true;
	double inelasticRadiance = 0.0, elasticRadiance = 0.0;
	ok = ok && CalculateRadiances(rSums, inelasticRadiance, elasticRadiance);
	measurement = elasticRadiance;
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedElastic::CalculateSecondaryVariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, double & variance) const
{
	bool ok = true;
	double inelasticRadiance = 0.0, elasticRadiance = 0.0;
	ok = ok && CalculateRadiances(rSums, inelasticRadiance, elasticRadiance);

	size_t cidx = 0, lidx = 0, uidx = 0;
	variance = 0.0;

	for (size_t vidx = 0; vidx < m_numCovarianceTerms; vidx++)
	{
		if (m_varToRamanOrder[vidx] == 0 || m_varToElasticRaman[vidx])
		{
			variance += rSums.varEstimate[vidx];
		}
	}

	for (size_t cidx = 0; cidx < m_numCovarianceTerms; cidx++)
	{
		lidx = m_covToLowerVarIdx[cidx]; 
		uidx = m_covToUpperVarIdx[uidx];
		if ((m_varToRamanOrder[lidx] == 0 || m_varToElasticRaman[lidx]) && (m_varToRamanOrder[uidx] == 0 || m_varToElasticRaman[uidx])) variance += rSums.covEstimate[cidx];
	}
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedFillingIn::CalculateDerivative(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, std::vector<double>& varDerivative, std::vector<double>& covDerivative) const
{
	bool ok = true;
	double inelasticRadiance = 0.0, elasticRadiance = 0.0;
	ok = ok && CalculateRadiances(rSums, inelasticRadiance, elasticRadiance);

	double inelasticDerivative = 1.0 / elasticRadiance;
	double elasticDerivative = - inelasticDerivative * inelasticDerivative * inelasticRadiance;
	double cabannesDerivative = inelasticDerivative + elasticDerivative;
	double lDeriv, uDeriv;

	size_t uVarIdx = 0, covIdx = 0;
	for (size_t lVarIdx = 0; lVarIdx < m_numVarianceTerms; lVarIdx++)
	{
		lDeriv = m_varToElasticRaman[lVarIdx] ? elasticDerivative : m_varToRamanOrder[lVarIdx] > 0 ? inelasticDerivative : cabannesDerivative;
		varDerivative[lVarIdx] = lDeriv * lDeriv;
		for (size_t j = 0; j < m_varToCovLowerSize[lVarIdx]; j++)
		{
			covIdx = m_varToCovLowerIdx[lVarIdx] + j;
			uVarIdx = m_covToUpperVarIdx[covIdx];
			uDeriv = m_varToElasticRaman[uVarIdx] ? elasticDerivative : m_varToRamanOrder[uVarIdx] > 0 ? inelasticDerivative : cabannesDerivative;
			covDerivative[covIdx] = lDeriv * uDeriv;
		}
	}

	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedFillingIn::CalculateSecondaryMeasurement(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, double & measurement) const
{
	bool ok = true;
	double inelasticRadiance = 0.0, elasticRadiance = 0.0;
	ok = ok && CalculateRadiances(rSums, inelasticRadiance, elasticRadiance);
	measurement = (inelasticRadiance - elasticRadiance) / elasticRadiance;
	return ok;
}

bool SKTRAN_OptimalScatterSequenceManager_OptimizedFillingIn::CalculateSecondaryVariance(const SKTRAN_MCThreadRadianceLogger::runningSums & rSums, double & variance) const
{
	bool ok = true;
	double inelasticRadiance = 0.0, elasticRadiance = 0.0;
	ok = ok && CalculateRadiances(rSums, inelasticRadiance, elasticRadiance);

	double inelasticDerivative = 1.0 / inelasticRadiance;
	double elasticDerivative = -1.0 / elasticRadiance;
	double cabannesDerivative = inelasticDerivative + elasticDerivative;

	size_t cidx = 0, lidx = 0, uidx = 0;
	double dl, du, dv;
	variance = 0.0;

	for (size_t vidx = 0; vidx < m_numVarianceTerms; vidx++)
	{
		dv = m_varToElasticRaman[vidx] ? elasticDerivative : m_varToRamanOrder[vidx] > 0 ? inelasticDerivative : cabannesDerivative;
		variance += rSums.varEstimate[vidx] * dv * dv;
	}

	for (size_t cidx = 0; cidx < m_numCovarianceTerms; cidx++)
	{
		lidx = m_covToLowerVarIdx[cidx];
		uidx = m_covToUpperVarIdx[cidx];
		dl = m_varToElasticRaman[lidx] ? elasticDerivative : m_varToRamanOrder[lidx] > 0 ? inelasticDerivative : cabannesDerivative;
		du = m_varToElasticRaman[uidx] ? elasticDerivative : m_varToRamanOrder[uidx] > 0 ? inelasticDerivative : cabannesDerivative;
		variance += rSums.covEstimate[cidx] * dl * du;
	}

	return ok;
}

