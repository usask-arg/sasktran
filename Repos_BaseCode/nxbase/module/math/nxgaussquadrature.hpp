
/*---------------------------------------------------------------------------
 *'					nxGaussQuadrature::Integrate		2003-9-26
 *	This integrates over the specified range using a function object
 *	that has a function operator () of double ()( double x)
 *	The function object evaluates the mathematical function at the value
 *	X and returns the value Y.
 *-------------------------------------------------------------------------*/

template <class FuncObj>
double	nxGaussQuadrature<FuncObj>::Integrate( FuncObj function )
{
	nxArrayIter<double> wi;
	nxArrayIter<double> xi;
	nxArrayIter<double>	we;
	double				sum;
	
	NXASSERT(m_N > 0);
	CheckDirty();
	if (m_N < 1)
	{
		nxLog::Record(NXLOG_WARNING, "nxGaussQuadrature::Integrate, the quadrature order (%d) is not valid, return 1.0E300 as integral", (int)m_N);
		sum = 1.0E300;
	}
	else
	{
		xi  = m_x.begin();
		wi  = m_weights.begin();
		we  = m_weights.end();
		sum = 0;
		while( wi != we)
		{
			sum += (*wi)*function(*xi);
			++wi;
			++xi;
		}
	}
	return sum;
}
/*---------------------------------------------------------------------------
 *'					nxTrapezoidalQuadrature::Integrate		2003-9-26
 *	This integrates over the specified range using a function object
 *	that has a function operator () of double ()( double x)
 *	The function object evaluates the mathematical function at the value
 *	X and returns the value Y.
 *-------------------------------------------------------------------------*/

template <class FuncObj>
double	nxTrapezoidalQuadrature<FuncObj>::Integrate( FuncObj function )
{
	double				h;
	double				xi;
	double				sum;
	nxBOOL				ok;
	int					N1;
	
	ok = ( m_N >= 2);
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING, "nxTrapezoidalQuadratureBase::Integrate, There are not enough points defined fro trapezoidal integration (%d), the integral is returned as 0", (int)m_N);
		sum = 0;
	}
	else
	{
		N1  = m_N - 1;
		h   = (m_x2-m_x1)/N1;
		sum = 0.5*function(m_x1);
		for (int i = 1; i < N1; i++)
		{
			xi = m_x1 + i*h;
			sum += function(xi);
		}
		sum += 0.5*function(m_x2);
		sum *= h;
	}
	return sum;
}
