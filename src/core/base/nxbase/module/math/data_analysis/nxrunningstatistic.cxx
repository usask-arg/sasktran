#include "nxbase_math.h"



//---------------------------------------------------------------------------
//						RunninsStraightLineFit::Constructor
//---------------------------------------------------------------------------

nxRunningStatistic::nxRunningStatistic()
{
	m_maxpoints = 1000;
	Reset();
}

//---------------------------------------------------------------------------
//						RunninsStraightLineFit::RemoveOldestPoint
//	Get the oldest point and remove it form the ongoing summations.
//---------------------------------------------------------------------------

void nxRunningStatistic::RemoveOldestPoint()
{
	double x = *m_x.begin();
	m_sx  -= x;
	m_sx2 -= (x*x);
	m_x.pop_front();
	if (x >= m_maxx)
	{
		iterator ptr1 = STL(max_element)(m_x.begin(), m_x.end() );
		m_maxx = *ptr1;
	}
	if (x <= m_minx)
	{
		iterator ptr2 = STL(min_element)(m_x.begin(), m_x.end() );
		m_minx = *ptr2;
	}
}

//---------------------------------------------------------------------------
//						RunninsStraightLineFit::Reset
//	Clear out any stored data points and reset the summation variables
//---------------------------------------------------------------------------

void nxRunningStatistic::Reset()
{
	m_sx  = 0;
	m_sx2 = 0;
	m_maxx = -1.0E100;
	m_minx =  1.0E100;
	m_x.erase( m_x.begin(), m_x.end() );
}


//---------------------------------------------------------------------------
//						RunninsStraightLineFit::Insert
//	insert an (x,y) data pair into the object.  If we have more than the 
//	permitted maximum, then also remove the oldest point.
//---------------------------------------------------------------------------

void nxRunningStatistic::Insert( double x )
{
	m_sx  += x;
	m_sx2 += (x*x);
	m_x.push_back(x);
	if (x < m_minx ) m_minx = x;
	if (x > m_maxx ) m_maxx = x;
	if ((int)m_x.size() > m_maxpoints) RemoveOldestPoint();
}

double nxRunningStatistic::SD()
{
	double sd;
	int n = (int)m_x.size();

	if (n > 1)
	{
		sd = (m_sx2 - m_sx*m_sx/n )/( n-1);
	}
	else
	{
		sd = 0.0;
	}
	return sd;
}

double nxRunningStatistic::Average()
{
	double avg;
	int n = (int)m_x.size();

	if (n > 0 ) avg = m_sx/n;
	else        avg = 0.0;
	return avg;
}
