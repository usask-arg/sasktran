#include "nxbase_math.h"


//---------------------------------------------------------------------------
//						RunninsStraightLineFit::Constructor
//---------------------------------------------------------------------------

RunningStraightLineFit::RunningStraightLineFit()
{
	Reset();
}

//---------------------------------------------------------------------------
//						RunninsStraightLineFit::RemoveOldestPoint
//	Get the oldest point and remove it form the ongoing summations.
//---------------------------------------------------------------------------

void RunningStraightLineFit::RemoveOldestPoint()
{
	double x = *m_x.begin();
	double y = *m_y.begin();

	m_sx  -= x;
	m_sy  -= y;
	m_sx2 -= (x*x);
	m_syx -= (y*x);
	m_sy2 -= (y*y);
	m_x.pop_front();
	m_y.pop_front();
}

//---------------------------------------------------------------------------
//						RunninsStraightLineFit::Reset
//	Clear out any stored data points and reset the summation variables
//---------------------------------------------------------------------------

void RunningStraightLineFit::Reset()
{
	m_sx  = 0;
	m_sy  = 0;
	m_sx2 = 0;
	m_syx = 0;
	m_sy2 = 0;
	m_x.erase( m_x.begin(), m_x.end() );
	m_y.erase( m_y.begin(), m_y.end() );
	m_maxpoints = 1000;
	m_origy = 0;
	m_origx = 0;
}


//---------------------------------------------------------------------------
//						RunninsStraightLineFit::Insert
//	insert an (x,y) data pair into the object.  If we have more than the 
//	permitted maximum, then also remove the oldest point.
//---------------------------------------------------------------------------

void RunningStraightLineFit::Insert( double x, double y )
{

	if ( m_x.empty())		// If the array is currently empty
	{						// then create the shift in origin
		m_origx = x;		// as it helps increase accuracy quite significantly
		m_origy = y;		// especially calculating chi-squared
	}
	x -= m_origx;
	y -= m_origy;

	m_sx  += x;
	m_sy  += y;
	m_sx2 += (x*x);
	m_syx += (y*x);
	m_sy2 += (y*y);

	m_x.push_back(x);
	m_y.push_back(y);

	if ((int)m_x.size() > m_maxpoints) RemoveOldestPoint();
}

//---------------------------------------------------------------------------
//						RunninsStraightLineFit::DoFit
//	Fit a straight line to the current data stored in memory.  All the
//	summations are done as points are added so this is a quick calculation
//---------------------------------------------------------------------------

nxBOOL RunningStraightLineFit::DoFit( double* gradient, double *intercept )
{
	double N = (double)m_x.size();
	nxBOOL ok;
	double m,c; //, chi;

	ok = (N > 1 );
	if (ok)
	{
 		m  = (N*m_syx - m_sy*m_sx)/( N*m_sx2 - m_sx*m_sx);						// Do the least squares fir
		c  = (m_sy - m*m_sx)/N;													// on the shifted origin
//		chi = (m_sy2 + c*(N*c-2*m_sy)  + 2*m*(c*m_sx -m_syx)  + m*m*m_sx2)/(N-1);	// Get chi squared

/*		STL(list)<double>::iterator px; 
		STL(list)<double>::iterator py = m_y.begin();
		double chi2 = 0;
		double yt;
		double dy, lx, ly;
		for ( px = m_x.begin();!(px == m_x.end()); px++)
		{
			lx = *px;
			ly = *py;
			
			yt = m*(*px) + c;
			dy = (*py)-yt;
			chi2 += dy*dy;
			py++;
		}
		chi2 /= N;
*/		
		c += m_origy - (m*m_origx);											// Then correct for the shift in origin

//		nxLog::Record(NXLOG_INFO, "DoFIt, Chi = %lf  %lf", (double)chi, (double) chi2 );
	}
	else
	{
		m = 0;
		c = 0;
	}

	*gradient = m;
	*intercept = c;
	return ok;
}

//---------------------------------------------------------------------------
//						RunninsStraightLineFit::Interpolate
//	Given a value of x estimate a value of Y.
//---------------------------------------------------------------------------

double RunningStraightLineFit::Interpolate( double x )
{
	double m;
	double c;
	double y;

	DoFit( &m, & c );
	y = m*x + c;
	return y;
}
