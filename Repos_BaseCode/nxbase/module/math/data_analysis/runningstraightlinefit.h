
/*-----------------------------------------------------------------------------
 *					class RunningStraightLineFit					2004-11-23*/
/** \ingroup Math_fitting
 *	Fits a straight line to an array of (x,y) data pairs.  Class was originally
 *	designed to fit Universal time to the Freja Satellite time word.  It uses a
 *	running average so it only fits to the last N points where N is defined by member
 *	m_maxpoints.
 **/
/*---------------------------------------------------------------------------*/

class  RunningStraightLineFit
{
	private:
		STL(list)<double>	m_x;			// List of individual x coords
		STL(list)<double>	m_y;			// List of individual y coords
		double				m_origx;		// Shifted X origin
		double				m_origy;		// Shifted Y origin
		double				m_sx;			// Sum of x
		double				m_sy;			// Sum of y
		double				m_sx2;			// Sum of x2
		double				m_syx;			// Sum of yx
		double				m_sy2;			// Sum of y2
		int					m_maxpoints;	// Maximum points in fit

	private:
		nxBOOL				DoFit( double* gradient, double *intercept );
		void				RemoveOldestPoint();

	public:
							RunningStraightLineFit();
		void				Reset();
		void				Insert( double x, double y );
		double				Interpolate( double x );
};


/*<---------------------------------------------------------------------------
 *'						class nxRunningStatistic
 *
 *>--------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------
 *					nxRunningStatistic		2004-11-23*/
/**	\ingroup Math_fitting
 *	Fits a straight line to an array of (x,y) data pairs.  Class is designed
 *	to fit Universal time to satellite time word.  It uses a running average
 *	so it only fits to the last N points where N is defined by member
 *	m_maxpoints.
 **/
/*---------------------------------------------------------------------------*/

class  nxRunningStatistic
{
	private:
		typedef STL(list)<double>::iterator iterator;
	private:
		STL(list)<double>	m_x;			// List of individual x coords
		double				m_sx;			// Sum of x
		double				m_sx2;			// Sum of x2
		double				m_minx;
		double				m_maxx;
		int					m_maxpoints;

	private:
		nxBOOL				DoFit( double* gradient, double *intercept );
		void				RemoveOldestPoint();

	public:
							nxRunningStatistic();
		void				Reset();
		void				Insert( double x );
		void				SetMaxPoints( int npts ) { Reset(); m_maxpoints = npts;}
		double				CurrentValue() { return m_x.back();}
		double				MinValue()	   { return m_minx; }
		double				MaxValue()     { return m_maxx; }
		double				SD();
		double				Average();
};

