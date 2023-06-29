
class CodeTimer
{
	private:
		nxTimeStamp m_t1;
		nxTimeStamp m_t2;
	public:
		CodeTimer( ){m_t1.MJD(0.0); m_t2.MJD(0.0);};
		~CodeTimer( ){ };

		void Start( )	{ m_t1.FromSystem( ); };
		void Stop( )	{ m_t2.FromSystem( ); };
		void Report( )	{ printf( "\t\tCodeTimer:: %lf seconds\n", (m_t2-m_t1).MJD( )*24.0*60.0*60.0 ); };
};

