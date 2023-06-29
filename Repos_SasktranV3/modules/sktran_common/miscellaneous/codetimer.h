/*-----------------------------------------------------------------------------
 *					SKTRAN_CodeTimer		2008-3-26*/
/** @ingroup diagnostics
 *	Class that can be used to measure execution times.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_CodeTimer
{
	private:
		size_t			m_t1;
		size_t			m_t2;

	public:
						SKTRAN_CodeTimer( );
		void			Start		( );
		void			Stop		( );
		void			Report		(const char * message );
};

