#if !defined( NXBASE_NXUNKNOWN_H)
#define NXBASE_NXUNKNOWN_H 1

/*---------------------------------------------------------------------------
 *					class nxUnknown                                2003-3-23*/
/**	\ingroup system_com
 *	A local class very similar to Microsofts pervasive IUnknown.  This class
 *	primarily provides lieftime management of objects via reference counting.
 *
 *	\p Multithreading 
 *	I have added added multithreading support for Windows implementation.
 *	The AddRef and Release are now mutlthread safe.  In addition I have added
 *	AddRef and Release methods for the as constant instant of nxUnknown.  While
 *	this does modify the constant object, contrary to the const declarator 
 *	it does enable const instances of derived classed to have their lifetime
 *	controlled using reference counting, which is very valuable.
 *
 *	\p Todo
 *	I still need to implement multithreading support on Linux using the
 *	Boost API or something similar.
 *	
 *-------------------------------------------------------------------------*/

class nxUnknown
{
	private:
		LONG			m_cref;
		nxBOOL			m_isstatic;

	public:
						nxUnknown();
		virtual		   ~nxUnknown();
		LONG			AddRef();
		LONG			Release();
		LONG			AddRef() const;		//!< Useful for reference counting on objects that are constant and should not be modified
		LONG			Release() const;
		void			SetStatic();
};

#endif

