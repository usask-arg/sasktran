#if !defined(NXBASE_NXTRACE_H)
#define NXBASE_NXTRACE_H 1
/*-----------------------------------------------------------------------------
 * 					nxTrace														*/
/**	\ingroup system_debug
 *  A source code trace debug function.  Normally it accessed via the NXTRACE
 *	and NXDEBUG macros. NXTRACE evaluates to nothing if NXDEBUG is not defined
 *
 *	\param formatstr
 *		string used to encode variable parameter list
**/
/*----------------------------------------------------------------------------*/

extern void nxTrace( const char* lpszFormat, ...);

#if defined(NXTRACE_AS_PRINT)			/* Do we explicitly want NXTRACE to be printf */
	#define NXTRACE(e) /*printf e*/			/* yup */
#elif defined(NXDEBUG)					/* otherwise are we in a debug mode */
	#define NXTRACE(e) nxTrace e		/* so use the implicit debug trace version */
#else									/* otherwise do nothing */
	#define NXTRACE(e)					/* so print a blank */
#endif

#define NXTRACE_ONCEONLY(firsttime,e) \
static bool firsttime=true;\
if (firsttime)\
{\
NXTRACE(e);\
firsttime = false;\
}\


#if defined(NXDEBUG)					/* if we are debugging */
	#if defined (_MSC_VER)				/* and if this is the MSC version */
		#define NXASSERT(expr)  do {if (!(expr) && (1 == _CrtDbgReport(_CRT_ASSERT, __FILE__, __LINE__, NULL, #expr))) _CrtDbgBreak(); } while (0)
	#else
		#define NXASSERT(expr) \
		if (!(expr) )\
		{\
			printf("**** ASSERTION FORCED in file (%s) at line %d\n", (const char*)__FILE__, (int)__LINE__);\
 			exit(0);\
		}
	#endif
#else									/* otherwise we are not debugging */
	#define NXASSERT(expr) 				/* so do nothing */
#endif

#endif

