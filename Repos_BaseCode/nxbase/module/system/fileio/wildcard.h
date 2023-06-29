
/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/

/*----------------------------------------------------------------------------
 *					class nxWildcard										*/
/**	\ingroup system_strings
 *   A class that may be used to compare a string against a pre-determined wildcard
 *   specification.  Used for selecting specific files from a list of files etc.
 *   The wildcard codes interprets both '*' and '?' in the normal manner. For example
 *
 *	\code
 *   nxWildcard  spec( "*345.*?0")				\\ Setup a wildcard specification
 *   spec.Matches( "P012345.TN0");				\\ returns TRUE
 *	\endcode
**/
/*---------------------------------------------------------------------------*/

class  nxWildcard
{

   private:
      char*			wildspec;			// The  current wild card specification

   public:
					nxWildcard( const char *wildspecstr = NULL );
				   ~nxWildcard();
      void			SetSpec( const char *wildspecstr );
      nxBOOL		Match  ( const char *ThisString );
};




