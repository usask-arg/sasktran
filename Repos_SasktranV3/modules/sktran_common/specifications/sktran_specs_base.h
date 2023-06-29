
/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Base		2013-9-27*/
/**  @ingroup specs
 *	Base class for specifications set by the user. The engine actually uses the
 *	information in derived User Specification classes to create Internal Specifications
 *	which are internally stored by the engine. The main issue handled is that the User based
 *	specifications has no requirements on the end-user to maintain the lifetime of the object as all of the
 *	information is essentially cloned into the internal specs. Early versions of the engine had problems
 *	when users trashed the specifications not realising that the engine was still using the data.
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_SpecsUser_Base
{
	public:
							SKTRAN_SpecsUser_Base() {};
		virtual			   ~SKTRAN_SpecsUser_Base() {};
};


/*-----------------------------------------------------------------------------
 *					class SKTRAN_SpecsInternal_Base		2013-9-27*/
/**  @ingroup specs
 *	Base class for storing internal engine specifications. This class is
 *	normally created by each engine to store all of the relevant specifications
 *	for the calculation. The object has lifetime management provided by nxUnknown.**/
/*---------------------------------------------------------------------------*/

class SKTRAN_SpecsInternal_Base : public nxUnknown
{
	public:
							SKTRAN_SpecsInternal_Base() {};
		virtual			   ~SKTRAN_SpecsInternal_Base() {};
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_Specifications_Base		2013-10-15*/
/**  @ingroup specs
 * Base class for Sasktran Version 3 Specifications
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_Specifications_Base : public SKTRAN_SpecsUser_Base
{
	public:
							SKTRAN_Specifications_Base() {};
		virtual			   ~SKTRAN_Specifications_Base() {};
};

