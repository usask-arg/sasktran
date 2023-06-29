#include <nxbase_core.h>
#include <nxbase_math.h>
#include <netcdf.h>
class nxNetcdfGroups;
class nxNetcdfGroup;
class nxNetcdfVar;

/*-----------------------------------------------------------------------------
 *					class nxNetcdfEntity		2011-7-15*/
/** A base class for NetCDF variable, dimensions, group or attribute. The NetCDF
 *  file appears to be a hierarchical system where an entity always has a parent
 *	and may have a number of children. The base class stores basic handle info
 *	about the parent and this object. Specific classes store information about
 *	children.
 **/
/*---------------------------------------------------------------------------*/

class nxNetcdfEntity
{
	private:
		nxNetcdfEntity*			m_parent;
		int						m_ncid;
		std::string				m_name;
	

	public:
								nxNetcdfEntity	() { m_ncid = 0; m_parent = NULL;}
		virtual				   ~nxNetcdfEntity	() {};
		bool					DeepCopy		( const nxNetcdfEntity& other);

		void					SetNCID			(int ncid)							{ m_ncid   = ncid;}
		void					SetName			( const char * name )				{ m_name   = name;}
		void					SetParent		( nxNetcdfEntity* parent )			{ m_parent = parent;}
		int						NCID			() const							{ return m_ncid;}
		int						ParentNCID		() const							{ return (m_parent != NULL) ? m_parent->NCID() : NC_GLOBAL;}
		const nxNetcdfGroup*	ParentGroup		() const;
		const nxNetcdfGroup*	RootGroup		() const;
		const std::string&		Name			() const							{ return m_name;}
		bool					operator <		( const nxNetcdfEntity& other ) const { return m_name <  other.m_name;}
		bool					operator ==		( const nxNetcdfEntity& other ) const { return m_name == other.m_name;}

	public:
		virtual int				Handler_Inq_NumIDs		( int /*ncid*/, int* /*numids*/, int* /*idarray*/) {return -1;}
		virtual int				Handler_Inq_Names		( int /*ncid*/, int  /*numids*/, char* /*name*/)  {return -1;}
 };

/*-----------------------------------------------------------------------------
 *					class nxNetcdfEntityArray		2011-7-15*/
/** A templated class that holds an array of NetCDFEntity based objects,
    either, Dimensions, Variables Attributes or Groups.  This class is used
	as a container for all the children of a given type belonging to a given 
	entity 
**/
/*---------------------------------------------------------------------------*/

template < class nxNetcdfEntityType>
class nxNetcdfEntityArray
{
	private:                  std::map< std::string, nxNetcdfEntityType>					m_array;
	public:  typedef typename std::map< std::string, nxNetcdfEntityType>::iterator			iterator;
	public:  typedef typename std::map< std::string, nxNetcdfEntityType>::const_iterator	const_iterator;
	public:  typedef typename std::map< std::string, nxNetcdfEntityType>::value_type		value_type;


	public:
		const nxNetcdfEntityType*						At				( const char* name ) const;
		const nxNetcdfEntityType*						Find			( int ncid ) const;
//		const nxNetcdfEntityType*						AddBlankEntry	( int ncid, const char* name) const;
		std::map< std::string, nxNetcdfEntityType>&		Array			()						{return m_array;} 
		void											Clear			()						{ m_array.clear();}
		bool											Load			( nxNetcdfGroup* groupparent);
		bool											LoadAttributes	( nxNetcdfVar*   variableparent );
};

/*-----------------------------------------------------------------------------
 *				struct nxNetcdfDimStruct	2011-7-19*/
/** **/
/*---------------------------------------------------------------------------*/

class nxNetcdfDimStruct
{
	public:
		std::string				name;
		std::string				units;
		std::vector<double>		dimensionpoints;

};

/*-----------------------------------------------------------------------------
 *					class nxNetcdfDim		2011-7-15*/
/** **/
/*---------------------------------------------------------------------------*/

class nxNetcdfDim : public nxNetcdfEntity
{
	public:
								nxNetcdfDim		()	{};
		virtual				   ~nxNetcdfDim		()	{};
		size_t					Length			() const;
		bool					GetDimensionInfo( nxNetcdfDimStruct* info ) const;

	public:
		virtual int				Handler_Inq_NumIDs		( int ncid, int* numids, int*idarray) { return nc_inq_dimids ( ncid, numids, idarray, 0 );}
		virtual int				Handler_Inq_Names		( int ncid, int  numids, char* name)  { return nc_inq_dimname( ncid, numids, name       );}

};


/*-----------------------------------------------------------------------------
 *					class nxNetcdfAtt		2011-7-15*/
/** **/
/*---------------------------------------------------------------------------*/

class nxNetcdfAtt : public nxNetcdfEntity
{
	public:
								nxNetcdfAtt(){};
		virtual				   ~nxNetcdfAtt(){};

};

class nxNetcdfVar;

/*---------------------------------------------------------------------------
 *                  Class nxNetcdfHyperSlabDefn                   2019-10-03 */
/** Used to define hyperslabs 
**/
/*---------------------------------------------------------------------------*/

class nxNetcdfHyperSlabDefn
{
	private:
		const nxNetcdfVar*			m_variable;
		nx1dArray<size_t>			m_rankspecs;
		std::vector<size_t>			m_start;				// THese are stored in reverse order so they come out in the correct order after netcdf
		std::vector<size_t>			m_count;
		std::vector<size_t>			m_stride;

	private:
		void						ClearSlabRange();
	public:
									nxNetcdfHyperSlabDefn		();
		bool						AttachToVariable			( const nxNetcdfVar* variable );
		bool						DefineRanges				( size_t npoints, size_t* start, size_t* count, size_t* stride= nullptr);
		bool						AllocateArrayToHoldHyperSlab( nxArrayLinear<double>* slabarray) const;
		const size_t*				Start						() const  { return &m_start.front();}
		const size_t*				Count						() const  { return &m_count.front();}
		const size_t*				Stride						() const  { return &m_stride.front();}
};

/*-----------------------------------------------------------------------------
 *					class nxNetcdfVar		2011-7-15*/
/** **/
/*---------------------------------------------------------------------------*/

class nxNetcdfVar : public nxNetcdfEntity
{
	private:
		nxNetcdfEntityArray<nxNetcdfAtt>		m_atts;

	private:
		bool								AdjustForScaleAndOffset	( nxArrayLinear<double>* vararray, double usermissingvalue ) const;
		bool								GetVarAttributes		( nxNetcdfDimStruct* info ) const;


	public:
											nxNetcdfVar(){}
		virtual							   ~nxNetcdfVar(){}
		bool								LoadAttributes			( );
		bool								LoadRankSpecs			( nx1dArray<size_t>*		     rankarray ) const;
		bool								LoadData				( nxArrayLinear<double>*	     vararray, double usermissingvalue   )const;
		bool								LoadDataSlice			( nxArrayLinear<double>*	     vararray, const nxNetcdfHyperSlabDefn& hyperslab, double usermissingvalue   )const;
		bool								LoadDimInfo				( std::vector<nxNetcdfDimStruct>*  dimsinfo  )const;
		bool								LoadVariable			( nxArrayLinear<double>* vararray, std::vector<nxNetcdfDimStruct>* dimsinfo, double usermissingvalue  )const;
		bool								AttributeString			( const char* attributename, std::string* value ) const;
		bool								AttributeDouble			( const char* attributename, double* value ) const;

	public:
		virtual int							Handler_Inq_NumIDs		( int ncid, int* numids, int*idarray) { return nc_inq_varids ( ncid, numids, idarray );}
		virtual int							Handler_Inq_Names		( int ncid, int  numids, char* name)  { return nc_inq_varname( ncid, numids, name  );}

};


class nxNetcdfGroups;

/*-----------------------------------------------------------------------------
 *					class nxNetcdfGroup		2011-7-15*/
/** **/
/*---------------------------------------------------------------------------*/

class nxNetcdfGroup : public nxNetcdfEntity
{
	private:
		nxNetcdfGroups*						m_groups;
		nxNetcdfEntityArray<nxNetcdfVar>		m_vars;
		nxNetcdfEntityArray<nxNetcdfDim>		m_dims;

	private:
		bool								LoadVariableAttributes();
		bool								LoadSubGroups	      ();


	public:
											nxNetcdfGroup();
											nxNetcdfGroup( const nxNetcdfGroup& other );
		virtual							   ~nxNetcdfGroup();
		bool								Load				();
		const nxNetcdfGroup*						At			( const char* groupname) const;
		const nxNetcdfVar*					VarAt				( const char* varname, bool recursegroups) const;
		const nxNetcdfDim*					DimAt				( const char* dimname ) const	{ return m_dims.At(dimname);}
		const nxNetcdfDim*					DimAt				( int ncid ) const				{ return m_dims.Find(ncid);}
		void								ReleaseResources	();
		const nxNetcdfVar*					CoordinateVariable	( const char* varname, bool searchfromroot = true ) const;

	public:
		virtual int							Handler_Inq_NumIDs		( int ncid, int* numids, int*idarray) { return nc_inq_grps   (  ncid, numids, idarray);}
		virtual int							Handler_Inq_Names		( int /*ncid*/, int  id,     char* name)  { return nc_inq_grpname(  id,   name   );}

};


/*-----------------------------------------------------------------------------
 *					class nxNetcdfGroups		2011-7-15*/
/** **/
/*---------------------------------------------------------------------------*/

class nxNetcdfGroups : public nxNetcdfEntity
{
	private:
		nxNetcdfEntityArray<nxNetcdfGroup>		m_groups;
	
	public:
											nxNetcdfGroups		(){};
										   ~nxNetcdfGroups		(){};
		void								ReleaseResources	();
		bool								LoadGroups_OLD		( nxNetcdfGroup* parent );
		bool								Load				( nxNetcdfGroup* parent )					{ return m_groups.Load( parent );}
		nxNetcdfEntityArray<nxNetcdfGroup>&	Groups				()											{ return m_groups;}
		const nxNetcdfGroup*				At					( const char* name )						{ return m_groups.At(name); }
		const nxNetcdfVar*					VarAt				( const char* varname, bool recursegroups );

};


/*-----------------------------------------------------------------------------
 *					class nxNetcdfFile		2011-7-15*/
/** **/
/*---------------------------------------------------------------------------*/

class nxNetcdfFile : public nxNetcdfEntity
{
	private:
		nxNetcdfGroup			m_topgroup;
		bool					m_isopen;

	public:
								nxNetcdfFile();
		virtual				   ~nxNetcdfFile();
		void					Close		();
		bool					OpenRead	( const char* filename );
		bool					IsOpen		() const { return m_isopen;}
		const nxNetcdfGroup*	At			( const char* name )			{ return m_topgroup.At(name);}
		const nxNetcdfVar*		VarAt		( const char* name )			{ return m_topgroup.VarAt(name, true);}
};

#include "netcdfarray.hpp"
