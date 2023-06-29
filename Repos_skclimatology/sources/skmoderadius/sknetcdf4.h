#include "nxnetcdfio.h"


/*-----------------------------------------------------------------------------
 *					sknetcdfVariable		2011-7-25*/
/** \internal
**/
/*---------------------------------------------------------------------------*/

class sknetcdfVariable
{
	private:
		std::string							m_name;
		nxArrayLinear<double>*				m_data;
		std::vector<nxNetcdfDimStruct>		m_diminfo;		// Information of the dimensions + 1 extra for info on the variable

	private:
		void								ReleaseResources	();


	public:
											sknetcdfVariable				();
											sknetcdfVariable				( const char* name);
											sknetcdfVariable				( const sknetcdfVariable& other );
										   ~sknetcdfVariable				();
		bool								LoadFromNetCDF				( nxNetcdfFile& ncfile, double missingvalue );
		const nxArrayLinear<double>*		Data						() const { return m_data;}
		const nxNetcdfDimStruct*			Dimension					( const char* name) const;
		const nxNetcdfDimStruct*			Dimension					( size_t idx ) const;
		bool								operator ==					( const sknetcdfVariable& other ) const		{ return m_name == other.m_name;}
		bool								operator <					( const sknetcdfVariable& other ) const		{ return m_name <  other.m_name;}
};


/*-----------------------------------------------------------------------------
 *					sknetcdf_IONetCDF		2005-3-30*/
/** \internal
 *	Class reads in ECMWF data from NETCDF4 files.
 **/
/*---------------------------------------------------------------------------*/

class sknetcdf_IONetCDF
{
	private:
		std::set<sknetcdfVariable>							m_fields;
		typedef std::set<sknetcdfVariable>::const_iterator	fielditerator;

	private:
		double						m_missingvalue;

	private:
		bool						AddEmptyField					( const char* name);
		bool						LoadFields						( const char* filename);

	public:
									sknetcdf_IONetCDF					();
		virtual					   ~sknetcdf_IONetCDF					();
		const sknetcdfVariable*		Field								( const char* name );
		bool						LoadFile							( const char* fullfilename, const char* fieldlist, double missingvalue);
		bool						VariableProfile_InterpolateLatLong	( const char* varname, double latitude, double longitude, std::vector<double>* profile);
		bool						IsLoaded							() const { return (m_fields.size() > 0);}
};

