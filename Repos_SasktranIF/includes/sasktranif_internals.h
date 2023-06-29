#pragma once

#include "sasktranif.h"
#include "sasktran_interfaces.h"

/*-----------------------------------------------------------------------------
 *					SasktranIF_ClassFactoryLocator		2014-2-8*/
 /** **/
 /*---------------------------------------------------------------------------*/

class SasktranIF_ClassFactoryLocator
{
private:
	bool			FindRegistrySetting	(const char*classname, const char* entityname, std::string* userdllname);
	void			AssignDLLname		(char** userdllname, const std::string& dllname);
	//		bool			LoadFunctionFromDLL( const char* dllname, const char* functionname, void** funcptr );
	//		bool			InitializeDLLLogger( void* funcptr );

public:
	bool			CreateISKEngine				(const char* enginename,      ISKEngine_Stub**           engine,      char** dllname);
	bool			CreateISKClimatology		(const char* climatologyname, ISKClimatology_Stub**      climatology, char** dllname);
	bool			CreateISKOpticalProperty    (const char* optpropname,     ISKOpticalProperty_Stub**  optprop,     char** dllname);
	bool			CreateISKGeodetic			(const char* geoidname,       ISKGeodetic_Stub**		 geoid,       char** dllname);
	bool			CreateISKSolarSpectrum		(const char* solarname,       ISKSolarSpectrum_Stub**	 solar,       char** dllname);
	bool			CreateISKEmission			(const char* emissionname,    ISKEmission_Stub**		 emission,    char** dllname);
	bool			CreateISKBrdf				(const char* emissionname,    ISKBrdf_Stub**			 brdf,        char** dllname);
	bool			CreateISKMie				(const char* miename,		  ISKMie_Stub**				 mie,		  char** dllname);
};

