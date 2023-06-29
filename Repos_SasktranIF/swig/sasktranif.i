/* File : sasktranif.i */
%module sasktranif

%{
#define SWIG_FILE_WITH_INIT
#include "sasktranif_internals.h"
#include "nxbase_math.h"

static PyObject *g_sasktranifError;
%}

%include "typemaps.i"
%include "numpy.i"

%init %{
import_array();
g_sasktranifError = PyErr_NewException("sasktranif.functionfail", NULL, NULL);
Py_INCREF(g_sasktranifError);
PyModule_AddObject(m, "functionfail", g_sasktranifError);
%}


%exception{
	try {
	$action
	}
	catch (const std::exception &exc)
	{
	    PyErr_SetString( PyExc_RuntimeError, exc.what() );		// catch anything thrown within try block that derives from std::exception
	}
}


//%{
/*-----------------------------------------------------------------------------
 *					CopyRadianceScalarArray		 2015- 9- 29*/
/** Copies a contiguous 2-D array from C++ to a standard python 2-D array
 *	represented as a list of lists. This could be upgraded to a 2-D numpy array
 *	in the future.
 **/
/*---------------------------------------------------------------------------*/

//static PyObject* CopyRadianceScalarArray( const double* radiance, int numlinesofsight, int numwavelens)
//{
//	PyObject*				full_object;
//	PyObject*				onewavelenobject;
//	const double*					radianceptr = radiance;
//
//	full_object = PyTuple_New( numwavelens );
//	for (Py_ssize_t iw = 0; iw < (Py_ssize_t)numwavelens; iw++)
//	{
//		onewavelenobject = PyTuple_New( numlinesofsight );
//		for (Py_ssize_t pos = 0; pos < (Py_ssize_t)numlinesofsight; pos++)
//		{
//			PyTuple_SetItem(onewavelenobject, pos, PyFloat_FromDouble( *radianceptr) );
//			radianceptr++;
//		}
//		PyTuple_SetItem( full_object, iw, onewavelenobject );
//	}
//	return full_object;
//}
//%}


//%{
//PyTypeObject*	m_geoid_typeobject= nullptr;
//
//PyObject* Create_GeodeticInstantObject(const GEODETIC_INSTANT& geoid) 
//{
//	PyObject*				geoid_object;
//
//	if (m_geoid_typeobject == nullptr)
//	{
//		PyStructSequence_Desc	desc;
//		PyStructSequence_Field	fields[5];
//		fields  [0].name = "latitude";
//		fields  [0].doc  = "The geodetic latitude in degrees";
//		fields  [1].name = "longitude";
//		fields  [1].doc  = "The geodetic longitude in degrees east.";
//		fields  [2].name = "heightm";
//		fields  [2].doc  = "The height in meters above sea level";
//		fields  [3].name = "mjd";
//		fields  [3].doc  = "The modified julian date in days" ;
//		fields  [4].name = nullptr;
//		fields  [4].doc  = nullptr;
//
//		desc.name   = "GEODETIC_INSTANT";
//		desc.doc    = "The 4-D coordinates of a point in space and time.";
//		desc.fields = fields;
//		desc.n_in_sequence = 4;
//		m_geoid_typeobject = PyStructSequence_NewType( &desc );
//	}
//	geoid_object = PyStructSequence_New( m_geoid_typeobject );
//	PyStructSequence_SetItem(geoid_object, 0, PyFloat_FromDouble( geoid.latitude) );
//	PyStructSequence_SetItem(geoid_object, 1, PyFloat_FromDouble( geoid.longitude) );
//	PyStructSequence_SetItem(geoid_object, 2, PyFloat_FromDouble( geoid.heightm) );
//	PyStructSequence_SetItem(geoid_object, 3, PyFloat_FromDouble( geoid.mjd) );
//	return geoid_object;
//}
//
//%}

/*-----------------------------------------------------------------------------
 *					Typemap for ISKEngine::CalculateRadiance		 2015- 10- 1*/
/**  Return the C++ array as a 2-D array
 **/
/*---------------------------------------------------------------------------*/

%typemap(in, numinputs=0) ( const double** radiance, int* numwavelens, int* numlinesofsight)  
						   (double *data_temp,  int numwave_temp, int numlos_temp)
{
  $1 = &data_temp;					// Let radiance point to local variable, data_temp
  $2 = &numwave_temp;				// Let numwavelengths point to local variable
  $3 = &numlos_temp;				// Let numlinesofsight point to local variable
}

%typemap(argout) ( const double** radiance, int* numwavelens, int* numlinesofsight)	
{
	PyObject*				radarray;									// Object to hold the numpy object copied form the radiance
	double*					outptr;										// Pointer to the data int the numpy object
	const double*			dataptr = *($1);
	int			            nw      = *($2);
	int			            nlos    = *($3);
	npy_intp				dims[2] = { nw, nlos };						// The dimensions of the numpy object
	int						numelem = nw*nlos;

	if (numelem <= 0)
	{
		radarray = Py_None;
		Py_INCREF(radarray);
	}
	else
	{
		radarray = PyArray_SimpleNew(2, dims, NPY_DOUBLE);				// Create the numpy 2-D array
		if (!radarray) SWIG_fail;										// see if it failed
		outptr = (double *)PyArray_DATA((PyArrayObject*)radarray);		// Get a pointer to the start of the array
		for (int wavidx = 0; wavidx < nw; wavidx++)						// c++ radiance array has lines of sight as the leading dimension
		{																// but numpy array has wavelength as the leading dimension
			for (int losidx = 0; losidx < nlos; losidx++)
			{
				outptr[wavidx + losidx *nw] = dataptr[wavidx + losidx *nw];
			}
		}
	}
	$result = SWIG_Python_AppendOutput($result,radarray);
}

/*-----------------------------------------------------------------------------
 *					Typemap for ISKEngine::CalculateRadiancePolarized		 2015- 10- 1*/
/**  Return the C++ array as a 2-D array of IKStokesVector objects
 **/
/*---------------------------------------------------------------------------*/

%typemap(in, numinputs=0) ( const ISKStokesVector** radiancep, int* numwavelens, int* numlinesofsight)  
						   (ISKStokesVector* data_temp,  int numwave_temp, int numlos_temp)
{
  $1 = &data_temp;				// Let radiance point to local variable, data_temp
  $2 = &numwave_temp;				// Let numwavelengths point to local variable
  $3 = &numlos_temp;				// Let numlinesofsight point to local variable
}

%typemap(argout) ( const ISKStokesVector** radiancep, int* numwavelens, int* numlinesofsight)	
{
	PyObject*				radarray;									// Object to hold the numpy object copied form the radiance
	PyObject**				outptr;										// Pointer to the data int the numpy object
	PyObject*				oldobject;									// Pointer to the data int the numpy object
	PyObject*				obj;										// Pointer to the data int the numpy object
	const ISKStokesVector*	dataptr = *($1);
	int			            nw      = *($2);
	int			            nlos    = *($3);
	npy_intp				dims[2] = { nw, nlos };						// The dimensions of the numpy object
	int						numelem = nw*nlos;

	if (numelem <= 0)
	{
		radarray = Py_None;
		Py_INCREF(radarray);
	}
	else
	{
		radarray = PyArray_SimpleNew(2, dims, NPY_OBJECT);					// Create the numpy 2-D array
		if (!radarray) SWIG_fail;											// see if it failed
		outptr = (PyObject **)PyArray_DATA((PyArrayObject*)radarray);		// Get a pointer to the start of the array
		for (int wavidx = 0; wavidx < nw; wavidx++)
		{
			for (int losidx = 0; losidx < nlos; losidx++)
			{
				const ISKStokesVector&	iquv = dataptr[wavidx + losidx*nw];
				obj = SWIG_NewPointerObj( new ISKStokesVector(iquv), SWIGTYPE_p_ISKStokesVector, SWIG_POINTER_OWN |  0 );
				oldobject = outptr[wavidx + losidx*nw];
				outptr[wavidx + losidx*nw] = obj;
				Py_XDECREF(oldobject);
			}
		}
	}
	$result = SWIG_Python_AppendOutput($result,radarray);
}

/*-----------------------------------------------------------------------------
 *					Typemap for ISKEngine::GetWeightingFunctions		 2015- 10- 1*/
/** 
**/
/*---------------------------------------------------------------------------*/

%typemap(in, numinputs=0) ( const double** wf, int* numwavelens, int* numlinesofsight, int* numwf)  
						   (double * data_temp,  int dim1_temp, int dim2_temp, int dim3_temp)
{
  $1 = &data_temp;				// Let wf point to local variable, data_temp
  $2 = &dim1_temp;				// Let numwavel point to local variable, dim1_temp
  $3 = &dim2_temp;				// Let numlinesofsight point to local variable, dim2_temp
  $4 = &dim3_temp;				// Let numwf point to local variable, dim3_temp
}

%typemap(argout) (const double** wf, int* numwavelens, int* numlinesofsight, int* numwf)
{
	PyObject*	wfarray;
	double*		outptr;
	npy_intp	dims[3] = { *($2), *($3), *($4) };
	int			numelem = (*$2)*(*$3)*(*$4);

	wfarray = PyArray_SimpleNew(3, dims, NPY_DOUBLE);
	if (!wfarray) SWIG_fail;
	outptr = (double *)PyArray_DATA((PyArrayObject*)wfarray);
	for (int i = 0; i < numelem; i++) outptr[i] = (*$1)[i];
	$result = SWIG_Python_AppendOutput($result, wfarray);
}

/*-----------------------------------------------------------------------------
 *					Typemap for ISKEngine::GetPropertyArray		 2015- 10- 1*/
/**  Return the C++ array as a 1-D array, 
 * In python:
 * ok,array = ISKEngine::GetPropertyArray( name)
 **/
/*---------------------------------------------------------------------------*/

%typemap(in, numinputs=0) (const double** value, int* numpoints)  
						  (double * data_temp,  int dim1_temp)
{
  $1 = &data_temp;				// Let radiance point to local variable, data_temp
  $2 = &dim1_temp;				// Let numlinesofsight point to local variable, dim1_temp
}

%typemap(argout) (const double** value, int* numpoints)	
{
	PyObject*	radarray;									// Object to hold the numpy object copied form the radiance
	double*		outptr;										// Pointer to the data int the numpy object
	npy_intp	dims[1] = { *($2)};					// The dimensions of the numpy object
	int			numelem = (*$2);

	radarray = PyArray_SimpleNew(1, dims, NPY_DOUBLE);				// Create the numpy 2-D array
	if (!radarray) SWIG_fail;										// see if it failed
	outptr = (double *)PyArray_DATA((PyArrayObject*)radarray);		// Get a pointer to the start of the array
	for (int i = 0; i < numelem; i++) outptr[i] = (*$1)[i];
	$result = SWIG_Python_AppendOutput($result,radarray);
}


/*-----------------------------------------------------------------------------
 *				Input Typemap for ISKOpticalProperty::CalculateCrossSections		 2015- 10- 1*/ 
/**  		
 *	Accepts a scalar or array
 *	Copy the input wavelength array into a locally allocated contiguous sequence and arrange for
 *	destruction in the [freearg] section. The absorption, scattering and extinction cross-sections
 *  are created locally and passed beck to the python code. 
 *	This typemap maps the C++ function to python function
 *	(ok, absxs, extxs, scattxs) = ISKOpticalProperty::CalculateCrossSectionsArray( wavenumber)
 *	$1 = wavenumber
 *	$2 = absxs
 *	$3 = extxs
 *  $4 = scattxs
 *	$5 = numortype
 **/
/*---------------------------------------------------------------------------*/
%typemap(in, numinputs=1) ( const double * wavenumber,  double *absxs, double* extxs, double* scattxs, int numortype)  
						   (double         scalarwavenumber, double         scalarabsxs, double         scalarextxs, double         scalarscatxs,
						    PyArrayObject* objectwavenumber, PyArrayObject* objectabsxs, PyArrayObject* objectextxs, PyArrayObject* objectscatxs, int isnewarray)
{
	isnewarray   = 0;
	if (PySequence_Check($input) || PyArray_Check($input))												// See if this is a Python sequence or numpy array
	{
		int numpoints;																											// If it is then
		npy_intp dims[1];

		objectwavenumber =  obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE, &isnewarray);				// Make a new contiguous array object for Python
		numpoints        =  (int)PyArray_Size( (PyObject*)objectwavenumber);
		dims[0]          = numpoints;
		objectabsxs      = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_DOUBLE);				// Create the absorption numpy array
		objectextxs      = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_DOUBLE);				// Create the extinction numpy array
		objectscatxs     = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_DOUBLE);				// Create the scattering numpy array
		if ((!objectabsxs) || !(objectextxs) || !(objectscatxs)) SWIG_fail;		// see if it failed

		$1 = (double *)PyArray_DATA( objectwavenumber);																	// Get a pointer to the start of the absorption array
		$2 = (double *)PyArray_DATA( objectabsxs);
		$3 = (double *)PyArray_DATA( objectextxs);
		$4 = (double *)PyArray_DATA( objectscatxs);
		$5 =  numpoints;
	}
	else if ( PyNumber_Check( $input))																				// Is this a scalar number
	{
		scalarwavenumber = PyFloat_AsDouble($input);
		$1 = &scalarwavenumber;
		$2 = &scalarabsxs;
		$3 = &scalarextxs;
		$4 = &scalarscatxs;
		$5 = -1;
	}
	else
	{
		$1 = nullptr;
		$2 = nullptr;
		$3 = nullptr;
		$4 = nullptr;
		$5 =  -9999;
		SWIG_Python_SetErrorMsg(PyExc_ValueError, "Expected a scalar or array "); 
		SWIG_fail;
	}
}

%typemap(argout) ( const double * wavenumber,  double *absxs, double* extxs, double* scattxs, int numortype)  
{
	int n = $5;
	if (n >= 0)
	{
		$result = SWIG_Python_AppendOutput($result, (PyObject*)objectabsxs$argnum);
		$result = SWIG_Python_AppendOutput($result, (PyObject*)objectextxs$argnum);
		$result = SWIG_Python_AppendOutput($result, (PyObject*)objectscatxs$argnum);
	}
	else
	{
		$result = SWIG_Python_AppendOutput($result, PyFloat_FromDouble(scalarabsxs$argnum));
		$result = SWIG_Python_AppendOutput($result, PyFloat_FromDouble(scalarextxs$argnum));
		$result = SWIG_Python_AppendOutput($result, PyFloat_FromDouble(scalarscatxs$argnum));

	}
}

%typemap(freearg)  (const double * wavenumber,  double *absxs, double* extxs, double* scattxs, int numortype)
{
  if (isnewarray$argnum && objectwavenumber$argnum)
    { Py_DECREF(objectwavenumber$argnum); }
}
/*-----------------------------------------------------------------------------
 *				Input Typemap for double* IN_ARRAY1, double*  OUTARRAY1, int NUMPOINTS		 2015- 10- 1*/ 
/**  		
 *	Accepts a scalar or array
 *	Copy the input (wavelength) array into a locally allocated contiguous sequence and arrange for
 *	destruction in the [freearg] section. The output array is created locally and passed back to the python code. 
 *	This typemap maps the C++ function to python function
 *	(ok, radiance) = ISKEmission::IsotropicEmiision( wavenumbers )
 *
 *	$1 = wavenumber
 *	$2 = isotrpoicradiance
 *	$3 = numortype
 **/
/*---------------------------------------------------------------------------*/
%typemap(in, numinputs=1) ( double* IN_ARRAY1, double*  OUTARRAY1, int NUMPOINTS) 
						   (double         scalarwavenumber, double         scalarrad, 
						    PyArrayObject* objectwavenumber, PyArrayObject* objectrad,int isnewarray)
{
	isnewarray   = 0;
	if (PySequence_Check($input) || PyArray_Check($input))												// See if this is a Python sequence or numpy array
	{
		int numpoints;																											// If it is then
		npy_intp dims[1];

		objectwavenumber =  obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE, &isnewarray);				// Make a new contiguous array object for Python
		numpoints        =  (int)PyArray_Size( (PyObject*)objectwavenumber);
		dims[0]          = numpoints;
		objectrad        = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_DOUBLE);				// Create the absorption numpy array
		if ((!objectrad)) SWIG_fail;		// see if it failed

		$1 = (double *)PyArray_DATA( objectwavenumber);																	// Get a pointer to the start of the absorption array
		$2 = (double *)PyArray_DATA( objectrad);
		$3 =  numpoints;
	}
	else if ( PyNumber_Check( $input))																				// Is this a scalar number
	{
		scalarwavenumber = PyFloat_AsDouble($input);
		$1 = &scalarwavenumber;
		$2 = &scalarrad;
		$3 = -1;
	}
	else
	{
		$1 = nullptr;
		$2 = nullptr;
		$3 =  -9999;
		SWIG_Python_SetErrorMsg(PyExc_ValueError, "Expected a scalar or array "); 
		SWIG_fail;
	}
}

%typemap(argout) ( double* IN_ARRAY1, double*  OUTARRAY1, int NUMPOINTS)  
{
	int n = $3;

	if (n >= 0)
	{
		$result = SWIG_Python_AppendOutput($result, (PyObject*)objectrad$argnum);
	}
	else
	{
		$result = SWIG_Python_AppendOutput($result, PyFloat_FromDouble(scalarrad$argnum));

	}
}

%typemap(freearg)  (double* IN_ARRAY1, double*  OUTARRAY1, int NUMPOINTS)
{
  if (isnewarray$argnum && objectwavenumber$argnum)
    { Py_DECREF(objectwavenumber$argnum); }
}
/*-----------------------------------------------------------------------------
 *				Input Typemap for const CLIMATOLOGY_HANDLE&			2015- 10- 1*/ 
/*---------------------------------------------------------------------------*/

%typemap(in) ( const CLIMATOLOGY_HANDLE& species)
{
	void*  objectptr;
	int	   result;
	nxString  name;	

	if (PyBytes_Check($input))
	{
		name = PyBytes_AsString($input);
		$1   = FindGlobalClimatologyHandle(name);
	}
	else if (PyUnicode_Check($input))
	{
		nxStringw		w;
		int			type;
		int				xlen = sizeof(wchar_t);
      
		PyUnicode_READY($input);
		type = PyUnicode_KIND($input);
	
		if      ( type == PyUnicode_1BYTE_KIND) name = (char *)PyUnicode_1BYTE_DATA($input);
		else if ( (type == PyUnicode_2BYTE_KIND) && (xlen ==2) ) 
		{
			w    = (wchar_t *)PyUnicode_2BYTE_DATA($input);
			name = w.ConvertToChar();
		}
		else if ( type == PyUnicode_4BYTE_KIND && (xlen ==2))
		{
			w    = (wchar_t *)PyUnicode_4BYTE_DATA($input);
			name = w.ConvertToChar();
		}
		else
		{
			SWIG_Python_SetErrorMsg(PyExc_ValueError, "Error converting string to simple ascii string. Current sasktranIF only handles simple strings. Error fetching CLIMATOLOGY_HANDLE"); 
			SWIG_fail;
		}
		name.RemoveWhiteSpace();
		$1   = FindGlobalClimatologyHandle(name);
	}
	else
	{

		result = SWIG_ConvertPtr($input, &objectptr, SWIGTYPE_p_GUID, 0 |  0 );
		if (!SWIG_IsOK(result))
		{
			SWIG_Python_SetErrorMsg(PyExc_ValueError, " error converting argument to CLIMATOLOGY_HANDLE const &"); 
			SWIG_fail;
		}
		$1 = reinterpret_cast< CLIMATOLOGY_HANDLE * >(objectptr);
	}
}

/*-----------------------------------------------------------------------------
 *	Input Typemap for a CLIMATOLOGY_HANDLE& 2019-09-18 */
/* 
 *	Typemap for.
 *		ISKCLimatology::GetParameter          ( const char * climatology_handle_name,  const GEODETIC_INSTANT& location, double* valueout );
 *		ISKCLimatology::GetHeightProfile      ( const char * climatology_handle_name,  GEODETIC_INSTANT location, const double* altitude, double *profile, int numalts );
 *		ISKCLimatology::SetPropertyUserDefined( const char * climatology_handle_name,  double* profilevalues, int numpoints);
 *		ISKEngine::AddSpecies                 ( const char*  climatology_handle_name,  ISKClimatology& climatology, ISKOpticalProperty& opticalproperty);
 *		ISKEngine::AddEmission		      ( const char*  climatology_handle_name,  ISKEmission&    emission );
 *
 *	It detects if the python has passed in a string or a GUID.  String is the proper format for the future
 *	but if a user has passed in a GUID then  try and convert it it to a string. There is a penalty hit in this
 *	as the conversion from handle to string  is not particularly efficient.
 *
 *	I have left the handle conversion in the code so we are backwardly compatible although I would like
 *	everyone to convert over, hence I have put a warning message in there.
*/
/*---------------------------------------------------------------------------*/

%typemap(in) ( const char* HANDLENAME )
{
	void*  		objectptr;
	int	   	result;

	if (PyBytes_Check($input))
	{
		$1 = (char*)(intptr_t)PyBytes_AsString($input);
		
	}
	else if (PyUnicode_Check($input))
	{
		$1   = (char*)(intptr_t)PyUnicode_AsUTF8($input);
	}
	else
	{
		static bool firsttime = true;
		if (firsttime)
		{
			nxLog::Record(NXLOG_WARNING,"Deprecated feature. We are dropping support for passing GUID handles through the sasktran interface\n."
			                            "Please replace the GUID with a corresponding string \n"
			                            "e.g. replace skif.SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3 with 'SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3'");
			firsttime = false;
		}
		result = SWIG_ConvertPtr($input, &objectptr, SWIGTYPE_p_GUID, 0 |  0 );
		if (!SWIG_IsOK(result))
		{
			SWIG_Python_SetErrorMsg(PyExc_ValueError, " error converting argument to CLIMATOLOGY_HANDLE const &"); 
			SWIG_fail;
		}
		$1 = (char*)(intptr_t)FindGlobalClimatologyNameOfHandle( (*(CLIMATOLOGY_HANDLE*)(objectptr)) );
	}
}

/*-----------------------------------------------------------------------------
 *					Typemap for nxVector conversions				2015- 10- 1*/
/*---------------------------------------------------------------------------*/

%typemap(out) const nxVector& 
{
	npy_intp dims[1] = { 3 };
	double*	 v;

   $result = PyArray_SimpleNew(1, dims, NPY_DOUBLE);				// Create the absorption numpy 2-D array
	v = (double *)PyArray_DATA((PyArrayObject*)$result);			// Get a pointer to the start of the absorption array
	v[0] = ($1)->X();
	v[1] = ($1)->Y();
	v[2] = ($1)->Z();

}

%typemap(in, numinputs=0) (nxVector*)
						  (nxVector tnxv)
{
  $1 = &tnxv;
}

%typemap(argout) nxVector*
{
	npy_intp dims[1] = { 3 };
	double* v$argnum;
	PyObject* tresult$argnum = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
	v$argnum = (double *)PyArray_DATA((PyArrayObject*)tresult$argnum);
	v$argnum[0] = (*$1).X();
	v$argnum[1] = (*$1).Y();
	v$argnum[2] = (*$1).Z();
	$result = SWIG_Python_AppendOutput($result, tresult$argnum);
}

%typemap(in) (const nxVector&)
			 (nxVector temp)
{
	if (PySequence_Check($input))																// See if this is a Python sequence
	{																							// If it is then
		if (PySequence_Length($input) != 3)												// Make sure its the correct size
		{
			SWIG_Python_SetErrorMsg(PyExc_ValueError,"Size mismatch. Expected 3 elements");
			SWIG_fail;
		}
		bool ok = true;
		PyObject *o1 = PySequence_GetItem($input,0);
		PyObject *o2 = PySequence_GetItem($input,1);
		PyObject *o3 = PySequence_GetItem($input,2);
		ok = ok && PyNumber_Check(o1) && PyNumber_Check(o2) && PyNumber_Check(o3);
		if( ok )
		{
			temp.SetCoords(PyFloat_AsDouble(o1), PyFloat_AsDouble(o2), PyFloat_AsDouble(o3));
		}
		else 
		{
			SWIG_Python_SetErrorMsg(PyExc_ValueError,"Sequence elements must be numbers");      
			SWIG_fail;
		}
		$1 = &temp;
	}
	else if (PyArray_Check($input))
	{
		PyArrayObject*	obj = (PyArrayObject* )($input);
		double*			d;

		if ( (PyArray_Size($input) != 3) || (PyArray_TYPE(obj) != NPY_DOUBLE))												// Make sure its the correct size
		{
			SWIG_Python_SetErrorMsg(PyExc_ValueError,"Size mismatch. Expected 3 elements of doubles");
			SWIG_fail;
		}
		d = (double*) PyArray_DATA(obj);
		temp.SetCoords(d[0], d[1], d[2]);
		$1 = &temp;
	}
	else
	{
		SWIG_Python_SetErrorMsg(PyExc_ValueError,"Expected Input as a list, numpy array");
		SWIG_fail;
	}
}

%typemap(out) nxVector
{
	npy_intp dims[1] = { 3 };
	double*	 v;

   $result = PyArray_SimpleNew(1, dims, NPY_DOUBLE);				// Create the absorption numpy 2-D array
	v = (double *)PyArray_DATA((PyArrayObject*)$result);			// Get a pointer to the start of the absorption array
	v[0] = ($1).X();
	v[1] = ($1).Y();
	v[2] = ($1).Z();
}

/*-----------------------------------------------------------------------------
 *					Typemap for GEODETIC_INSTANT conversions				2015- 10- 1*/
/*---------------------------------------------------------------------------*/

%typemap(in, numinputs=1) (const GEODETIC_INSTANT&)
						  (GEODETIC_INSTANT temp)
{
	if (PySequence_Check($input))																// See if this is a Python sequence
	{																							// If it is then
		if (PySequence_Length($input) != 4)												// Make sure its the correct size
		{
			SWIG_Python_SetErrorMsg(PyExc_ValueError,"Size mismatch. Expected 4 elements");
			SWIG_fail;
		}
		bool ok = true;
		PyObject *o1 = PySequence_GetItem($input,0);
		PyObject *o2 = PySequence_GetItem($input,1);
		PyObject *o3 = PySequence_GetItem($input,2);
		PyObject *o4 = PySequence_GetItem($input,3);
		ok = ok && PyNumber_Check(o1) && PyNumber_Check(o2) && PyNumber_Check(o3) && PyNumber_Check(o4);
		if( ok )
		{
			temp.latitude = PyFloat_AsDouble(o1);
			temp.longitude = PyFloat_AsDouble(o2);
			temp.heightm = PyFloat_AsDouble(o3);
			temp.mjd = PyFloat_AsDouble(o4);
		}
		else 
		{
			SWIG_Python_SetErrorMsg(PyExc_ValueError,"Sequence elements must be numbers");      
			SWIG_fail;
			return NULL;
		}
		$1 = &temp;
	}
	else if (PyArray_Check($input))
	{
		PyArrayObject*	obj = (PyArrayObject* )($input);
		double*			d;

		if ( (PyArray_Size($input) != 4) || (PyArray_TYPE(obj) != NPY_DOUBLE))												// Make sure its the correct size
		{
			SWIG_Python_SetErrorMsg(PyExc_ValueError,"Size mismatch. Expected 3 elements of doubles");
			SWIG_fail;
		}
		d = (double*) PyArray_DATA(obj);
		temp.latitude = d[0];
		temp.longitude = d[1];
		temp.heightm = d[2];
		temp.mjd = d[4];
		$1 = &temp;
	}
	else
	{
		SWIG_Python_SetErrorMsg(PyExc_ValueError,"Expected Input as a list, numpy array");
		SWIG_fail;
	}
}

%typemap(out) GEODETIC_INSTANT
{
	npy_intp dims[1] = { 4 };
	double*	 v;

   $result = PyArray_SimpleNew(1, dims, NPY_DOUBLE);				// Create the absorption numpy 2-D array
	v = (double *)PyArray_DATA((PyArrayObject*)$result);			// Get a pointer to the start of the absorption array
	v[0] = ($1).latitude;
	v[1] = ($1).longitude;
	v[2] = ($1).heightm;
	v[3] = ($1).mjd;
}

%typemap(out) const IQUV&
{
	npy_intp dims[1] = { 4 };
	double*	 v;

   $result = PyArray_SimpleNew(1, dims, NPY_DOUBLE);				// Create the absorption numpy 2-D array
	v = (double *)PyArray_DATA((PyArrayObject*)$result);			// Get a pointer to the start of the absorption array
	v[0] = ($1)->I;
	v[1] = ($1)->Q;
	v[2] = ($1)->U;
	v[3] = ($1)->V;
}

/*-----------------------------------------------------------------------------
 *				Typemap for nxVector::FromSequence( const double fixedarray[3])	2015- 10- 1 
 *				and GEODETIC_INSTANT::FromSequence( const double fixedarray[4]) */
/**  		
 **/
/*---------------------------------------------------------------------------*/

%typemap(in) double fixedarray[ANY] (double temp[$1_dim0]) 
{
	int i;

	if (PySequence_Check($input))																// See if this is a Python sequence
	{																							// If it is then
		if (PySequence_Length($input) != $1_dim0)												// Make sure its the correct size
		{
			SWIG_Python_SetErrorMsg(PyExc_ValueError,"Size mismatch. Expected $1_dim0 elements");
			SWIG_fail;
		}
		for (i = 0; i < $1_dim0; i++) 
		{
			PyObject *o = PySequence_GetItem($input,i);
			if (PyNumber_Check(o)) 
			{
				temp[i] = PyFloat_AsDouble(o);
			} 
			else 
			{
				SWIG_Python_SetErrorMsg(PyExc_ValueError,"Sequence elements must be numbers");      
				SWIG_fail;
			}
		}
	}
	else if (PyArray_Check($input))
	{
		PyArrayObject*	obj = (PyArrayObject* )($input);
		double*			d;

		if ( (PyArray_Size($input) != $1_dim0) || (PyArray_TYPE(obj) != NPY_DOUBLE))												// Make sure its the correct size
		{
			SWIG_Python_SetErrorMsg(PyExc_ValueError,"Size mismatch. Expected $1_dim0 elements of doubles");
			SWIG_fail;
		}
		d = (double*) PyArray_DATA(obj);
		for (i = 0; i < $1_dim0; i++) 
		{
			temp[i] = d[i];
		}
	}
	else
	{
		SWIG_Python_SetErrorMsg(PyExc_ValueError,"Expected a sequence");
		SWIG_fail;
	}

  $1 = temp;
}

/*-----------------------------------------------------------------------------
 *				Typemap for ISKModuleBase::SetProperty			2015- 10- 1*/ 
/**  		
 *	Accepst a scalar, array or sequence or an ISKModuleBase object 
 **/
/*---------------------------------------------------------------------------*/
  
%typemap(in, numinputs=1) ( void* valueorobject, int numpoints_or_type)  
						   (void* objectptr,  void* arrayptr, double scalarvalue, PyArrayObject* newarray, int isnewarray)
{
	isnewarray   = 0;
	newarray = nullptr;
	if (PyUnicode_Check($input))																				// See if we have a string coming in
	{
		$1 = (void *)(intptr_t)PyUnicode_AsUTF8($input);
		$2 = -3;
	}
	else if (PySequence_Check($input) || PyArray_Check($input))													// See if this is a Python sequence or numpy array
	{																											// If it is then
		newarray =  obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE, &isnewarray);					// Make a new contiguous array object for Python
		arrayptr = 	PyArray_DATA( newarray);																	// Get a pointer to the start of the absorption array
		$1       =  arrayptr;
		$2       =  (int)PyArray_Size( (PyObject*)newarray);
	}
	else if ( PyNumber_Check( $input))																			// Is this a scalar number
	{
		scalarvalue = PyFloat_AsDouble($input);
		$1 = &scalarvalue;
		$2 = -1;
	}

	else if ( SWIG_IsOK( SWIG_ConvertPtr($input, &objectptr, SWIGTYPE_p_ISKModuleBase, 0 |  0 ) ) )				// IS this an instance of ISKModuleBase or a derived class
	{
		$1 = objectptr;
		$2 = -2;
	}
	else
	{
		$1   = nullptr;
		$2   = -9999;
		SWIG_Python_SetErrorMsg(PyExc_ValueError, "Expected a scalar, array or SasktranIF object derived from ISKModuleBase"); 
		SWIG_fail;
	}
}


%typemap(freearg)  ( void* valueorobject, int numpoints_or_type)
{
  if (isnewarray$argnum && newarray$argnum)
    { Py_DECREF(newarray$argnum); }
}


%typemap(out) (bool)
{
    if (!$1)
	{
        PyErr_SetString(g_sasktranifError, "Sasktran Interface Function returned NOT OKAY status");
        return NULL;
    }
	$result = PyBool_FromLong( ($1) ? 1 : 0);
}


/*-----------------------------------------------------------------------------
 *				Typemap for ISKModuleBase::GetProperty			2015- 10- 1*/ 
/**  		
 *	Accepst a scalar or double array  from ISKModuleBase::GetProperty
 **/
/*---------------------------------------------------------------------------*/
  
%typemap(in, numinputs=0) ( const double** propertyvalue, int* numpoints)  
						  ( double* propptr,  int numpts)
{
	$1 = (double**)&propptr;
	$2 = &numpts;
}

%typemap(argout) ( const double** propertyvalue, int* numpoints)	
{
	PyObject*				radarray;									// Object to hold the numpy object copied form the radiance
	double*					outptr;										// Pointer to the data int the numpy object
	const double*			dataptr = *($1);
	int			            npts    = *($2);
	npy_intp				dims[1] = { npts };						// The dimensions of the numpy object

	if ( dataptr == NULL)
	{
		radarray = Py_None;
		Py_INCREF(radarray);
	}
	else
	{
		if (npts == 0)															// IF we are returning a scalar value
		{
			radarray = PyFloat_FromDouble(*dataptr);							// Create the scalar value numpy 2-D array
			if (!radarray) SWIG_fail;											// see if it failed
		}
		else
		{
			radarray = PyArray_SimpleNew(1, dims, NPY_DOUBLE);					// Create the numpy array
			if (!radarray) SWIG_fail;											// see if it failed
			outptr = (double *)PyArray_DATA((PyArrayObject*)radarray);			// Get a pointer to the start of the array
			for (int i = 0; i < npts;  i++)
			{
				outptr[i] = dataptr[i];
			}
		}
	}
	$result = SWIG_Python_AppendOutput($result,radarray);
}


%apply const char* HANDLENAME {( const char* climatology_handle_name )};


%apply double *OUTPUT	{ (double* valueout),	        // ISKClimatology::GetParameter
			  (double* minwavelength_nm),	// ISKSolarSpectrum::MinValidWavelength
			  (double* maxwavelength_nm),	// ISKSolarSpectrum::MaxValidWavelength
			  (double* resolution_nm_fwhm),	// ISKSolarSpectrum::NanometerResolutionFWHM
			  (double* sample_spacing),	// ISKSolarSpectrum::SampleSpacing
              (double* return_brdf),	// ISKBrdf::BRDF
			  (double* real),			// ISKMie:: multiple
			  (double* imag)            // ISKMie:: multiple
			};

%apply      ( double ARGOUT_ARRAY1[ANY] )                           {   (double phasematrix[16])             // ISKOpticalProperty::CalculatePhaseMatrix
																	}


%apply int *OUTPUT		{ (int* losindex) };

%apply		(double* IN_ARRAY1, double*  OUTARRAY1, int NUMPOINTS ) {	( const double* wavenumber,        double* isotropicradiance,  int numorscalar),	// ISKEmission::IsotropicEmission
																		( const double* wavelen_nm_vacuum, double* irradiance,         int numpoints),		// ISKSolarSpectrum::Irradiance and Irradianceat1AU
																		( const double* wavelen_nm_vacuum, double* resolution_nm_fwhm, int numpoints),		// ISKSolarSpectrum::NanometerResolutionFWHMArray
																		( const double* wavelen_nm_vacuum, double* sample_spacing,     int numpoints ),		// ISKSolarSpectrum::SampleSpacingArray
																		( const double* altitude,          double *profile,            int numalts)			// ISKClimatology::GetHeightProfile
																	}

%apply		(double* IN_ARRAY1, int DIM1 ) 							{	(double* profilevalues,  int numpoints),		// ISKClimatology::SetPropertyUserDefined
																		(double* wavelen_nm,     int numwave),			// ISKOpticalProperty::AddUserDefined
																		(double* crosssection,   int numcross),			// ISKOpticalProperty::AddUserDefined
																		(double* pressure,       int numpressure),      // ISKOpticalProperty::AddUserDefinedPressure
																		(double* wavelen_nm,     int numwavel),         // ISKOpticalProperty::AddUserDefinedPressure
																		(double* temperature,    int numtemperature),   // ISKOpticalProperty::AddUserDefinedPressure
																		(const double* value,    int numpoints),		// SetPropertyArray
																		(const double* wavelen,  int numwavelen)		// ISKEngine::SetWavelengths
																	}

%typemap(in, numinputs=0) (std::complex<double>** s, int* numpoints)
(std::complex<double> *data_temp, int nump)
{
$1 = &data_temp;
$2 = &nump;
}

%typemap(argout) ( std::complex<double>** s, int* numpoints)
{
    PyObject*				out_real_array;
    PyObject*				out_imag_array;
    double*					out_real_ptr;
    double*                 out_imag_ptr;
    std::complex<double>*	dataptr = *($1);
    int			            np      = *($2);
    npy_intp				dims[1] = { np };						// The dimensions of the numpy object

    if (np <= 0)
    {
        out_real_array = Py_None;
        out_imag_array = Py_None;
        Py_INCREF(out_real_array);
        Py_INCREF(out_imag_array);
    }
    else
    {
        out_real_array = PyArray_SimpleNew(1, dims, NPY_DOUBLE);				// Create the numpy 2-D array
        out_imag_array = PyArray_SimpleNew(1, dims, NPY_DOUBLE);				// Create the numpy 2-D array
        if (!out_real_array) SWIG_fail;										// see if it failed
        if (!out_imag_array) SWIG_fail;										// see if it failed

        out_real_ptr = (double *)PyArray_DATA((PyArrayObject*)out_real_array);		// Get a pointer to the start of the array
        out_imag_ptr = (double *)PyArray_DATA((PyArrayObject*)out_imag_array);		// Get a pointer to the start of the array

        for (int idx = 0; idx < np; idx++)
        {
            out_real_ptr[idx] = dataptr[idx].real();
            out_imag_ptr[idx] = dataptr[idx].imag();
        }
    }

    $result = SWIG_Python_AppendOutput($result,out_real_array);
    $result = SWIG_Python_AppendOutput($result,out_imag_array);
}


%typemap(in, numinputs=0) (double** pmom, int* numlegendre)
(double *data_temp, int nump)
{
$1 = &data_temp;
$2 = &nump;
}

%typemap(argout) ( double** pmom, int* numlegendre)
{
    PyObject*				out_array;
    double*					out_ptr;
    double*	                dataptr = *($1);
    int			            np      = *($2);
    npy_intp				dims[2] = { 4, np };						// The dimensions of the numpy object

    if (np <= 0)
    {
        out_array = Py_None;
        Py_INCREF(out_array);
    }
    else
    {
        out_array = PyArray_SimpleNew(2, dims, NPY_DOUBLE);				// Create the numpy 2-D array
        if (!out_array) SWIG_fail;										// see if it failed

        out_ptr = (double *)PyArray_DATA((PyArrayObject*)out_array);		// Get a pointer to the start of the array

        for (int idx = 0; idx < np * 4; idx++)
        {
            out_ptr[idx] = dataptr[idx];
        }
    }

    $result = SWIG_Python_AppendOutput($result,out_array);
}



%include "../../Repos_BaseCode/nxbase/module/math/nxvector.h"
%include "../../Repos_BaseCode/nxbase/module/sktran_core/geodetic_instant.h"
%include "../../Repos_SasktranIF/includes/climatology_handles.h"
%include "../../Repos_SasktranIF/includes/sasktran_polarization.h"
%include "../../Repos_SasktranIF/includes/sasktran_interfaces.h"
