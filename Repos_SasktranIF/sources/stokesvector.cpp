#include "sasktranif_internals.h"


/*-----------------------------------------------------------------------------
 *					ISKBasisDirection::Assign		 2015- 10- 29*/
/** **/
/*---------------------------------------------------------------------------*/

void ISKBasisDirection::Assign( const nxVector& prop, const nxVector& theta, const nxVector& phi)
{
	m_propagation.SetCoords( prop.X(), prop.Y(), prop.Z());
	m_theta.SetCoords( theta.X(), theta.Y(), theta.Z());
	m_phi.SetCoords( phi.X(), phi.Y(), phi.Z());
}

/*-----------------------------------------------------------------------------
 *					ISKStokesVector::ISKStokesVector		 2015- 10- 28*/
/** **/
/*---------------------------------------------------------------------------*/

ISKStokesVector::ISKStokesVector()
{
	double nan = std::numeric_limits<double>::quiet_NaN();

	m_stokes.I = nan;
	m_stokes.Q = nan;
	m_stokes.U = nan;
	m_stokes.V = nan;
}

/*-----------------------------------------------------------------------------
 *					ISKStokesVector::ISKStokesVector		 2015- 9- 28*/
/** **/
/*---------------------------------------------------------------------------*/

ISKStokesVector::ISKStokesVector(  const IQUV& stokes , const ISKBasisDirection& new_basis )
{
	m_stokes = stokes;
	m_basis  = new_basis;
}


/*-----------------------------------------------------------------------------
 *					ISKStokesVector::Assign		 2015- 10- 29*/
/** **/
/*---------------------------------------------------------------------------*/

void ISKStokesVector::Assign(  const IQUV& stokes , const ISKBasisDirection& new_basis )
{
	m_stokes = stokes;
	m_basis  = new_basis;
}

/*-----------------------------------------------------------------------------
 *					to_new_basis		 2015- 9- 28*/
/** Converts the stokes vector to a new basis. A rotation matrix between the new basis and old basis is constructed
 *	and applied to the stokes vector.  Note that this process overrides the old basis
 *
 *	Parameters
 *	----------
 *	new_basis : The new basis.
 **/
/*---------------------------------------------------------------------------*/

void ISKStokesVector::to_new_basis( const ISKBasisDirection& new_basis)
{
	double		cos_eta;
	double		sin_eta;
	double		cos_two_eta;
	double		sin_two_eta;
	IQUV	S = m_stokes;

	bool ok = true;
	double cosPropDirAngle = m_basis.Propagation().Dot( new_basis.Propagation() );
	ok = ok && 0.999 < cosPropDirAngle;

	if(ok){
	    cos_eta =  (m_basis.Theta() & new_basis.Theta());
	    sin_eta = -(m_basis.Theta() & new_basis.Phi());
	    cos_two_eta = cos_eta*cos_eta - sin_eta*sin_eta;
	    sin_two_eta = 2.0*cos_eta*sin_eta;
	
		//rot = np.array([[ 1,           0,            0,  0],
					//	  [ 0, cos_two_eta, -sin_two_eta,  0],
		//                [ 0, sin_two_eta,  cos_two_eta,  0],
					//	  [ 0,           0,            0,  1]]);

		m_stokes.I = S.I;			//		 = rot.dot(m_stokes);
		m_stokes.Q = S.Q*cos_two_eta - S.U*sin_two_eta;
		m_stokes.U = S.Q*sin_two_eta + S.U*cos_two_eta;
		m_stokes.V = S.V;
		m_basis    = new_basis;
	} else{
		nxLog::Record( NXLOG_WARNING, "ISKStokesVector::to_new_basis, Propagation directions must be the same: dot(old,new)=%7.4e", cosPropDirAngle );
	}
}


/*-----------------------------------------------------------------------------
 *					ISKStokesVector::to_new_basis		 2015- 11- 25*/
/** **/
/*---------------------------------------------------------------------------*/

void ISKStokesVector::to_new_basis( const nxVector& prop, const nxVector& theta, const nxVector& phi)
{
    ISKBasisDirection   basis;

    basis.Assign( prop, theta, phi);
    to_new_basis(  basis );
}

#if 0
#include <Python.h>

/*-----------------------------------------------------------------------------
 *					class ISKPythonHelperFuncs		 2015- 9- 29*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKPythonHelperFuncs
{
	private:
		PyTypeObject*	m_istokes_typeobject;
		PyTypeObject*	m_iquv_typeobject;
		PyTypeObject*	m_basis_typeobject;

	public:
						ISKPythonHelperFuncs();
					   ~ISKPythonHelperFuncs			();
		PyObject*		Create_nxVectorObject			( const nxVector* v);
		PyObject*		Create_ISKBasisObject			(const ISKBasisDirection* basis);
		PyObject*		Create_ISK_IQUVObject			(const IQUV* radiance); 
		PyObject*		Create_ISKStokesVectorObject	(const ISKStokesVector* radiance);
		PyObject*		CopyRadianceArray				( ISKStokesVector* radiance, int numlinesofsight, int numwavelens);
};

/*-----------------------------------------------------------------------------
 *					ISKPythonHelperFuncs::ISKPythonHelperFuncs		 2015- 9- 29*/
/** **/
/*---------------------------------------------------------------------------*/

ISKPythonHelperFuncs::ISKPythonHelperFuncs()
{
	m_istokes_typeobject = nullptr;
	m_iquv_typeobject    = nullptr; 
	m_basis_typeobject     = nullptr;
}


/*-----------------------------------------------------------------------------
 *					ISKPythonHelperFuncs::~ISKPythonHelperFuncs		 2015- 9- 29*/
/** **/
/*---------------------------------------------------------------------------*/

ISKPythonHelperFuncs::~ISKPythonHelperFuncs()
{
}

/*-----------------------------------------------------------------------------
 *					CreatenxVectorObject		 2015- 9- 29*/
/** **/
/*---------------------------------------------------------------------------*/

PyObject* ISKPythonHelperFuncs::Create_nxVectorObject( const nxVector* v)
{
	PyObject*	x;
	PyObject*	y;
	PyObject*	z;

	x = PyFloat_FromDouble(v->X() );
	y = PyFloat_FromDouble(v->Y() );
	z = PyFloat_FromDouble(v->Z());
	return PyTuple_Pack( 3, x, y, z);
}


/*-----------------------------------------------------------------------------
 *					Create_ISKBasisObject		 2015- 9- 29*/
/** **/
/*---------------------------------------------------------------------------*/

PyObject* ISKPythonHelperFuncs::Create_ISKBasisObject(const ISKBasisDirection* basis) 
{
	PyObject*				basis_object;

	if (m_basis_typeobject == nullptr)
	{
		PyStructSequence_Desc	desc;
		PyStructSequence_Field	fields[4];
		fields  [0].name = "propagation";
		fields  [0].doc  = "The unit vector propagation direction of the ray";
		fields  [1].name = "theta";
		fields  [1].doc  = "The theta direction unit vector perpendicular to the ray direction";
		fields  [2].name = "phi";
		fields  [2].doc  = "The phi direction unit vector perpendicular to the ray direction";
		fields  [3].name = nullptr;
		fields  [3].doc  = nullptr;

		desc.name   = "ISKBasisDirection";
		desc.doc    = "The three unit vectors defining the ray. Note that (theta x phi => propagation).";
		desc.fields = fields;
		desc.n_in_sequence = 3;
		m_basis_typeobject = PyStructSequence_NewType( &desc );
	}
	basis_object = PyStructSequence_New( m_basis_typeobject );
	PyStructSequence_SetItem(basis_object, 0, Create_nxVectorObject( &basis->Propagation()) );
	PyStructSequence_SetItem(basis_object, 1, Create_nxVectorObject( &basis->Theta()) );
	PyStructSequence_SetItem(basis_object, 2, Create_nxVectorObject( &basis->Phi()) );
	return basis_object;
}


/*-----------------------------------------------------------------------------
 *					Create_ISK_IQUVObject		 2015- 9- 29*/
/** **/
/*---------------------------------------------------------------------------*/

PyObject* ISKPythonHelperFuncs::Create_ISK_IQUVObject(const IQUV* radiance) 
{
	PyObject*				iquv_object;

	if (m_iquv_typeobject == nullptr)
	{
		PyStructSequence_Desc	desc;
		PyStructSequence_Field	fields[5];
		fields  [0].name = "I";
		fields  [0].doc  = "The I component of the stokes vector";
		fields  [1].name = "Q";
		fields  [1].doc  = "The Q component of the stokes vector";
		fields  [2].name = "U";
		fields  [2].doc  = "The U component of the stokes vector";
		fields  [3].name = "V";
		fields  [3].doc  = "The V component of the stokes vector";
		fields  [4].name = nullptr;
		fields  [4].doc  = nullptr;

		desc.name   = "IQUV";
		desc.doc    = "The full radiance 4 element vector.";
		desc.fields = fields;
		desc.n_in_sequence = 4;
		m_iquv_typeobject = PyStructSequence_NewType( &desc );
	}
	iquv_object = PyStructSequence_New( m_iquv_typeobject );
	PyStructSequence_SetItem(iquv_object, 0, PyFloat_FromDouble(radiance->I) );
	PyStructSequence_SetItem(iquv_object, 1, PyFloat_FromDouble(radiance->Q) );
	PyStructSequence_SetItem(iquv_object, 2, PyFloat_FromDouble(radiance->U) );
	PyStructSequence_SetItem(iquv_object, 3, PyFloat_FromDouble(radiance->V) );
	return iquv_object;
}


/*-----------------------------------------------------------------------------
 *					Create_ISKStokesVectorObject		 2015- 9- 29*/
/** **/
/*---------------------------------------------------------------------------*/

PyObject* ISKPythonHelperFuncs::Create_ISKStokesVectorObject(const ISKStokesVector* radiance) 
{
	PyObject*				istokes_object;

	if (m_istokes_typeobject == nullptr)
	{
		PyStructSequence_Desc	desc;
		PyStructSequence_Field	fields[3];
		fields  [0].name = "stokes";
		fields  [0].doc  = "The IQUV stokes vector";
		fields  [1].name = "basis";
		fields  [1].doc  = "The basis vectors of the stokes vector";
		fields  [2].name = nullptr;
		fields  [2].doc  = nullptr;

		desc.name   = "ISKStokesVector";
		desc.doc    = "The Stokes Vector object";
		desc.fields = fields;
		desc.n_in_sequence = 2;
		m_istokes_typeobject = PyStructSequence_NewType( &desc );
	}
	istokes_object = PyStructSequence_New( m_istokes_typeobject );
	PyStructSequence_SetItem(istokes_object, 0, Create_ISK_IQUVObject( &radiance->Stokes()) );
	PyStructSequence_SetItem(istokes_object, 1, Create_ISKBasisObject( &radiance->Basis())  );
	return istokes_object;
}


/*-----------------------------------------------------------------------------
 *					CopyRadianceArray		 2015- 9- 29*/
/** **/
/*---------------------------------------------------------------------------*/

PyObject* ISKPythonHelperFuncs::CopyRadianceArray( ISKStokesVector* radiance, int numlinesofsight, int numwavelens)
{
	PyObject*				full_object;
	PyObject*				onewavelenobject;
	ISKStokesVector*		radianceptr = radiance;

	full_object = PyTuple_New( numwavelens );
	for (Py_ssize_t iw = 0; iw < (Py_ssize_t)numwavelens; iw++)
	{
		onewavelenobject = PyTuple_New( numlinesofsight );
		for (Py_ssize_t pos = 0; pos < (Py_ssize_t)numlinesofsight; pos++)
		{
			PyTuple_SetItem(onewavelenobject, pos, Create_ISKStokesVectorObject( radianceptr));
			radianceptr++;
		}
		PyTuple_SetItem( full_object, iw, onewavelenobject );
	}
	return full_object;
}


/*-----------------------------------------------------------------------------
 *					CopyRadianceArray		 2015- 9- 29*/
/** **/
/*---------------------------------------------------------------------------*/

PyObject* ISKPythonHelperFuncs::CopyRadianceScalarArray( const double* radiance, int numlinesofsight, int numwavelens)
{
	PyObject*				full_object;
	PyObject*				onewavelenobject;
	const double*					radianceptr = radiance;

	full_object = PyTuple_New( numwavelens );
	for (Py_ssize_t iw = 0; iw < (Py_ssize_t)numwavelens; iw++)
	{
		onewavelenobject = PyTuple_New( numlinesofsight );
		for (Py_ssize_t pos = 0; pos < (Py_ssize_t)numlinesofsight; pos++)
		{
			PyTuple_SetItem(onewavelenobject, pos, PyFloat_FromDouble( radianceptr) );
			radianceptr++;
		}
		PyTuple_SetItem( full_object, iw, onewavelenobject );
	}
	return full_object;
}

#endif


