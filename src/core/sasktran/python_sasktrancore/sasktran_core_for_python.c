
#include "Python.h"


/*
 This is a full-fledged python extension which is used to "fool" the setup.py
 to create the proeprp name of the binary wheel dirstributon. The real DLL/Shareable
 object is pre-built by scripts and installed as package data in setup.py.
 This avoids complications on Windows as Python really,really wants to build with the same
 compiler version as the Python interpreter (eg python 3.5 is Visual Studio 2015) but this
 creates problems for new sasktran engines which dont have VS2015 solutions.

 The real DLL has no python dependency and uses a C/C++ ABI to interface with sasktranif. It can be built
 with any compiler as long as its run-times are on the target installation python machine.

*/
/*---------------------------------------------------------------------------
 *						sasktran_core_methods
 *-------------------------------------------------------------------------*/

static struct PyMethodDef sasktran_core_methods[]=
{
	{ NULL,                      NULL, 0, NULL}
};


static struct PyModuleDef moduledef = {
		PyModuleDef_HEAD_INIT,
		"_dummy_sasktran_core",
		"Sasktran Core C++ Internals",
		-1,
		sasktran_core_methods
};


/*---------------------------------------------------------------------------
 *						initosl1
 *-------------------------------------------------------------------------*/

PyMODINIT_FUNC PyInit__dummy_sasktran_core(void)
{

    return PyModule_Create(&moduledef);
}

