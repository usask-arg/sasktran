
import os
import sys
from shutil import copyfile, rmtree
import versioneer
from setuptools import setup
from setuptools import find_packages
from setuptools.dist import Distribution
from distutils.core import Extension
import numpy as np  
import glob
from copy import copy
from pathlib import Path

def get_cplus_paths():
    cplus_include_path = []
    cplus_library_path = []
    cplus_include_pathstr = os.environ.get("SASKTRAN_IF_INCLUDE_PATH")
    cplus_library_pathstr = os.environ.get("SASKTRAN_IF_LIBRARY_PATH")

    print("\n**************************** READ THIS NOTICE ****************************************\n")
    print("CPLUS_INCLUDE_PATH = {}".format(os.environ.get("SASKTRAN_IF_INCLUDE_PATH")))
    print("CPLUS_LIBRARY_PATH = {}".format(os.environ.get("SASKTRAN_IF_LIBRARY_PATH")))
    print("If these folders are incorrect you should manually edit your environment variables")
    print("\n**************************************************************************************\n")

    cplus_include_pathstr = cplus_include_pathstr.replace('"', '')
    cplus_library_pathstr = cplus_library_pathstr.replace('"', '')

    if (cplus_include_pathstr is not None): cplus_include_path = list( filter( bool, cplus_include_pathstr.split(';')) )
    if (cplus_library_pathstr is not None): cplus_library_path = list( filter( bool , cplus_library_pathstr.split(';')) )

    cplus_include_path = [Path(p.strip()).as_posix() for p in cplus_include_path if not p.isspace()]
    cplus_library_path = [Path(p.strip()).as_posix() for p in cplus_library_path if not p.isspace()]

    # Because of a cmake bug these paths could actually be full paths to a file, so check for that
    for idx, p in enumerate(cplus_include_path):
        if Path(p).is_file():
            cplus_include_path[idx] = Path(p).parent.as_posix()
        for idx, p in enumerate(cplus_library_path):
            if Path(p).is_file():
                cplus_library_path[idx] = Path(p).parent.as_posix()


    return cplus_include_path, cplus_library_path


package_data = {}
data_files = []
extra_objects = []

cplus_include_path, cplus_library_path = get_cplus_paths()

include_dirs      =[  r'../Repos_SasktranIF/includes',                                                              # setup the include folders for the compilation
                      np.get_include(),
                      r'../Repos_BaseCode/nxbase'
                      ] +  cplus_include_path

if sys.platform =='win32':                                                                              # if we are on a windows Machine
    library_dirs      =[ r'..\fortran_libraries\lib\Windows_x64',
                    ] +  cplus_library_path
    libraries         =['sasktranif_Release','User32','Advapi32','Ole32','OleAut32','nxbase_Release', 'yaml-cpp', 'bcrypt']
    extra_compile_args=['/MD']
    extra_link_args   =[]
else:
    lib_path = os.environ.get('LIBRARY_PATH')
    library_dirs = cplus_library_path

    if sys.platform == 'darwin':
        libraries = ['sasktranif', 'nxbase', 'yaml-cpp', 'boost_system', 'boost_thread', 'boost_filesystem']
    else:
        libraries = ['sasktranif', 'nxbase', 'yaml-cpp', 'boost_system', 'boost_thread', 'boost_filesystem', 'rt']

    # We want to statically link every library we can, so check to see if any static lib exists in the search paths
    for directory in library_dirs:
        for lib in copy(libraries):
            possible_static_lib = Path(directory).joinpath('lib' + lib + '.a')
            # print(possible_static_lib.as_posix())
            if possible_static_lib.exists():
                extra_objects.append(possible_static_lib.as_posix())
                libraries.remove(lib)

    extra_compile_args=['-std=gnu++11', '-fvisibility=hidden','-fPIC' ]
    if sys.platform == 'darwin':
        # mac
        extra_link_args = []
    else:
        # Linux
        extra_link_args   =[]
    data_files        =[]

try:
    os.mkdir('sasktranif')
except OSError:
    # Directory already exists
    pass

copyfile( r'../Repos_SasktranIF/python_sasktranif/sasktranif/__init__.py', 'sasktranif/__init__.py' )
copyfile( r'../Repos_SasktranIF/python_sasktranif/sasktranif/sasktranif_registry.py', 'sasktranif/sasktranif_registry.py' )
copyfile( r'../Repos_SasktranIF/swig/sasktranif.py', 'sasktranif/sasktranif.py')                               # Copy the swig generated file
copyfile( r'../Repos_SasktranIF/sources/pythonlogger.cpp','pythonlogger.cpp')                                            # Copy the sources over. Note that relative links to the source code
copyfile( r'../Repos_SasktranIF/swig/sasktranif_wrap.cxx','sasktranif_wrap.cxx')                                         # failed when building on Linux

# Sasktran Core Components parts

try:
    os.mkdir('sasktran_core')
except OSError:
    # Directory already exists
    pass

data_files = [ ('share/usask-arg/sasktran/installed_modules',['sasktran_core.sktran'])]

copyfile(r'../Repos_SasktranV3/python_sasktrancore/sasktran_core/sasktran_core_firsttime.sktran', 'sasktran_core/sasktran_core_firsttime.sktran')
copyfile(r'../Repos_SasktranV3/python_sasktrancore/sasktran_core/__init__.py', 'sasktran_core/__init__.py')
copyfile(r'../Repos_SasktranV3/python_sasktrancore/sasktran_core/hapi.py', 'sasktran_core/hapi.py')
copyfile(r'../Repos_SasktranV3/python_sasktrancore/sasktran_core/hitran_manager.py', 'sasktran_core/hitran_manager.py')
copyfile(r'../Repos_SasktranV3/python_sasktrancore/sasktran_core/update_settings.py', 'sasktran_core/update_settings.py')

copyfile(r'../Repos_SasktranV3/python_sasktrancore/sasktran_core.sktran', 'sasktran_core.sktran')



if sys.platform =='win32':                                                                              # if we are on a windows Machine
    print('Copying Fortran and Sasktran core DLLs into package')
    copyfile(r'..\Repos_SasktranV3\python_sasktrancore\sasktran_core\_sasktran_core_internals.dll', 'sasktran_core/_sasktran_core_internals.dll')

    dllnames = glob.glob(  os.path.join( 'sasktran_core','_sasktran_core_internals.dll') )
    dllname = os.path.basename(dllnames[0])                                                                                    # This file is built before calling setup.py by the SasktranCoreComponents/buildarglibraries.bat
    copyfile( r'..\fortran_libraries\lib\Windows_x64\wiscombemie.dll',        r'sasktran_core\wiscombemie.dll')
    copyfile( r'..\fortran_libraries\lib\Windows_x64\tmatrixrandomep.dll',    r'sasktran_core\tmatrixrandomep.dll')
    copyfile( r'..\fortran_libraries\lib\Windows_x64\msis90e.dll',            r'sasktran_core\msis90e.dll')
    copyfile( r'..\fortran_libraries\lib\Windows_x64\hitran_tips.dll',        r'sasktran_core\hitran_tips.dll')
    copyfile( r'..\fortran_libraries\lib\Windows_x64\netlib.dll',             r'sasktran_core\netlib.dll')
    copyfile( r'..\fortran_libraries\lib\Windows_x64\blaslapack.dll',         r'sasktran_core\blaslapack.dll')

    package_data  = {'sasktran_core': [ 'wiscombemie.dll',                             # these are necessary DLL's for the windows build
                                         'tmatrixrandomep.dll',
                                         'msis90e.dll',
                                         'hitran_tips.dll',
                                         'netlib.dll',
                                         'blaslapack.dll',
                                         dllname,
                                         'sasktran_core_firsttime.sktran'               # this file tells sasktran_core that it is the firsttime it is loaded. This will trigger initialization of the module.
                                       ]
                     }
    print('Copying completed')

else:
    if sys.platform == 'darwin':
        soname = 'lib_sasktran_core_internals.dylib'
    else:
        soname = 'lib_sasktran_core_internals.so'

    copyfile('../Repos_SasktranV3/python_sasktrancore/sasktran_core/{}'.format(soname), 'sasktran_core/{}'.format(soname))
    package_data      ={'sasktran_core': [ 'sasktran_core_firsttime.sktran', soname ] }                # this file tells sasktran_core that it is the firsttime it is loaded. This will trigger initialization of the module.

extension_module = Extension(                                                                           # specify that we are creating a python extension
     name              = 'sasktranif._sasktranif',
     sources           =['pythonlogger.cpp',
                         'sasktranif_wrap.cxx'],
     include_dirs      =include_dirs,
     library_dirs      =library_dirs,
     libraries         =libraries,
     extra_compile_args=extra_compile_args,
     extra_link_args   =extra_link_args,
     extra_objects=extra_objects
)

setup(
    name='sasktran',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    url='https://arg.usask.ca/docs/sasktran/',
    license='MIT',
    author='Daniel Zawada',
    author_email='daniel.zawada@usask.ca',
    description='',
    include_package_data=True,
    package_data={**package_data, **{'sasktran': ['aband/data/*']}},
    data_files=data_files,
    ext_modules=[extension_module],
    install_requires=[
        'numpy',
        'pyyaml',
        'urllib3',
        'appdirs',
        'packaging'
    ]
)


os.remove('pythonlogger.cpp')                                                                           
os.remove('sasktranif_wrap.cxx')
os.remove('sasktran_core.sktran')
rmtree('sasktran_core')
rmtree('sasktranif')