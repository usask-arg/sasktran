import os
import os.path
import sys
import sysconfig
import logging
from . import sasktranif
from .sasktranif import ISKClimatology, ISKOpticalProperty, ISKEngine, ISKBrdf, ISKSolarSpectrum, ISKEmission, ISKGeodetic, SKTRAN_IFSetRegistryDirectory
from .sasktranif_registry import SasktranifRegistry


#-----------------------------------------------------------------------------------------
#                       Start of main init scrip
#  This initializes (just) the SasktranIF module (_sasktranif.pyd) to use the python
#  environment registry. All other Sasktran C++ modules will be automatically initialized to use
#  the python environment registry when they are loaded by the _sasktranif.pyd code. For this to work
#  properly each module must invoke sasktranif.sasktranif.SKTRAN_IFCreateRegistryEntriesForDLL(dllname, paramstr)
#  as part of the installation procedure. An example of configuring the regsitry can be found in installcorecomponents.py
#
#-----------------------------------------------------------------------------------------

registry_directory = os.path.join(sysconfig.get_path('data'), 'share', 'usask-arg', 'registry', 'sasktranif')
package_directory  = os.path.join(sysconfig.get_path('data'), 'share', 'usask-arg', 'sasktran', 'installed_modules')
SKTRAN_IFSetRegistryDirectory(registry_directory)                                                                              # Setup each sasktran module so it uses a text based registry within its own area. This value is stored internally in the sasktranif C++ code
registry = SasktranifRegistry()
if os.path.exists(package_directory):
    for install_file in os.listdir(package_directory):                                                                          # Now look for all of the installed sasktran modules
            name,ext = os.path.splitext( install_file)                                                                          # so get the name and extension of each file in that directory
            if (ext == '.sktran'):                                                                                              # if the extension is .sktran
                try:
                    module = __import__(name)                                                                                   # then import the module given by the name. This usually does very little, (it does not load teh DLL for example)
                    ok = 'check_installation_called_from_sasktranif' in dir(module)                                                 # then we need to check the module has completed its initialization, by writing entries into the registry through sasktranif
                    if (ok):                                                                                                        # so if we have that method
                        module.check_installation_called_from_sasktranif( sasktranif.SKTRAN_IFCreateRegistryEntriesForDLL, registry )         # check the module has completed its installation. Pass the special function in , as sasktranif is still in __init__func and we will get issues if we import form the module
                    else:                                                                                                           # if the method is not available then log a message.
                        logging.info("Module {:s} not properly imported by sasktranif as it does not have a method <check_installation_called_from_sasktranif> to check the modules initialization".format(name))
                except Exception as e:
                    print( 'Module {:s} was not loaded and threw an exception [{}]. That indicates an installation problem'.format(name, e) )
                    logging.info( 'Module {:s} was not loaded and threw an exception [{}]. That indicates an installation problem'.format(name,e))
    #except:
        #    raise ValueError('Error loading module from {%s}'.format(install_file))
else:
   logging.warning("Note: SasktranIF cannot find any Sasktran components to load. Thats unusual and suspects an installation error has occurred")
