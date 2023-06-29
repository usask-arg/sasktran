import os
import os.path
import sys
import glob
import logging
from appdirs import AppDirs

def _createDefaultRegistryEntries(registry):

    defaulthapi       = r'\\UTLS\SasktranFiles\hitranhapi'
    defaultno2        = r'\\UTLS\SasktranFiles\pratmo'
    defaultbaum       = r'\\UTLS\SasktranFiles\BaumIceCrystals'
    defaultera        = r'\\UTLS\ecmwf'
    defaultgribdef    = r'\\UTLS\SasktranFiles\grib\definitions'
    defaultgribsamp   = r'\\UTLS\SasktranFiles\grib\samples'

    dirs              = AppDirs('sasktran','usask-arg')
    cachedir          = dirs.user_data_dir
    mieaerosolcache   = os.path.join( cachedir, 'mie_aerosol')                                 # Caches go into the local python environment
    iceaerosolcache   = os.path.join( cachedir, 'ice_aerosol')
    baumaerosolcache  = defaultbaum     if os.path.exists(defaultbaum)     else os.path.join( cachedir, 'baum')
    no2file           = defaultno2      if os.path.exists(defaultno2)      else os.path.join( cachedir, 'pratmo')
    hitran            = defaulthapi     if os.path.exists(defaulthapi)     else os.path.join( cachedir, 'hitran')
    hitrancache       =                                                         os.path.join( cachedir, 'hitrancache')
    erainterim        = defaultera      if os.path.exists(defaultera)      else os.path.join( cachedir, 'erainterim')
    gribdef           = defaultgribdef  if os.path.exists(defaultgribdef)  else os.path.join( cachedir, 'grib','definitions')
    gribsamp          = defaultgribsamp if os.path.exists(defaultgribsamp) else os.path.join( cachedir, 'grib','samples')

    registry.setkey( 'Software/USask-ARG/skOpticalProperties/MieAerosol/Aerosol_Cache_Directory',          mieaerosolcache)
    registry.setkey( 'Software/USask-ARG/skOpticalProperties/Mie_IceCrystals/Ice_Crystal_Cache_Directory', iceaerosolcache)
    registry.setkey( 'Software/USask-ARG/skOpticalProperties/Baum_IceCrystals/Storage/2014Database',       baumaerosolcache)
    registry.setkey( 'Software/USask-ARG/skOpticalProperties/Hitran/BaseDirectory',                        hitran)
    registry.setkey( 'Software/USask-ARG/skOpticalProperties/Hitran/SpectralLineCacheDir',                 hitrancache)
    registry.setkey( 'Software/USask-ARG/Climatology/Pratmo/Storage',                                      no2file)
    registry.setkey( 'Software/USask-ARG/Climatology/ECMWF/Storage/EraInterimDir',                         erainterim)
    registry.setkey( 'Software/USask-ARG/Grib/definitions_folder',                                         gribdef)
    registry.setkey( 'Software/USask-ARG/Grib/samples_folder',                                             gribsamp)
    registry.save()

#------------------------------------------------------------------------------
#           check_installation_called_from_sasktranif
#------------------------------------------------------------------------------

def check_installation_called_from_sasktranif( SKTRAN_IFCreateRegistryEntriesForDLL, registry ):

    dirname = os.path.dirname(__file__)
    firstimefile = os.path.join( dirname, 'sasktran_core_firsttime.sktran' )
    if os.path.exists(firstimefile):
        logging.info("sasktran_core, first time setup: Initiliazing Sasktran Core Components registry settings")
        _createDefaultRegistryEntries(registry)
        dllnames = glob.glob( dirname + os.sep + "*_sasktran_core_internals*")            # Search that directory for the core components DLL or shareable object
        dllname = dllnames[0]                                                        # We must find at least 1 (or fail)
        SKTRAN_IFCreateRegistryEntriesForDLL(dllname, dllname)                       # param[0] =  full path name of this DLL.
        os.unlink(firstimefile)
        logging.info("First Time Setup: Sasktran Core Components default initialization is complete.")
        logging.info("To modify sasktran_core settings and download various databases, please use:\n\n >>> import sasktran_core.update_settings\n >>> sasktran_core.update_settings.modify()\n\n ")
 