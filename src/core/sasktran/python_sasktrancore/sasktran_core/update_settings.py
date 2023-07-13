import os
import os.path
import glob
import sys
import urllib.request
import tarfile
from . import hapi
from . hitran_manager import Hitran_HAPI_Manager
from sasktranif.sasktranif_registry import  SasktranifRegistry

#------------------------------------------------------------------------------
#           _rlinput
#------------------------------------------------------------------------------

def _rlinput(prompt, default =''):

    question = '{} [{}] ?:'.format(prompt,default)
    line = input(question)
    value = line.strip()
    if (len(value) < 1):
        line = default
    return line

#------------------------------------------------------------------------------
#           update_mie_aerosol_cache
#------------------------------------------------------------------------------

def _update_folder( prompt : str, skregistry: SasktranifRegistry, keyindexstring : str, response = None ):

    key,idx = keyindexstring.rsplit('/', maxsplit=1)
    idx = idx.lower()
    folder   = skregistry.subkey(key)[idx]
    if (response is None):
        response = _rlinput("Enter the {} ".format(prompt), default = folder)
    if sys.platform == 'win32':
        response = response. replace('\\', '/')
    skregistry.subkey(key)[idx] = response

# ------------------------------------------------------------------------------
#           update_baum_folder
# ------------------------------------------------------------------------------

def update_baum_folder( folder: str, skregistry: SasktranifRegistry):
    _update_folder( '', skregistry, 'Software/USask-ARG/skOpticalProperties/Baum_IceCrystals/Storage/2014Database', response=folder)

#------------------------------------------------------------------------------
#           update_hitran_folder
#------------------------------------------------------------------------------

def update_hitran_folder( folder: str, skregistry: SasktranifRegistry):
    _update_folder( '', skregistry, 'Software/USask-ARG/skOpticalProperties/Hitran/BaseDirectory', response=folder)

#------------------------------------------------------------------------------
#           _download_hitran_with_hapi
#------------------------------------------------------------------------------

def _download_hitran_with_hapi( hapifolder:str):

    if (sys.platform == 'win32'):
        hapifolder = hapifolder.replace('/', '\\')
    hitran = Hitran_HAPI_Manager( hapifolder )
    hitran.make_molparam_file()
    hitran.load_hitran_lines()

# ------------------------------------------------------------------------------
#           _delete_hitran_cache
# ------------------------------------------------------------------------------

def _delete_hitran_cache(hapifolder:str):

    if (sys.platform == 'win32'):
        hapifolder = hapifolder.replace('/', '\\')

    cachefiles = os.path.join( hapifolder, 'cache', "*.bin")
    listofiles = glob.glob( cachefiles )
    for f in listofiles:
        print('Deleting HITRAN cache file: {:s}'.format(f))
        os.remove(f)
    print('Hitran Cache files are deleted. They will be automaticaly regenerated as needed')

#------------------------------------------------------------------------------
#           _download_hitran_with_hapi
#------------------------------------------------------------------------------

def _download_pratmo_from_arg( pratmofolder:str):

    if (sys.platform == 'win32'):
        pratmofolder = pratmofolder.replace('/', '\\')

    if not os.path.exists(pratmofolder): os.makedirs(pratmofolder)
    filename = os.path.join( pratmofolder, 'skClimatology_PratmoNO2.bin')
    try:
        file,mes = urllib.request.urlretrieve( 'https://arg.usask.ca/sasktranfiles/skClimatology_PratmoNO2.bin', filename=filename )
        ok = os.path.exists(filename)
        print('SUCCESS. Pratmo climatology file [{}] has been downloaded from arg.usask.ca'.format(filename))
    except Exception as e:
        ok = False
        print(e)
        print('Error downloading file skClimatology_PratmoNO2.bin from arg.usask.ca')

#------------------------------------------------------------------------------
#           _download_baum_from_arg
#------------------------------------------------------------------------------

def _download_baum_from_arg( baumfolder:str):

    if (sys.platform == 'win32'):
        pratmofolder = baumfolder.replace('/', '\\')

    if not os.path.exists(baumfolder): os.makedirs(baumfolder)
    filename = os.path.join( baumfolder, 'AggregateSolidColumns_SeverelyRough_AllWavelengths_FullPhaseMatrix.nc')
    names = ['AggregateSolidColumns_SeverelyRough_AllWavelengths_FullPhaseMatrix.nc',
             'GeneralHabitMixture_SeverelyRough_AllWavelengths_FullPhaseMatrix.nc',
             'SolidColumns_SeverelyRough_AllWavelengths_FullPhaseMatrix.nc',
             'Baum_jsqrt_2014.pdf']
    try:
        for name in names:
            filename = os.path.join( baumfolder, name)
            print('Downloading {} to {}'.format(name, filename))
            file,mes = urllib.request.urlretrieve( 'https://arg.usask.ca/sasktranfiles/BaumIceCrystals/'+name, filename=filename)
            ok = os.path.exists(filename)
        print('SUCCESS. The BAUM Ice Crystal database has been downloaded from arg.usask.ca')
    except Exception as e:
        ok = False
        print(e)
        print('Error downloading the BAUM Ice Crystal database from arg.usask.ca')

#------------------------------------------------------------------------------
#           _download_grib_from_arg
#------------------------------------------------------------------------------

def _download_grib_from_arg(  skregistry, definitions_folder_key: str, samples_folder_key : str ):

    definitions_folder = skregistry.subkey(definitions_folder_key)
    grib_base_folder, defname = os.path.split(definitions_folder)
    if not os.path.exists(grib_base_folder): os.makedirs(grib_base_folder)

    definitions_folder = os.path.join( grib_base_folder, 'definitions' )
    samples_folder     = os.path.join( grib_base_folder, 'samples')

    tarfilename ='grib_definitions.tar.gz'
    if  os.path.exists(tarfilename):
        os.unlink(tarfilename)
    try:
        file, mes = urllib.request.urlretrieve('https://arg.usask.ca/sasktranfiles/grib_definitions.tar.gz',  filename=tarfilename)
        ok = os.path.exists(tarfilename)
        print('SUCCESS. The grib definitions and sample file [{}] has been downloaded from arg.usask.ca'.format(tarfilename))
    except Exception as e:
        ok = False
        print(e)
        print('Error downloading file grib_definitions.tar.gz from arg.usask.ca')
    if (ok):

        print('Extracting grib definitions and samples to folder [{}] '.format(grib_base_folder))
        if not os.path.exists(grib_base_folder): os.makedirs(grib_base_folder)
        m_tar = tarfile.open(tarfilename)
        m_tar.extractall(grib_base_folder)
        m_tar.close()
        _update_folder( '',skregistry, definitions_folder_key, response= definitions_folder)
        _update_folder( '',skregistry, samples_folder_key,     response= samples_folder)
    if  os.path.exists(tarfilename):
        os.unlink(tarfilename)


#------------------------------------------------------------------------------
#           display_menu
#------------------------------------------------------------------------------

def _display_menu(skregistry:SasktranifRegistry):

    id = 1

    folderarray = [  ['folder for the Mie aerosol cache',           'Software/USask-ARG/skOpticalProperties/MieAerosol/Aerosol_Cache_Directory'],
                     ["folder for the Mie ice crystal cache",       'Software/USask-ARG/skOpticalProperties/Mie_IceCrystals/Ice_Crystal_Cache_Directory'],
                     ["folder for the Baum 2014 ice crystal cache", 'Software/USask-ARG/skOpticalProperties/Baum_IceCrystals/Storage/2014Database'],
                     ["folder for Hitran cross-sections",           'Software/USask-ARG/skOpticalProperties/Hitran/BaseDirectory'],
                     ["folder for Hitran caches",                   'Software/USask-ARG/skOpticalProperties/Hitran/SpectralLineCacheDir'],
                     ["folder for the Pratmo climatology",          'Software/USask-ARG/Climatology/Pratmo/Storage'],
                     ["folder for the ECMWF Era Interim database",  'Software/USask-ARG/Climatology/ECMWF/Storage/EraInterimDir'],
                     ["folder containing GRIB definitions",         'Software/USask-ARG/Grib/definitions_folder'],
                     ["folder containing GRIB samples",             'Software/USask-ARG/Grib/samples_folder']
                  ]

    while (id > 0):
        print("\nPackage options for the Sasktran Core Components")
        print("------------------------------------------------\n")
        for optionid in range(len(folderarray)):
            print("     {:2d}) {:<44s}-> {}".format( optionid+1, folderarray[optionid][0], skregistry.subkey(folderarray[optionid][1]) ))

        hitranoption = len(folderarray) + 1
        hitrancache  = hitranoption + 1
        pratmooption = hitranoption + 2
        griboption   = hitranoption + 3
        baumoption   = hitranoption + 4
        print("     {:2d}) download Hitran database from hitran.org using HAPI. This may take a long time (~ 1 hour) ".format(hitranoption))
        print("     {:2d}) delete and reset the current Hitran cache".format(hitrancache))
        print("     {:2d}) download pratmo climatology from arg.usask.ca ".format(pratmooption))
        print("     {:2d}) download Grib definitions and samples from arg.usask.ca ".format(griboption))
        print("     {:2d}) download Baum Ice Crystal database from arg.usask.ca ~500 MB".format(baumoption))
        print("     0) Save and Exit")
        print("    -1) Quit and do not save")
        print("\n")
        selection = input('Enter your selection: ')
        try:
            id = int( selection)
        except:
            id = 1
        if (id > 0):
            if id <= len(folderarray):
                _update_folder( folderarray[id-1][0], skregistry, folderarray[id-1][1])
            elif id == hitranoption:
                _download_hitran_with_hapi(skregistry.subkey('software/usask-arg/skopticalproperties/hitran/basedirectory'))
            elif id == hitrancache:
                _delete_hitran_cache(skregistry.subkey('software/usask-arg/skopticalproperties/hitran/basedirectory'))
            elif id == pratmooption:
                _download_pratmo_from_arg(skregistry.subkey('software/usask-arg/climatology/pratmo/storage'))
            elif id == griboption:
                _download_grib_from_arg(skregistry, 'software/usask-arg/grib/definitions_folder','software/usask-arg/grib/samples_folder')
            elif id == baumoption:
                _download_baum_from_arg(skregistry.subkey('software/usask-arg/skopticalproperties/baum_icecrystals/storage/2014database'))

            else:
                print("Unknown option {:d}".format(id) )
        else:
            if (id == 0):
                skregistry.save()
                print("The registry settings have been saved.")
            else:
                print("The registry settings have NOT been saved.")
    return id

#------------------------------------------------------------------------------
#           test_hapi
#------------------------------------------------------------------------------

def _test_hapi():
    hapi.db_begin(r'C:\temp')
    hapi.fetch_by_ids('H2O', [1, 2, 3, 4, 5, 6, 7], 0.0, 1.0E6)
    #miekeys = skregistry.subkey('software/usask-arg/skopticalproperties/mieaerosol')
    #miekeys['aerosol_cache_directory'] = input('Enter folder for the Mie aerosol cache: ')
    #print('The value is ', miekeys['aerosol_cache_directory'] )
    
    #hitran = Hitran_HAPI_Manager( r'\\UTLS\OsirisDataProducts\hitranhapi')
    #hitran.make_molparam_file()
    #hitran.load_hitran_lines()
    
    
    
    #for mol_id,entries in molecules.items():
    #    outputfile =  "{:02d}_hit".format(mol_id)
    #    fullname  = hitranfolder + os.sep + outputfile
    #    if not os.path.exists(fullname+".data"):
    #        print( "Fetching molecule ", mol_id, " to file " ,fullname+".data", " with entries ", entries )
    #        hapi.fetch_by_ids(outputfile, entries,0.0, 1.0E6) # H2O
    #print("Done")

#------------------------------------------------------------------------------
#           def modify_sasktran_core_registry_settings():
#------------------------------------------------------------------------------

def modify():
    skregistry = SasktranifRegistry()
    id = _display_menu(skregistry)
    print('Finished modifying sasktran_core registry settings')

if __name__ == "__main__":
    skregistry = SasktranifRegistry()
    args = sys.argv
    if (len(args) > 1):
        command = args[1].lower()
        if command == 'download_pratmo':
            location = skregistry.subkey('software/usask-arg/skopticalproperties/hitran/basedirectory')
            print('Downloading Pratmo climatology to ',location)
            _download_pratmo_from_arg(location)
        elif command=='download_hitran':
            location = skregistry.subkey('software/usask-arg/skopticalproperties/hitran/basedirectory')
            print('Downloading HITRAN database to ', location)
            _download_hitran_with_hapi(location)
        elif command=='download_grib':
            print('Downloading Grib definitions and samples database to ')
            _download_grib_from_arg( skregistry, 'software/usask-arg/grib/definitions_folder', 'software/usask-arg/grib/samples_folder')
        elif command =='set_baum_folder':
            folder = args[2]
            update_baum_folder(folder, skregistry)
            skregistry.save()
            print('Updated Baum folder to <{:s}>'.format(folder))
        elif command =='set_hitran_folder':
            folder = args[2]
            update_hitran_folder(folder, skregistry)
            skregistry.save()
            print('Updated HITRAN folder to <{:s}>'.format(folder))
        else:
            print('Unrecognized command')
            exit(1)
