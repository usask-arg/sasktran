from typing import List
import sys
import os
import os.path
from . import hapi

#------------------------------------------------------------------------------
#           class Hitran_HAPI_Manager
#------------------------------------------------------------------------------

class Hitran_HAPI_Manager:

    #------------------------------------------------------------------------------
    #           __init__
    #------------------------------------------------------------------------------

    def __init__(self, hitran_base_directory : str = None):

        self.m_hapi_initialized = False
        if (hitran_base_directory is None):
            if sys.platform == 'win32':
                hitran_base_directory = os.path.normpath( os.environ['LOCALAPPDATA'] + os.sep + '.HitranHapi' )
            else:
                hitran_base_directory = os.path.normpath(os.environ['HOME'] + os.sep + '.HitranHapi')
        self.m_hitran_base_directory = hitran_base_directory

    #------------------------------------------------------------------------------
    #           start_hapi
    #------------------------------------------------------------------------------

    def _start_hapi(self):

        if not self.m_hapi_initialized:
            basedir = os.path.join(self.m_hitran_base_directory, 'By-Molecule' , 'Uncompressed-files')
            if not os.path.exists(basedir): os.makedirs(basedir)
            hapi.db_begin(basedir)
            self.m_hapi_initialized = True

    #------------------------------------------------------------------------------
    #           load_molecules_and_isotopologues
    #------------------------------------------------------------------------------

    def _load_molecules_and_isotopologues(self):

        self._start_hapi()
        molecules     = {}
        isotopologues = {}
        MIDX = hapi.ISO_ID_INDEX['M']                                   # get the index of the molecule ID in ISO_ID        (typically this is always 0)
        IIDX = hapi.ISO_ID_INDEX['I']                                   # get the index of the Isotopologue ID in ISO_ID    (typically this is always 1)
        for global_id,k in hapi.ISO_ID.items():                         # for each molecule entry is global dictionary ISO_ID
            mol_id = k[MIDX]                                            # Get the molecule ID
            iso_id = k[IIDX]                                            # Getthe ISO ID. This number has chnaged from earlier version sof HITRAN
            entries = molecules.get(mol_id)                             # Get the listof isotopologues currently processed for this molecule
            abundance      = hapi.abundance(mol_id, iso_id)              # Get the abundance
            molecularmass  = hapi.molecularMass(mol_id, iso_id)          # and molecular mass
            qt             = 0.0
            gj             = 0
            thisentry      = [ iso_id, abundance, molecularmass, qt, gj, global_id]
            if entries is None:
                molecules    [mol_id] = [global_id]
                isotopologues[mol_id] = [thisentry]
            else:
                molecules    [mol_id].append(global_id)
                isotopologues[mol_id].append( thisentry )
        molecules = list(molecules.keys())
        molecules.sort()
        return molecules, isotopologues

    #------------------------------------------------------------------------------
    #           make_molparam_file
    #------------------------------------------------------------------------------

    def make_molparam_file( self, overwrite_existing = False ):
        """
        Makes a copy of molparam.txt using the hitran hapi.py interface. The sasktran C++ actually use this molparam.txt
        file rather than try to use the python hapi interface.

        Parameters
        ----------
        hitran_base_directory : str
            The base directory where Hitran files are stored.

        """
        molecules,isotopologues = self._load_molecules_and_isotopologues()
        basedir = os.path.normpath(self.m_hitran_base_directory  + os.sep + 'Global_Data')
        if not os.path.exists(basedir): os.makedirs(basedir)
        filename = os.path.normpath( basedir + os.sep + 'molparam.txt')
        if (os.path.exists(filename)) and not overwrite_existing:
            print('Skipping creating {} as file already exists'.format(filename))
        else:
            with open(filename, 'wt') as f:
                f.write( 'Molecule # Iso Abundance     Q(296K)      gj    Molar Mass(g)\n')
                for mol_id in molecules:
                    name = hapi.moleculeName(mol_id)
                    f.write('{:>7s} ({:1d})\n'.format(name,mol_id))
                    for entry in isotopologues[mol_id]:
                        iso_id, abundance, molecularmass, qt, gj, global_id = entry
                        f.write('          {:2d}  {:12.6e}  0.0000E+00    0     {:9.6f}\n'.format( iso_id, abundance,molecularmass))
            print("Created HITRAN file: ", filename)

    #------------------------------------------------------------------------------
    #           load_hitran_lines
    #------------------------------------------------------------------------------

    def load_hitran_lines(self, mol_id : List[int]  = None, overwrite_existing = False):

        molecules,isotopologues = self._load_molecules_and_isotopologues()
        basedir = 'By-Molecule' +os.sep +'Uncompressed-files'
        fulldir = os.path.normpath(self.m_hitran_base_directory + os.sep + basedir)
        if not os.path.exists( fulldir): os.makedirs(fulldir)
        if (mol_id is not None): molecules = mol_id
        for mol_id in molecules:
            isoid = []
            for entry in isotopologues[mol_id]:
                iso_id, abundance, molecularmass, qt, gj, global_id = entry
                isoid.append( global_id)
            name = hapi.moleculeName(mol_id)
            outputfile =  basedir + os.sep + name
            fullname = os.path.normpath(self.m_hitran_base_directory + os.sep+outputfile + '.data')

            if (os.path.exists(fullname)) and not overwrite_existing:
                print('Skipping HITRAN HAPI data for {} ({:1d}) as file {} already exists'.format(name,mol_id, outputfile))
            else:
                print("Fetching HITRAN HAPI data for {} ({:1d}) to file {}".format(name,mol_id, outputfile))
                for e in isoid:
                    print(name, e)
                try:
                    hapi.fetch_by_ids(name, isoid, 0.0, 1.0E6)
                except Exception as e:
                    print(e)



