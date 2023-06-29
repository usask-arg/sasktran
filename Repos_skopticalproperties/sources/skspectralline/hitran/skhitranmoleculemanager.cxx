#include <skopticalproperties21.h>
#include "nxbase_threads.h"


/*-----------------------------------------------------------------------------
 *					skHitranMoleculeManager::ToKey		2013-3-12*/
/** **/
/*---------------------------------------------------------------------------*/

size_t skHitranPartitionTableEntry::ToKey()
{
	return (m_isotopeid) + (m_moleculeid*1000000);
}


/*-----------------------------------------------------------------------------
 *					skHitranPartitionTableEntry::ReserveSpace		2013-3-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool skHitranPartitionTableEntry::ReserveSpace(size_t numtemperature)
{
	m_T.clear();
	m_Q.clear();
	m_T.reserve(numtemperature);
	m_Q.reserve(numtemperature);
	return true;
}


/*-----------------------------------------------------------------------------
 *					skHitranPartitionTableEntry::AddEntry		2013-3-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool skHitranPartitionTableEntry::AddEntry( double t, double q)
{
	m_T.push_back(t);
	m_Q.push_back(q);
	return true;
}


/*-----------------------------------------------------------------------------
 *					skHitranPartitionTableEntry::InternalPartition		2013-3-19*/
/** **/
/*---------------------------------------------------------------------------*/


double skHitranPartitionTableEntry::InternalPartition	( double T ) const
{
	double	value = 0.0;
	if (m_hapicompliant)
	{
		static std::mutex				fortran_mutexlock;
		std::unique_lock<std::mutex>	lock(fortran_mutexlock);	// The fortran code is not thread safe so protect
		int MOL = (int) m_moleculeid;								// The Hitran molecule number (typically 1-41) , eg 1 is H2O
		int ISO = (int) m_isotopeid;								// The isotope id. this is typically 3 or 4 digits like 121 but can be up to 4 digits. we allow up to 7
		double temp = T;
		double gi;
		double QT;
		BD_TIPS_2017( &MOL, &temp, &ISO, &gi, &QT);					// Call the TIPS fortran code supplied with HITRAN .
		value = QT;
	}
	else
	{
		value = nxLinearInterpolate::EvaluateYatX(T, m_T, m_Q, nxLinearInterpolate::ENUM_MISSINGVALUE, std::numeric_limits<double>::quiet_NaN());
	}
	return value;
}

/*-----------------------------------------------------------------------------
 *					skHitranMoleculeManagerNamePred					2013-3-13*/
/** \ingroup hitranoptpropinternals
 *	Small Predicate class used to search for molecules and isotopes in the
 *	skHitranMoleculeManager m_molecules array.  The class can match both
 *	chemical name and isotope id. It will ignore isotope id  is it is set to 0.
 **/
/*---------------------------------------------------------------------------*/

class skHitranMoleculeManagerNamePred
{
	private:
		std::string		m_name;
		size_t			m_isotopeid;

	public:
						skHitranMoleculeManagerNamePred( const char* name, size_t isotopeid)		{ m_name = name; m_isotopeid = isotopeid;}
		bool			operator()					   ( skHitranMoleculeManager::value_type pr)	{ return ((m_isotopeid ==0) || (pr.second.m_isotopeid == m_isotopeid)) && (pr.second.m_chemicalname == m_name);}
};

/*-----------------------------------------------------------------------------
 *					skHitranMoleculeManager::FindMoleculeEntry		2013-3-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool skHitranMoleculeManager::FindMoleculeEntry( const char* chemicalname, size_t isotopeid, skHitranPartitionTableEntry** entry )
{
	iterator	iter;
	bool		ok;
	nxString	buffer(chemicalname);

	buffer.MakeUpper();
	iter = std::find_if( m_molecules.begin(), m_molecules.end(), skHitranMoleculeManagerNamePred( buffer, isotopeid) );
	ok = !(iter == m_molecules.end());
	if (!ok)
	{
		*entry = NULL;
	}
	else
	{
		*entry = &(*iter).second;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skHitranMoleculeManager::FindMoleculeEntry		2013-3-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool skHitranMoleculeManager::FindMoleculeEntry( size_t moleculeid, size_t isotopeid, const skHitranPartitionTableEntry** entry ) const
{
	const_iterator						iter;
	bool							ok;
	skHitranPartitionTableEntry		dummy;
	size_t							key;

	dummy.m_moleculeid = moleculeid;
	dummy.m_isotopeid  = isotopeid;
	key                = dummy.ToKey();

	iter = m_molecules.find( key );
	ok = !(iter == m_molecules.end());
	if (!ok)
	{
		*entry = NULL;
	}
	else
	{
		*entry = &(*iter).second;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					IsotopeOrderBinaryPred		2013-3-14*/
/** A binary predicate function used to sort the isotope orders of a molecule
**/
/*---------------------------------------------------------------------------*/

static bool	IsotopeOrderBinaryPred( const skHitranPartitionTableEntry* a, const skHitranPartitionTableEntry* b )
{
	return (a->m_isotopeorder < b->m_isotopeorder);
}

/*-----------------------------------------------------------------------------
 *					skHitranMoleculeManager::FetchAllIsotopeEntries		2013-3-14*/
/** Fetch all of the isotope entries associated with a given molecule
 *	and sorted into ascending isotope abundance. (In reality they are
 *	in the same order as written in file molparam.txt). 
 **/
/*---------------------------------------------------------------------------*/

bool skHitranMoleculeManager::FetchAllIsotopeEntries( size_t moleculeid, std::vector<const skHitranPartitionTableEntry*> *table ) const
{
	const_iterator						iter;
	const skHitranPartitionTableEntry*	entry;

	table->clear();
	table->reserve(40);
	for (iter = m_molecules.begin(); !(iter == m_molecules.end()); ++iter)
	{
		entry = &(*iter).second;
		if ( entry->m_moleculeid == moleculeid)
		{
			
			table->push_back( entry );
		}
	}
	std::sort( table->begin(), table->end(), IsotopeOrderBinaryPred );
	return true;
}

/*-----------------------------------------------------------------------------
 *					skHitranMoleculeManager::FindMoleculeId		2013-3-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool skHitranMoleculeManager::FindMoleculeId( const char* chemicalname, size_t*	moleculeid) const
{
	const_iterator	iter;
	bool		ok;
	nxString	buffer(chemicalname);

	buffer.MakeUpper();
	iter = std::find_if( m_molecules.begin(), m_molecules.end(), skHitranMoleculeManagerNamePred( buffer, 0) );
	ok = !(iter == m_molecules.end());
	if (!ok)
	{
		*moleculeid = 0;
	}
	else
	{
		*moleculeid = (*iter).second.m_moleculeid;
	}
	return ok;
}





static std::mutex					g_skHitranMoleculeManager_refcountlock;				// then add a mutex for AddRef and Release locking purposes.
static skHitranMoleculeManager*		g_manager2008 = NULL;
static skHitranMoleculeManager*		g_managerhapi = NULL;
/*-----------------------------------------------------------------------------
 *					skHitranMoleculeManager::CreateManagerInstance		2013-3-13*/
/** A static member used to manager **/
/*---------------------------------------------------------------------------*/

const skHitranMoleculeManager* skHitranMoleculeManager::CreateManagerInstance(bool hapicompliant)
{
	skHitranMoleculeManager* result = nullptr;

	std::lock_guard<std::mutex>	lock( g_skHitranMoleculeManager_refcountlock );		// lock the mutex. This needs to be thrrad safe.
	{
		if (hapicompliant)
		{
			if (g_managerhapi == NULL) g_managerhapi = new skHitranMoleculeManager(true);	// Allocate the single instance
			g_managerhapi->AddRef();														// And add a reference
			result = g_managerhapi;

		}
		else
		{
			if (g_manager2008 == NULL) g_manager2008 = new skHitranMoleculeManager(false);							// Allocate the single instance
			g_manager2008->AddRef();																	// And add a reference
			result = g_manager2008;
		}
	}																							// and that is that
	return result;																			// return the manager instance
}

/*-----------------------------------------------------------------------------
 *					skHitranMoleculeManager::skHitranMoleculeManager		2013-3-13*/
/** **/
/*---------------------------------------------------------------------------*/

skHitranMoleculeManager::skHitranMoleculeManager( bool hapicompliant)
{
	bool	ok;

	m_hapicompliant = hapicompliant;
	ok =       LoadMoleculeDefinitions();													//
	ok = ok && LoadPartitionDefinitions();
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skHitranMoleculeManager, Error initializing the Hitran molecule database");
	}

}


/*-----------------------------------------------------------------------------
 *					skHitranMoleculeManager::~skHitranMoleculeManager		2013-3-13*/
/** **/
/*---------------------------------------------------------------------------*/

skHitranMoleculeManager::~skHitranMoleculeManager()
{
	m_molecules.clear();
	{
		std::lock_guard<std::mutex>	lock( g_skHitranMoleculeManager_refcountlock );		// lock the mutex. This needs to be therad safe.
		if (m_hapicompliant)
		{
			if (g_managerhapi == this) g_managerhapi = nullptr;
		}
		else
		{
			if (g_manager2008 == this) g_manager2008 = nullptr;
		}
	}
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_IceCrystalCached::LoadParametersFromRegistrySetting		2009-6-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skHitranMoleculeManager::LoadBaseDirectoryNameFromRegistry( nxString* filename )
{
	nxRegistryConfiguration		config( "USask-ARG", "skOpticalProperties/Hitran/",nxRegistryConfiguration::GLOBAL_INI, true);
	bool						ok;

	ok = config.LocateDirectoryFromKey( "BaseDirectory", filename, true, true, "Enter location of HITRAN base directory");
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skSpectralLineCollection_HitranChemical::LoadDirectoryNameFromRegistry, error loading directory name from registry");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranChemical::FindHitranMoleculeDirectory		2013-3-12*/
/** **/
/*---------------------------------------------------------------------------*/

bool skHitranMoleculeManager::FindHitranMoleculeDirectory( nxString* moleculedir) const
{
	nxString	filespec;
	nxString	basedir;
	bool		ok;
	
	ok = LoadBaseDirectoryNameFromRegistry( &basedir );
	if (ok)
	{
		if (m_hapicompliant)
		{
			moleculedir->sprintf("%s/By-Molecule/Uncompressed-files", (const char*)basedir);
			moleculedir->MakeDirectorySeparatorsOSConsistent();
			ok = nxDirectory::FileExists(*moleculedir);
		}
		else
		{
			basedir.EnsureLastCharIsDirectoryChar();
			filespec.sprintf("HITRAN2*");

			nxDirectory	files;
			files.ScanDirectory(filespec, true, basedir);

			ok = (files.List().GetSize() == 1);
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "skHitranMoleculeManager::FindHitranMoleculeDirectory, There was not exactly one entry that matched %s", (const char*)filespec);
			}
			else
			{
				moleculedir->sprintf("%s/By-Molecule/Uncompressed-files", (const char*)files.List().GetAt(0));
				moleculedir->MakeDirectorySeparatorsOSConsistent();
				ok = nxDirectory::FileExists(*moleculedir);
			}
		}
		if (!ok)
		{
				nxLog::Record(NXLOG_WARNING,"skHitranMoleculeManager::FindHitranMoleculeDirectory, Directory %s does not exist. Thats not good", (const char*)(*moleculedir) );
		}
		moleculedir->EnsureLastCharIsDirectoryChar();
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					skHitranMoleculeManager::FindHitranGlobalFileFile		2013-3-12*/
/** **/
/*---------------------------------------------------------------------------*/

bool skHitranMoleculeManager::FindHitranGlobalFileFile( const char* filename, nxString* fullfilename )
{
	bool		ok;
	nxString	basedir;

	ok = LoadBaseDirectoryNameFromRegistry( &basedir );
	if (ok)
	{
		basedir.EnsureLastCharIsDirectoryChar();
		fullfilename->sprintf("%sGlobal_Data/%s",(const char*)basedir, (const char*) filename);
		fullfilename->MakeDirectorySeparatorsOSConsistent();
		ok = nxDirectory::FileExists( *fullfilename);
		if (!ok)
		{
			fullfilename->sprintf("%sGlobal-Data/%s",(const char*)basedir, (const char*) filename);
			fullfilename->MakeDirectorySeparatorsOSConsistent();
			ok = nxDirectory::FileExists( *fullfilename);
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"skHitranMoleculeManager::FindHitranMolparamFile, File %s does not exist. Thats not good", (const char*)(filename) );
			}
		}
	}
	if (!ok) fullfilename->Empty(false);
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skHitranMoleculeManager::FindHitranMolparamFile		2013-3-12*/
/** **/
/*---------------------------------------------------------------------------*/

bool skHitranMoleculeManager::FindHitranMolparamFile( nxString* filename)
{
	bool ok;
	ok = FindHitranGlobalFileFile( "molparam.txt", filename);
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skHitranMoleculeManager::FindHitranParsumFile		2013-3-12*/
/** **/
/*---------------------------------------------------------------------------*/

bool skHitranMoleculeManager::FindHitranParsumFile( nxString* filename)
{
	bool	ok;

	ok = FindHitranGlobalFileFile( "parsum.dat", filename);
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skHitranMoleculeManager::LoadMoleculeDefinitions		2013-3-12*/
/** **/
/*---------------------------------------------------------------------------*/

bool skHitranMoleculeManager::LoadMoleculeDefinitions( )
{
	bool							ok;
	bool							ok1;
	size_t							numgood = 0;
	size_t							numbad  = 0;
	char							line[1024];
	std::ifstream					inputfile;
	nxString						filename;
	nxString						name;
	nxStringArray					tokens;
	skHitranPartitionTableEntry		entry;
	size_t							isotopeorder = 0;

	m_molecules.clear();

	ok = FindHitranMolparamFile(&filename);
	if (ok)
	{
		inputfile.open(	filename);
		ok = !inputfile.fail();
		if (ok)
		{
			inputfile.getline( line, N_ELEMENTS(line), '\n' );			// Skip the first line

			while (!inputfile.eof())
			{
				inputfile.getline( line, N_ELEMENTS(line), '\n' );			// Skip the first line
				tokens.Strtok( line, " ()\r");
				if (tokens.GetSize() == 2)
				{
					name                 = tokens.GetAt(0);
					name.MakeUpper();
					entry.m_chemicalname = (const char*)name;
					entry.m_moleculeid   = atoi(tokens.GetAt(1));
					entry.m_isotopeid    = 0;
					entry.m_isotopeorder = 0;
					entry.m_abundance    = 0.0;
					entry.m_gj           = 0.0;
					entry.m_molarmass    = 0.0;		// The molar mass as givn in HITRAN file molparam.txt
					entry.m_hapicompliant = m_hapicompliant;
					ok1 = (entry.m_moleculeid > 0) && (entry.m_moleculeid < 1000);
					isotopeorder         = 0;
				}
				else if (tokens.GetSize() == 5)
				{
					entry.m_isotopeorder = ++isotopeorder;						// The isotope order follows the order as stored in the molparam.txt file
					entry.m_isotopeid = atoi( tokens.GetAt(0) );
					entry.m_abundance = atof( tokens.GetAt(1) );
					entry.m_gj        = atof( tokens.GetAt(3) );
					entry.m_molarmass = atof( tokens.GetAt(4) );
					entry.m_hapicompliant = m_hapicompliant;

					ok1 =      (entry.m_moleculeid >  0  ) && (entry.m_moleculeid < 1000)
							&& (entry.m_abundance  >  0.0) && (entry.m_abundance  <= 1.0)
							&& (entry.m_isotopeid  >  0  ) && (entry.m_isotopeid  < 999999)
							&& (entry.m_molarmass  >  0.0) && (entry.m_molarmass  < 1000.0);
					if (!ok1)
					{
						numbad++;
					}
					else
					{
						numgood++;
						value_type	value( entry.ToKey(), entry );
						m_molecules.insert( value );
					}
				}
				else
				{
					if (tokens.GetSize() > 0)
					{
						numbad++;
					}
				}
			}
		}
		inputfile.close();
		ok = ok && (numbad == 0);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skHitranMoleculeManager::LoadMoleculeDefinitions, There were (%d) errors loading in values form hitran file molparam.txt", (int)numbad);
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skHitranMoleculeManager::LoadMoleculeDefinitions		2013-3-12*/
/** Reads in a decodes the lines of parsum.dat. An brief exmaple of the file
 *	format (as of 2013-03-13) is given below.
 *
 * Temp(K)            H2O_161                    H2O_181                    H2O_171
 *  70.0           21.000000                  20.171000                 120.535679 
 *
**/
/*---------------------------------------------------------------------------*/

bool skHitranMoleculeManager::LoadPartitionDefinitions( )
{
	bool										ok = true;
	bool										ok1;
	bool										ok2;
	char										line[4096];
	std::ifstream								inputfile;
	nxString									filename;
	nxString									name;
	nxStringArray								tokens;
	size_t										i;
	std::vector<skHitranPartitionTableEntry*>	partitionmolecules;
	skHitranPartitionTableEntry*				entry;
	const char*									chemicalname;
	size_t										isotopeid;
	size_t										numtokens;
	size_t										nummolecules;
	double										temperature;
	double										q;

	if (!m_hapicompliant)
	{
		ok = FindHitranParsumFile(&filename);
		if (ok)
		{
			inputfile.open(filename);
			ok = !inputfile.fail();
			if (ok)
			{
				inputfile.getline(line, N_ELEMENTS(line), '\n');			// Skip the first line
				tokens.Strtok(line, "_ \r");
				numtokens = tokens.GetSize();
				ok = (numtokens > 2) && ((numtokens & 1) == 1);				// Make sure we have at least one moluecule and an even number of tokens
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING, "skHitranMoleculeManager::LoadPartitionDefinitions, The first line of the parsum.dat is not as expected");
				}
				else
				{
					nummolecules = (numtokens - 1) / 2;
					partitionmolecules.reserve(nummolecules);
					for (i = 0; i < nummolecules; i++)
					{
						chemicalname = (const char*)(tokens.GetAt(2 * (int)i + 1));
						isotopeid = atoi(tokens.GetAt(2 * (int)i + 2));
						ok2 = FindMoleculeEntry(chemicalname, isotopeid, &entry);
						partitionmolecules.push_back(entry);
						if (!ok2)
						{
							nxLog::Verbose(NXLOG_WARNING, "skHitranMoleculeManager::LoadPartitionDefinitions, Could not find an entry in molparam.txt for checmical %s isotope %d", (const char*)chemicalname, isotopeid);
						}
					}

					for (i = 0; i < nummolecules; i++)
					{
						if (partitionmolecules.at(i) != NULL)
						{
							partitionmolecules.at(i)->ReserveSpace(3000);
						}
					}
					ok1 = true;
					while (!inputfile.eof())
					{
						inputfile.getline(line, N_ELEMENTS(line), '\n');			// Skip the first line
						tokens.Strtok(line, " ");
						ok2 = tokens.GetSize() == (nummolecules + 1);
						if (ok2)
						{
							temperature = atof(tokens.GetAt(0));
							for (i = 0; i < nummolecules; i++)
							{
								if (partitionmolecules.at(i) != NULL)
								{
									q = atof(tokens.GetAt((int)i + 1));
									partitionmolecules.at(i)->AddEntry(temperature, q);
								}
							}
						}
						else
						{
							ok2 = tokens.GetSize() == 0;
						}
						ok1 = ok1 && ok2;
					}
					if (!ok1)
					{
						nxLog::Record(NXLOG_WARNING, "skHitranMoleculeManager::LoadPartitionDefinitions, Error loading in the internal partition data from parsum.data. Thats not good");
					}
					ok = ok && ok1;
				}
			}
			inputfile.close();
		}
	}
	return ok;
}

