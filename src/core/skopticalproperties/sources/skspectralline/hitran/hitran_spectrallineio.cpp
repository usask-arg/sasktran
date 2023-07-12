#include <skopticalproperties21.h>

/*---------------------------------------------------------------------------
 *             HitranIsotopeCache::HitranIsotopeCache             2019-11-07 */
/** **/
/*---------------------------------------------------------------------------*/

HitranIsotopeCache::HitranIsotopeCache( int isotopeid)
{
	m_isotopeid = isotopeid;
}

/*---------------------------------------------------------------------------
 *                           RecordSize                           2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/

int HitranIsotopeCache::RecordSize()
{
	return (int)sizeof(HitranLineStruct );
}


/*---------------------------------------------------------------------------
 *              HitranLineStructCache::LoadSpectralLines              2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/

bool HitranIsotopeCache::LoadSpectralLines (FILE* f, double min_wavenum, double max_wavenum )
{
	bool					ok= true;
	std::string				filename;
	std::vector<double>		allnu;
	int						maxrecords;
	long					nextpos;
	size_t					idx0;
	size_t					idx1;
	size_t					numrecords;
	size_t					offset;

	ok = ok && (fread( &maxrecords,  sizeof(int), 1, f) == 1);			// Read in the maxrecords.
	ok = ok && (maxrecords > 0) && (maxrecords < 100000000);
	//printf("lines = %d lines\n", (int)maxrecords);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"HitranIsotopeCache::LoadSpectralLines, binary file cache appears to be corrupted as key parameters are invalid,  maxrecords = %d", (int) maxrecords);
	}
	else
	{
		allnu.resize( maxrecords );
		ok = ok &&  (fread( allnu.data(), sizeof(double), maxrecords, f) == maxrecords);		// Load in the wavenumber of all the lines
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING," HitranIsotopeCache::LoadSpectralLines, There were errors loading in the wave numbers of the spectral lines from file [%s]",(const char*)filename.c_str());
		}
		else
		{
			auto iter0 = std::lower_bound( allnu.begin(), allnu.end(), min_wavenum );			// Find the start and end of tthe
			auto iter1 = std::upper_bound( allnu.begin(), allnu.end(), max_wavenum );			// required spectral lines. They are in sorted  order

			//if ( (iter1 != allnu.begin()) && (iter0 != iter1 )) --iter1;						// step back one to get bounding lines and make sure we dont go below zero
			idx0 = iter0 - allnu.begin();														// Get the index of the first
			idx1 = iter1 - allnu.begin();														// and last record													
			numrecords = idx1 - idx0;															// get the number of subset records required to be read from the file

			m_lines.resize( numrecords );																// allocate contiguous space for the records
			nextpos = ftell( f) + maxrecords*RecordSize();											// Get the current file position and work out the start of the next isotope 
			ok = (numrecords == 0);																		// If we have zero records then we are done
			if (!ok)																					// otherwise we have records to read in
			{																							// so
				offset  = idx0*RecordSize();															// Get the offset in bytes from the current file position to the first record
				ok      =       (fseek( f, (long)offset, SEEK_CUR ) == 0);								// Go to the start of the first record
				ok      = ok && (fread( &m_lines.front(), RecordSize(), numrecords, f) == numrecords);	// suck in all the required spectral lines in one gulp
				m_nu00.assign( allnu.begin()+idx0, allnu.begin()+idx1 );								// and assign the subset of wavenumbers
			}
			else
			{
				m_nu00.clear();
			}
			//long thispos = ftell( f);
			fseek( f, (long)nextpos, SEEK_SET );													// Position the file to the start of the next  isotope (this may hit EOF)
			//long apos = ftell(f);
			//printf("Seeking = %d (%d)\n", (int)nextpos, (int)apos);
		}
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *              HitranLineStructCache::HitranLineStructCache              2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/

HitranLineStructCache::HitranLineStructCache()
{
	m_versionid = 0xDEADBEEF;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_IceCrystalCached::LoadParametersFromRegistrySetting		2009-6-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool HitranLineStructCache::LoadBaseDirectoryNameFromRegistry( nxString* filename )
{
	nxRegistryConfiguration		config( "USask-ARG", "skOpticalProperties/Hitran/",nxRegistryConfiguration::GLOBAL_INI, true);
	bool						ok;

	ok = config.LocateDirectoryFromKey( "SpectralLineCacheDir", filename, true, true, "Enter location of HITRAN spectral line cache directory");
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"HitranLineStructCache::LoadBaseDirectoryNameFromRegistry, error loading directory name from registry");
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *                  HitranLineStructCache::FindFile                   2019-11-07 */
/** **/
/*---------------------------------------------------------------------------*/

bool HitranLineStructCache::FindFile( int molNum, std::string* filename, bool* exists)
{
	nxString basedir;
	nxString	name;
	bool ok;

	ok = LoadBaseDirectoryNameFromRegistry( &basedir );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"HitranLineStructCache::FindFile, error fetching base directory for the spectral line cache for molecule id (%d)", (int)molNum);
		filename->clear();
		*exists = false;
	}
	else
	{
		basedir.EnsureLastCharIsDirectoryChar();
		if (!nxDirectory::FileExists(basedir) && *exists)
		{
			nxDirectory::CreateADirectory(basedir);
		}
		name.sprintf("%slinecache2_mol%03d.bin",(const char*)basedir, (int)molNum);
		name.MakeDirectorySeparatorsOSConsistent();
		*filename = (const char*)name;
		*exists   = nxDirectory::FileExists(name);
	
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *            HitranLineStructCache::IsValidIsotopeID             2019-11-07 */
/** **/
/*---------------------------------------------------------------------------*/

bool HitranLineStructCache::IsValidIsotopeID( int isotopeid, const std::vector<const skHitranPartitionTableEntry*>&	isotopetable)
{
	bool	ok = false;
	auto iter = isotopetable.begin();

	while (!ok && !(iter == isotopetable.end()))
	{
		ok = (*iter++)->m_isotopeid == isotopeid;
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *              HitranLineStructCache::LoadSpectralLines              2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/

bool HitranLineStructCache::LoadSpectralLines (skSpectralLineCollection_HitranChemical* chemical,  const std::vector<const skHitranPartitionTableEntry*>&	isotopetable, double min_wavenum, double max_wavenum )
{
	bool					ok, ok1;
	std::string				filename;
	bool					exists=false;
	FILE*					f;
	int						recsize;
	int						numisotopes;
	int						versionid;
	int						syncinteger;
	int						isotopeid;
	int						molnum;
	HitranIsotopeCache*		isotopecache;

	molnum = chemical->MoleculeNumber();
	m_isotopes.clear();

	ok = FindFile(molnum, &filename, &exists);
	ok = ok && exists;
	if (ok)
	{
		f = fopen( filename.c_str(), "rb");
		ok = (f != nullptr);
		if (ok)
		{
			ok = ok && (fread( &versionid,   sizeof(int), 1, f) == 1);			// Read in the file version
			ok = ok && (fread( &recsize,     sizeof(int), 1, f) == 1);			// Read in the record size
			ok = ok && (fread( &numisotopes, sizeof(int), 1, f) == 1);			// Read in the maxrecords.
			ok = ok && (versionid == m_versionid);								// Check that the versionid is as expected.
			ok = ok && (recsize == HitranIsotopeCache::RecordSize());			// Check that record size is as expected.
			ok = ok && (numisotopes > 0) && (numisotopes < 100);				// Check to make sure the number of isotopes are not crazy.
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"HitranLineStructCache::LoadSpectralLines, Error loading spectral line cache header. You should delete file [%s] so it can be regenerated", (const char*)filename.c_str());
			}
			else
			{
				for ( int i = 0;  i < numisotopes; i++ )
				{
					long pos = ftell(f);
					ok1 =        (fread( &isotopeid,   sizeof(int), 1, f) == 1);				// Read in the isotopeid
					ok1 = ok1 && (fread( &syncinteger, sizeof(int), 1, f) == 1);				// Read in the integer sync byte
					ok1 = ok1 && (syncinteger == m_versionid );									// Make sure the sync value is as expected.
					//printf("Reading isotopeid = %d sync integer = %08x, pos = %d ", (int)isotopeid, (int)syncinteger, (int)pos);

					if (!ok1)
					{
						nxLog::Record(NXLOG_WARNING," HitranLineStructCache::LoadSpectralLine,Errors keeping synchronization with the sync word, expected %08x but got %08x. Thats a problem and needs debugging", (int)m_versionid, (int)syncinteger);

					}
					else
					{
						ok1 = ok1 && IsValidIsotopeID( isotopeid, isotopetable );
						ok1 = ok1 && m_isotopes.find( isotopeid) == m_isotopes.end();				// Make sure we dont already have entries for this isotope
						if (!ok1)
						{
							nxLog::Record(NXLOG_WARNING," HitranLineStructCache::LoadSpectralLine, Errors with the isotopeid read from file (%d). It is either diplicated or is not a valid hitran isotopeid.. Both cases are problems and need debugging", (int) isotopeid);
						}
						else
						{
							auto status = m_isotopes.insert( value_type( isotopeid, HitranIsotopeCache(isotopeid) ) );
							isotopecache =  &(status.first->second);
							ok1 = status.second;
							ok1 = ok1 && isotopecache->LoadSpectralLines( f, min_wavenum, max_wavenum );
							ok1 = ok1 && chemical->InsertAllSpectralLines( isotopeid, isotopecache->Lines() );
							ok = ok && ok1;
						}

					}
				}
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING," HitranLineStructCache::LoadSpectralLine, There were errors loading in the spectral lines for one or more isotopes");
				}
			}
			fclose(f);
		}
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *              HitranLineStructCache::LoadSpectralLines              2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/

bool HitranLineStructCache::WriteSpectralLines (skSpectralLineCollection_HitranChemical* chemical )
{
	bool											ok, ok1;
	std::string										filename;
	bool											exists = true;
	FILE*											f;
	int												recsize;
	int												numisotopes;
	int												versionid   = m_versionid;
	int												syncinteger = m_versionid;
	int												isotopeid;
	size_t											numlines;
	std::vector<double>								nu00;
	std::vector<HitranLineStruct>					lineentries;
	double											lastnu;
	bool											isascending;
	const skSpectralLineCollection_HitranIsotope*	entry;
	const std::vector< skSpectralLineEntry*  >*		spectral_lines;
	const skSpectralLine_HitranLine*			    hitran_line;
	int												molnum;

	molnum = chemical->MoleculeNumber();
	m_isotopes.clear();

	ok = FindFile(molnum, &filename, &exists);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"HitranLineStructCache::WriteSpectralLines, Error fetching filename for molecule [%d].",(int)molnum);
	}
	else
	{
		f = fopen( filename.c_str(), "wb");
		ok = (f != nullptr);
		if (ok)
		{
			numisotopes = (int) chemical->Isotopes().size();
			recsize     = (int) HitranIsotopeCache::RecordSize();

			ok = ok && (fwrite( &m_versionid,  sizeof(int), 1, f) == 1);			// Read in the file version
			ok = ok && (fwrite( &recsize,     sizeof(int), 1, f) == 1);			// Read in the record size
			ok = ok && (fwrite( &numisotopes, sizeof(int), 1, f) == 1);			// Read in the maxrecords.
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"HitranLineStructCache::WriteSpectralLines, Error writing spectral line cache header to file [%s]. It might be locked by another application", (const char*)filename.c_str());
			}
			else
			{
				for ( auto iter = chemical->Isotopes().begin(); !(iter == chemical->Isotopes().end()); ++iter)
				{
					long pos = ftell(f);
					isotopeid = (int)(iter->first);
					entry     = &(iter->second);
					ok1 =        (fwrite( &isotopeid,   sizeof(int), 1, f) == 1);				// Read in the isotopeid
					ok1 = ok1 && (fwrite( &syncinteger, sizeof(int), 1, f) == 1);				// Read in the integer sync byte
					spectral_lines = &entry->SpectralLines();									// Fetch the spectral lines read in by the HITRAN code

					numlines = spectral_lines->size();											// Get the number of lines for this isotope
					nu00.resize(numlines);														// resize the wavenumber array as we prepare to copy
					lineentries.resize(numlines);												// resize the list of line entries as we prepare to copy
					isascending = true;															// reset flag to check that wavenumbers are ascending
					lastnu      = -9999.0;
					//printf("Writing isotopeid = %d sync integer = %08x, pos = %d, numlines = %d\n", (int)isotopeid, (int)syncinteger, (int)pos, (int)numlines);

					for (size_t i = 0; i < numlines; i++)
					{
						hitran_line       = dynamic_cast<const skSpectralLine_HitranLine*>(spectral_lines->at(i)->SpectralLine());
						lineentries.at(i) = hitran_line->HitranEntry();
						nu00.at(i)        = spectral_lines->at(i)->Nu00();
						isascending       = (nu00.at(i) >= lastnu);
						lastnu            = nu00.at(i);
					}
					ok = isascending;
					if (!ok)
					{
						nxLog::Record(NXLOG_WARNING,"HitranLineStructCache::WriteSpectralLines, The wavenumbers for isotope(%d) are not in ascending order. Thats a problem that needs debugging and fixing",(int) isotopeid);
					}
					else
					{
						int maxrecords = (int)numlines;
						ok = ok && (fwrite( &maxrecords,         sizeof(int), 1, f) == 1);												// Write the number of records.
						ok = ok && (fwrite( nu00.data(),         sizeof(double), maxrecords, f) == maxrecords);							// Write the wavenumber of all the lines
						ok = ok && (fwrite( lineentries.data(),  HitranIsotopeCache::RecordSize(), maxrecords, f) == maxrecords);		// Write the HITRAN entry of all the lines
						if (!ok)
						{
							nxLog::Record(NXLOG_WARNING,"HitranLineStructCache::WriteSpectralLines, Error writing %d records for isotope(%d).Thats a problem that needs debugging and fixing", (int)maxrecords, (int) isotopeid);
						}
					}
				}
			}
			fclose(f);

		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING," HitranIsotopeCache::WriteSpectralLines, There were errors writing hitran spectral spectral lines to cache file file [%s]",(const char*)filename.c_str());
	}
	return ok;
}



