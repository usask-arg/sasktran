#pragma once
/*---------------------------------------------------------------------------
 *                    Class HitranIsotopeCache                    2019-11-07 */
/** **/
/*---------------------------------------------------------------------------*/

class HitranIsotopeCache
{
	private:
		int												m_isotopeid;
		std::vector<HitranLineStruct>					m_lines;			// Cache of lines in the file, one ste of lines for each isotope
		std::vector<double>								m_nu00;				// The wavenumber of all the lines in the file

	public:
														HitranIsotopeCache( int isotopeid);
		bool											LoadSpectralLines	(FILE* f, double min_wavenum, double max_wavenum );
		static int										RecordSize			();
		const std::vector<HitranLineStruct>&			Lines				() const { return m_lines;}

};

/*---------------------------------------------------------------------------
 *                    Class HitranLineStructCache                     2019-11-06 */
/** A class that caches the Hitran spectral lines for each molecule in a binary
 *	file. The binary file has the format::
 *
 *				int 								versionid				// Version ID of this file
 *			    int 								recsize					// The binary size fo each record. This may be very compielr/machine specific
 *				int									numisotopes				// The number of isotopes for this molecule
 *				for each isotope
 *				{
 *				    int								 isotopeid				// The isotopeid
 *					int								 syncinteger			// Equals 0xDEADBEEF = 3735928559
 *				    int								 numlines				// The total number of lines
 *					double							 nu00[numlines]			// The wavenumber of each line in the database. This is in ascending wavenumber order
 *					skSpectralLine_HitranLineStruct	 lines[numlines]		// The spectral line info for each line in teh database. It is in the same order as nu00
 *				}
 *	The class will only load the spectral lines contained between user supplied min and max wavenumbers
 *	which can significantly increase file I/O performance.
 **/
/*---------------------------------------------------------------------------*/

class HitranLineStructCache
{
	private:
		int														m_versionid;
				std::map< int, HitranIsotopeCache >				m_isotopes;			// Cache of lines in the file, one ste of lines for each isotope
		typedef std::map< int, HitranIsotopeCache >::value_type	value_type;
		typedef std::map< int, HitranIsotopeCache >::iterator	iterator;

	private:
		bool											LoadBaseDirectoryNameFromRegistry	( nxString* filename );
		bool											FindFile							( int molNum, std::string* filename, bool* exists);
		bool											IsValidIsotopeID					( int isotopeid, const std::vector<const skHitranPartitionTableEntry*>&	isotopetable);
	public:
														HitranLineStructCache				();
		bool											LoadSpectralLines					(skSpectralLineCollection_HitranChemical* chemical,  const std::vector<const skHitranPartitionTableEntry*>&	isotopetable, double min_wavenum, double max_wavenum );
		bool											WriteSpectralLines					(skSpectralLineCollection_HitranChemical* chemical );
};


