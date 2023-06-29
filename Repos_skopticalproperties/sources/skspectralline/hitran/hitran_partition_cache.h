

class skHitranPartitionTableEntry;

/*---------------------------------------------------------------------------
 *                 Class HitranPartitioTableCache                 2019-11-08 */
/** **/
/*---------------------------------------------------------------------------*/

class HitranPartitionTableCache
{
	private:
		const skHitranPartitionTableEntry*	m_parent;
		int									m_versionid;
		double								m_deltaT;
		double								m_Tstart;
		double								m_Tend;
		nx1dArray<double>					m_Q;

	private:
		void							CreateTable();
		double							InterpolateTable(double T);
		bool							LoadBaseDirectoryNameFromRegistry( nxString* filename );
		bool							FindFile( std::string* filename, bool* exists);
		bool							LoadCache();
		bool							CreateCache();
		bool							CheckAndLoadCache();

	public:
										HitranPartitionTableCache( const skHitranPartitionTableEntry*	parent );
		double							InternalPartition( double T );
};

