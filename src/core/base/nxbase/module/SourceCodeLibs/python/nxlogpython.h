//---------------------------------------------------------------------------
//						class  nxLogMATLAB
//---------------------------------------------------------------------------

class  nxLogPython : public nxLog
{
	protected:
		virtual void		DisplayEntry( const nxLogEntry& entry);

	public:
							nxLogPython();
};

