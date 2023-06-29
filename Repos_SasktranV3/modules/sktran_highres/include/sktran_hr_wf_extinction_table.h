


class SKTRAN_HR_WF_Extinction_Table
{
	private:
		nx2dArray<double>		m_table;	
		std::vector<double>		m_wavenumber;
	public:
		
		double				ExtinctionAt( double wavelen, size_t pertidx );
		void				FillTable( const SKTRAN_HR_WF_Store& pertlist,
									   skOpticalProperties& optprop,
									   const SKTRAN_CoordinateTransform_V2& coords,
									   skClimatology& neutral,
									   const std::vector<double>& wavelen,
									   double mjd );
};