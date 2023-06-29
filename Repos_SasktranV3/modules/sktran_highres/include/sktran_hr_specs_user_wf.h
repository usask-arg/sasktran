

class SKTRAN_HR_Specs_User_wf
{
	private:
		bool					m_doextinction;			// 
		bool					m_doscatextinction;		// not supported yet, no option to change it implemented

		SKTRAN_HR_wf_Mode		m_wfmode;
		SKTRAN_HR_wf_aerosol_Mode m_aerosolmode;
		double					m_manualwfresolution;
		double					m_maxwfheight;
		double					m_wfinterpwidth;
		std::vector<CLIMATOLOGY_HANDLE>		m_wfspecies;
		std::vector<SKTRAN_HR_wf_Species_Mode> m_wfspeciesmode;
		SKTRAN_HR_wf_precision	m_wfprecision;

		std::vector<double>		m_manualwfheights;
		std::vector<double>		m_manualwfwidth;
		std::vector<double>		m_manualwfwidthleft;
		std::vector<double>		m_manualwfwidthright;
		
		double					m_aerosolsizepercentchange;

	public:
		SKTRAN_HR_Specs_User_wf ();

		void			SetWeightingFunctionMode ( const SKTRAN_HR_wf_Mode& mode ) { m_wfmode = mode; }
		void			SetWeightingFunctionprecision ( const SKTRAN_HR_wf_precision& prec ) { m_wfprecision = prec; }
		void			SetWeightingFunctionSpecies( const std::vector<CLIMATOLOGY_HANDLE>& species ) { m_wfspecies = species; }
		void			SetWeightingFunctionSpeciesString(const std::vector<std::string>& species);
		void			SetWeightingFunctionHeight( const std::vector<double>& height) { m_manualwfheights = height; }
		void			SetWeightingFunctionWidth( const std::vector<double>& width ) { m_manualwfwidth = width; }
		void			SetWeightingFunctionLeftWidth(const std::vector<double>& width) { m_manualwfwidthleft = width; }
		void			SetWeightingFunctionRightWidth(const std::vector<double> & width) { m_manualwfwidthright = width; }
		void			SetAerosolWeightingFunctionMode(const SKTRAN_HR_wf_aerosol_Mode& mode) { m_aerosolmode = mode;  }
		void			SetAerosolSizePercentChange(double change) { m_aerosolsizepercentchange = change; }

		friend class SKTRAN_HR_Specs_Internal_wf;
		friend class SKTRAN_HR_Engine;
};