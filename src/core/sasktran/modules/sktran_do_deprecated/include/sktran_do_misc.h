#pragma once

#if 0
#include <Windows.h>

namespace sktran_do_detail
{
	// DISORT function pointer prototype
	typedef void(__stdcall *DISORT_FPTR)(int* NLYR, double* DTAUC, double* SSALB, int* NMOM, double* PMOM, double* TEMPER, double* WVNMLO,
										 double* WVNMHI, int* USRTAU, int* NTAU, double* UTAU, int* NSTR, int* USRANG, int* NUMU,
										 double* UMU, int* NPHI, double* PHI, int* IBCND, double* FBEAM, double* UMU0, double* PHI0,
										 double* FISOT, int* LAMBER, double* ALBEDO, double* BTEMP, double* TTEMP, double* TEMIS,
										 int* PLANK, int* ONLYFL, double* ACCUR, int* PRNT, char* HEADER, int* MAXCLY,
										 int* MAXULV, int* MAXUMU, int* MAXPHI, int* MAXMOM, double* RFLDIR, double* RFLDN,
										 double* FLUP, double* DFDT, double* UAVG, double* UU, double* ALBMED, double* TRNMED);

	// Load a DISORT-2.0 DLL.	** WINDOWS ONLY **
	bool LoadDISORTDLL(LPCWSTR path, LPCSTR procname, DISORT_FPTR* fptr);


	/**
	 * Object which suppresses all console output for the scope it is in.
	 * DISORT has a lot of internal console output so this is used to suppress
	 * that.
	 */
	class SuppressConsoleOutputThisScope
	{
	public:
		SuppressConsoleOutputThisScope()
		{
//			freopen("nul", "w", stdout);
			m_suppressed = true;
		}
		SuppressConsoleOutputThisScope(bool suppress)
		{
			if(suppress) {
			//	freopen("nul", "w", stdout);
				m_suppressed = true;
			}
		}
		~SuppressConsoleOutputThisScope()
		{
		//	if(m_suppressed) freopen("CON", "w", stdout);
		}
	private:
		bool m_suppressed;
	};

}

#endif