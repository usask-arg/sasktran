#include "modules/sktran_do_deprecated/include/sktran_do.h"
#if 0
#include <windows.h>
bool sktran_do_detail::LoadDISORTDLL(LPCWSTR path, LPCSTR procname, DISORT_FPTR* fptr) {
	bool ok = true;
	// Attempt to load dll
	HINSTANCE handle = LoadLibrary((LPCSTR)path);
	// Check if it was successfully opened
	if (handle == nullptr) {
		nxLog::Record(NXLOG_ERROR, "sktran_do_detail::LoadDISORTDLL, failed to load dll.");
		ok = false;
	}
	// Get the function and check we were successful
	*fptr = (DISORT_FPTR) GetProcAddress(handle, procname);
	if (fptr == nullptr) {
		nxLog::Record(NXLOG_ERROR, "SKTRAN_DO_UserSpec::LoadDISORT, dll was loaded but the DISORT function could not be resolved.");
		ok = false;
	}
	return ok;
}
#endif
