# Copy the fortran DLLs to the local .exe location for testing/debugging
if(WIN32)

	file(GLOB FORTRAN_DLLS "lib/Windows_x64/*.dll")
	file(COPY ${FORTRAN_DLLS} DESTINATION ../Repos_SasktranV3)

endif()