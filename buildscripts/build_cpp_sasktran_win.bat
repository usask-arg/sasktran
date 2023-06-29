
REM ****************************************************************
REM *	Setup defaults
REM *
REM *	CFG    = the configuration to build (typically Release or Debug)
REM *	PLTFRM = the platform to build
REM ****************************************************************


IF /i %PROCESSOR_ARCHITECTURE%==x86 (SET PLTFRM=win32) ELSE (SET PLTFRM=x64)

REM ****************************************************************
REM *	Process Command line arguments
REM *	%1 = VS2015, VS2017
REM *	%2 = Configuration (eg Debug, Release)
REM *	%3 = Platform ( eg x64, win32 )
REM *   %4 = rebuild (optional) 
REM ****************************************************************

SET SK_SVN_DIR="DONT USE SK_SVN_DIR DIRECTLY"
SET SK_OBJ_DIR="DONT USE SK_OBJ_DIR DIRECTLY"

SET VSXXXX=%1
SET VAR1=%2
SET VAR2=%3
SET VAR4=%4

IF DEFINED VAR2 SET PLTFRM=%VAR2%
IF DEFINED VAR1 SET CFG=%VAR1%

:BUILDVS2012
SET RB=
SET RBLD=%VAR4%
IF A%RBLD%B EQU AB GOTO BUILDVS2012A
ECHO RBLD IS DEFINED = %RBLD%
IF /i %RBLD% EQU /rebuild GOTO BUILDVS2012B
ECHO
ECHO ********* ERROR *************************************************************
ECHO * Invalid 4th Parameter. You can either leave it blank or only use /rebuild *
ECHO *****************************************************************************
ECHO
GOTO FAILEDBUILDERROR

:BUILDVS2012B
SET RB=/t:rebuild

:BUILDVS2012A

MSBUILD %RB% /m /p:Platform=%PLTFRM% /p:Configuration=%CFG%  Repos_BaseCode\nxbase\CompilerIDE\%VSXXXX%\nxbase.sln                   &  IF ERRORLEVEL 1 (GOTO FAILEDBUILDERROR)
MSBUILD %RB% /m /p:Platform=%PLTFRM% /p:Configuration=%CFG%  Repos_BaseCode\nxhdf\CompilerIDE\%VSXXXX%\nxhdfonyx.sln                 &  IF ERRORLEVEL 1 (GOTO FAILEDBUILDERROR)
MSBUILD %RB% /m /p:Platform=%PLTFRM% /p:Configuration=%CFG%  Repos_SasktranIF\CompilerIDE\%VSXXXX%\SasktranIF\sasktranIF.sln         &  IF ERRORLEVEL 1 (GOTO FAILEDBUILDERROR)
MSBUILD %RB% /m /p:Platform=%PLTFRM% /p:Configuration=%CFG%  Repos_skclimatology\CompilerIDE\%VSXXXX%\skclimatology.sln              &  IF ERRORLEVEL 1 (GOTO FAILEDBUILDERROR)
MSBUILD %RB% /m /p:Platform=%PLTFRM% /p:Configuration=%CFG%  Repos_skopticalproperties\CompilerIDE\%VSXXXX%\skopticalproperties.sln  &  IF ERRORLEVEL 1 (GOTO FAILEDBUILDERROR)
MSBUILD %RB% /m /p:Platform=%PLTFRM% /p:Configuration=%CFG%  Repos_SasktranV3\CompilerIDE\%VSXXXX%\sasktranv301.sln                  &  IF ERRORLEVEL 1 (GOTO FAILEDBUILDERROR)
copy locallibs\%VSXXXX%\x64\library\_sasktran_core_internals*.dll Repos_SasktranV3\python_sasktrancore\sasktran_core\*.*

:SUCCESSBUILD
ECHO.
ECHO  	    The script to build the sasktran libraries worked!
EXIT /B 0

:FAILEDBUILDERROR
ECHO.
ECHO ************************** FAILURE ******************************************************
ECHO.
ECHO * The script to build the sasktran libraries failed. Thats a problem!
ECHO *
ECHO ************************** FAILURE *******************************************************
EXIT /B 1
