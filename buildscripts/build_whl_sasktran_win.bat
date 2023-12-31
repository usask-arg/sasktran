ECHO OFF
SET PYENV=%1
SET SASKTRAN_IF_INCLUDE_PATH=%2
SET SASKTRAN_IF_LIBRARY_PATH=%3

ECHO %3

call conda activate %PYENV%
pushd pysasktran 
python -m pip install setuptools-scm
python -m pip wheel . -w ..\wheelhouse --no-deps --no-build-isolation
IF %ERRORLEVEL% NEQ 0 EXIT /B 1
popd
call conda deactivate
