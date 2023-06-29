echo on
SET PYVER=%1
SET ERRCODE=1
call conda env create --file pysasktran\ci\%PYVER%.yml
IF %ERRORLEVEL% NEQ 0 GOTO exitprog
call conda activate %PYVER%_test
IF %ERRORLEVEL% NEQ 0 GOTO exitprog
pip install sasktran -f .\wheelhouse\
IF %ERRORLEVEL% NEQ 0 GOTO exitprog
python -c "import sasktran"
python -m sasktran_core.update_settings set_baum_folder //datastore/valhalla/data/BaumIceCrystals/
python -m sasktran_core.update_settings set_hitran_folder //datastore/valhalla/data/hitranhapi/
python -m sasktran.config.update_settings set_glossac_file \\datastore\valhalla\data\GloSSAC\GloSSAC-V1.nc
IF %ERRORLEVEL% NEQ 0 GOTO exitprog
pytest --pyargs sasktran
SET ERRCODE=%ERRORLEVEL%

:exitprog
call conda deactivate
call conda remove -y -q --name %PYVER%_test --all
echo 'Exiting with error %ERRCODE%'
exit /b %ERRCODE% 
 

