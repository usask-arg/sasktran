
void CompareSingleScatter();

void SK_SingleScatter( const std::vector<double>& wavelen, const SKTRAN_LineOfSightArray_V21& linesofsight, SKTRAN_AtmosphericOpticalState_V21& opticalstate, const nxVector& sun, nx2dArray<double>& radiance );
void HR_SingleScatter( const std::vector<double>& wavelen, const SKTRAN_LineOfSightArray_V21& linesofsight, SKTRAN_AtmosphericOpticalState_V21& opticalstate, const nxVector& sun, nx2dArray<double>& radiance, bool usecurve, std::string filename );
void MC_SingleScatter( const std::vector<double>& wavelen, const SKTRAN_LineOfSightArray_V21& linesofsight, SKTRAN_AtmosphericOpticalState_V21& opticalstate, const nxVector& sun, nx2dArray<double>& radiance, const double& mc_precision );