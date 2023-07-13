#include "include/sktran_hr_internals.h"

void SKTRAN_HR_Specs_User_Diffuse::ConfigureDefaults()
{
	m_numprofiles = 1;
	m_numofflook = 1;
	m_numbeforehoriz = 6;
	m_numhoriz = 8;
	m_numafterhoriz = 10;
	m_numazi = 12;
	m_numoutgoing = 169;
	m_angleofflook = 2.0;
	m_heightres = 1000.0;
	m_maxdiffuseheight = std::numeric_limits<double>::quiet_NaN();
//	m_surfaceheight = 0.0;
	//m_forcelinearplacement = true;
	m_forcev21incomingsphere = false;
	m_horizonsize = 20.0;
	m_incomingspheretype = SKTRAN_HR_DiffuseIncomingType_Default;
	m_diffuseplacement = SKTRAN_HR_DiffuseProfilePlacement_LinearLOS;
    m_diffuseMatrixStorageMethod = SKTRAN_HR_DiffuseMatrixStorageMethod::scalar;
}

