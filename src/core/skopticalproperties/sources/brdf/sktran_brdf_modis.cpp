#include <skopticalproperties21.h>
#include <sources/sasktranif_opticalimpl/skbrdf_stubs.h>



/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_MODIS_RossThickLiSparseReciprocal::ISKBrdf_Stub_MODIS_RossThickLiSparseReciprocal		 2017-07-28*/
/** **/
/*---------------------------------------------------------------------------*/

ISKBrdf_Stub_MODIS_RossThickLiSparseReciprocal::ISKBrdf_Stub_MODIS_RossThickLiSparseReciprocal(SKTRAN_BRDF_MODIS_RossThickLiSparseReciprocal* modisbrdf)
	: ISKBrdf_Stub_Base(modisbrdf)
{
	m_modisbrdf = modisbrdf;
}


/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_MODIS_RossThickLiSparseReciprocal::~ISKBrdf_Stub_MODIS_RossThickLiSparseReciprocal		 2017-07-28*/
/** **/
/*---------------------------------------------------------------------------*/

ISKBrdf_Stub_MODIS_RossThickLiSparseReciprocal::~ISKBrdf_Stub_MODIS_RossThickLiSparseReciprocal()
{}

/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_MODIS_RossThickLiSparseReciprocal::SetPropertyArray		 2017-07-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKBrdf_Stub_MODIS_RossThickLiSparseReciprocal::SetPropertyArray(const char* propertyname, const double* value, int numpoints)
{
	nxString					name(propertyname);
	bool						ok = m_modisbrdf != nullptr;

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "ISKBrdf_Stub_MODIS_RossThickLiSparseReciprocal::SetProperty, the internal C++ object is undefined. Thats not good");
	}
	else
	{
		if (name == "BRDFParameters")
		{
			ok = (numpoints == 3);
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "ISKBrdf_Stub_MODIS_RossThickLiSparseReciprocal::SetProperty(BRDFParameters), accepts only 3 parameters. You passed in %d parameters", (int)numpoints);
			}
			else
			{
				double f_iso = value[0]; // isotropic (lambertian) component
				double f_vol = value[1]; // linear weight of volume scattering kernel (Ross-Thick)
				double f_geo = value[2]; // linear weight of geometric scattering kernel (Li-Sparse)


				ok = m_modisbrdf->SetBRDFParameters(f_iso, f_vol, f_geo);
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING, "ISKBrdf_Stub_MODIS_RossThickLiSparseReciprocal::SetProperty(BRDFParameters), failed to set parameters to f_iso=%e, f_vol=%e, f_geo=%e", (double)f_iso, (double)f_vol, (double)f_geo);
				}
			}
		}
		else
		{
			ok = ISKBrdf_Stub_Base::SetPropertyArray(propertyname, value, numpoints);
		}
	}
	return ok;
}



/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_MODIS_RossThickLiSparseReciprocal::SKTRAN_BRDF_MODIS_RossThickLiSparseReciprocal		 2017-07-28*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_BRDF_MODIS_RossThickLiSparseReciprocal::SKTRAN_BRDF_MODIS_RossThickLiSparseReciprocal()
{
	m_f.emplace_back(std::numeric_limits<double>::quiet_NaN());
	m_f.emplace_back(std::numeric_limits<double>::quiet_NaN());
	m_f.emplace_back(std::numeric_limits<double>::quiet_NaN());
	
	SKTRAN_BRDF_Lambertian* k_iso = new SKTRAN_BRDF_Lambertian(1.0);
	SKTRAN_BRDF_RossThick_Kernel* k_vol = new SKTRAN_BRDF_RossThick_Kernel();
	SKTRAN_BRDF_LiSparseReciprocal_Kernel* k_geo = new SKTRAN_BRDF_LiSparseReciprocal_Kernel();
	k_geo->SetBRDFParameters(1.0, 2.0);

	m_kernels.emplace_back(k_iso);
	m_kernels.emplace_back(k_vol);
	m_kernels.emplace_back(k_geo);
}


/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_MODIS_RossThickLiSparseReciprocal::SetBRDFParameters		 2017-07-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_BRDF_MODIS_RossThickLiSparseReciprocal::SetBRDFParameters(double f_iso, double f_vol, double f_geo)
{
	m_f[0] = f_iso;
	m_f[1] = f_vol;
	m_f[2] = f_geo;
	return true;
}
