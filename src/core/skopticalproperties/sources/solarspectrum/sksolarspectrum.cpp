#include <skopticalproperties21.h>

/*-----------------------------------------------------------------------------
 *					skSolarSpectrum::skSolarSpectrum		 2015- 1- 6*/
/** **/
/*---------------------------------------------------------------------------*/

skSolarSpectrum::skSolarSpectrum()
{
	m_solardistanceinAU = 1.0;
}

/*-----------------------------------------------------------------------------
 *					skSolarSpectrum::~skSolarSpectrum		 2015- 1- 6*/
/** **/
/*---------------------------------------------------------------------------*/

skSolarSpectrum::~skSolarSpectrum()
{
}


/*-----------------------------------------------------------------------------
 *					skSolarSpectrum::SetSolarDistanceFromMjd		 2015- 1- 6*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSolarSpectrum::SetSolarDistanceFromMjd( double mjd )
{
	PlanetSun		sun;
	nxTimeStamp		t(mjd);
	const double	AU = 149597871000.0;

	sun.UpdateECIPosition(t );
	m_solardistanceinAU = sun.Location().Magnitude()/AU;
	return true;
}


/*-----------------------------------------------------------------------------
 *					skSolarSpectrum::SetSolarDistanceFromAU		 2015- 1- 6*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSolarSpectrum::SetSolarDistanceFromAU( double au )
{
	m_solardistanceinAU = au;
	return true;
}


/*-----------------------------------------------------------------------------
 *					skSolarSpectrum::Irradiance		 2015- 1- 6*/
/** **/
/*---------------------------------------------------------------------------*/

double skSolarSpectrum::Irradiance( double wavelen_nm_vacuum )
{
	double irradiance;

	irradiance = IrradianceAt1AU(wavelen_nm_vacuum)/(m_solardistanceinAU*m_solardistanceinAU);
	return irradiance;
}


/*-----------------------------------------------------------------------------
 *					skSolarSpectrum_SAO2010::skSolarSpectrum_SAO2010		 2015- 1- 6*/
/** **/
/*---------------------------------------------------------------------------*/

skSolarSpectrum_TabulatedWavelength::skSolarSpectrum_TabulatedWavelength()
{
}


/*-----------------------------------------------------------------------------
 *					skSolarSpectrum_SAO2010::~skSolarSpectrum_SAO2010		 2015- 1- 6*/
/** **/
/*---------------------------------------------------------------------------*/

skSolarSpectrum_TabulatedWavelength::~skSolarSpectrum_TabulatedWavelength()
{
}


/*-----------------------------------------------------------------------------
 *					skSolarSpectrum_TabulatedWavelength::AttachToTable		 2015- 1- 6*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSolarSpectrum_TabulatedWavelength::AttachToTable	( const double* usertable, size_t numelements )
{
	size_t				numwavelen;
	bool				ok;
	nx2dArray< double>	irradiancetable;
	nx1dArray< double>	w;
	nx1dArray< double>	f;
	double*				table = (double *)usertable;	


	numwavelen  = numelements/2;
	
	ok = irradiancetable.Attach(2,numwavelen, table );
	irradiancetable.YSlice( 0, &w);
	irradiancetable.YSlice( 1, &f);
	m_wavelen_nm_invacuum.DeepCopy(w);
	m_irradiance.DeepCopy(f);
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skSolarSpectrum_TabulatedWavelength::AttachToTable		 2015- 1- 6*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSolarSpectrum_TabulatedWavelength::AttachToTable	( const double* userflux, const double* userwavelength, size_t numelements )
{
	bool	ok;
	nx1dArray< double>	w;
	nx1dArray< double>	f;
	double* flux        = (double *)userflux;
	double* wavelength  = (double *)userwavelength;

	ok =       w.Attach( numelements, wavelength);
	ok = ok && f.Attach(numelements, flux  );
	ok = ok && m_wavelen_nm_invacuum.DeepCopy( w );
	ok = ok && m_irradiance.DeepCopy( f );
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skSolarSpectrum_SAO2010::MinValidWavelength		 2015- 1- 6*/
/** **/
/*---------------------------------------------------------------------------*/

double skSolarSpectrum_TabulatedWavelength::MinValidWavelength() const
{
	return m_wavelen_nm_invacuum.front();
}


/*-----------------------------------------------------------------------------
 *					skSolarSpectrum_SAO2010::MaxValidWavelength		 2015- 1- 6*/
/** **/
/*---------------------------------------------------------------------------*/

double skSolarSpectrum_TabulatedWavelength::MaxValidWavelength() const
{
	return m_wavelen_nm_invacuum.back();
}
		


/*-----------------------------------------------------------------------------
 *					skSolarSpectrum_SAO2010::IrradianceAt1AU		 2015- 1- 6*/
/** **/
/*---------------------------------------------------------------------------*/

double skSolarSpectrum_TabulatedWavelength::IrradianceAt1AU( double wavelen_nm_vacuum )
{
	bool				ok;
	nxArrayIter<double>	iter;
	size_t				lowercell;
	size_t				uppercell;
	double				x0;
	double				x1;
	double				y[2];
	double				irradiance = std::numeric_limits<double>::quiet_NaN();

	ok =    (wavelen_nm_vacuum >= MinValidWavelength())
		 && (wavelen_nm_vacuum <= MaxValidWavelength());
	if (ok)
	{
		ok = nxLinearInterpolate::FindBoundingIndicesAscending<double, nxArrayIter<double> > ( m_wavelen_nm_invacuum.begin(),
																							   m_wavelen_nm_invacuum.end(),
																							   wavelen_nm_vacuum,
																							   &lowercell,
																							   &uppercell,
																							   &x0,
																							   &x1);
		if (ok)
		{
			y[0] = m_irradiance.At(lowercell);
			y[1] = m_irradiance.At(uppercell);
			irradiance = nxLinearInterpolate::FromTwoPoints( wavelen_nm_vacuum, x0, x1, y);
		}
	}
	return irradiance;
}

