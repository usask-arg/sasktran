// msis90.cpp : Implementation of nxmsis90
#include <skclimatology21.h>
#include <nxbase_threads.h>


static std::mutex	g_mutex;
//static nxMutex			g_mutex;


#if defined(NX_WINDOWS)
	#define FORTRAN_GTD6	GTD6
#else
	#define FORTRAN_GTD6	gtd6_
#endif

extern "C" void  FORTRAN_GTD6(  nxLONG*		IYD,		// 4 byte integer
							double*		SEC,
							double*		ALT,
							double*		GLAT,
							double*		GLONG,
							double*		STL,
							double*		F107A,
							double*		F107,
							double*		AP,
							nxLONG*		MASS,			// 4 byte integer
							double*		D,
							double*		T);

/*---------------------------------------------------------------------------
 *'					nxmsis90::nxmsis90                                2002-7-17
 *-------------------------------------------------------------------------*/

nxmsis90::nxmsis90()
{
	m_mjd.MJD(51000.0);
	m_altitude  = 80.0;			// altitude in kms
	m_latitude  = 52.0;			// Geodetic latitude
	m_longitude = -106.0;		// Geodetic longitude
	m_f107a     = 150.0;		// F107, F107A, and AP effects are not large below 80 km
	m_f107      = 150.0;		// and these can be set to 150., 150., and 4. respectively.
	for (int i = 0; i < 7; i++) m_ap[i] = 4.0;

	m_HE               = 0;
	m_O                = 0;
	m_N2               = 0;
	m_O2               = 0;
	m_AR               = 0;
	m_TOTALMASS        = 0;
	m_H                = 0;
	m_N                = 0;
	m_exosphere_temp   = 0;
	m_altitude_temp    = 0;
	m_isdirty          = nxTRUE;
}


/*-----------------------------------------------------------------------------
 *					nxmsis90::DeepCopy		2008-5-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxmsis90::DeepCopy( const nxmsis90& other )
{
	m_mjd				= other.m_mjd;
	m_altitude			= other.m_altitude;			// altitude in kms
	m_latitude			= other.m_latitude;			// Geodetic latitude
	m_longitude			= other.m_longitude;		// Geodetic longitude
	m_f107a				= other.m_f107a;
	m_f107				= other.m_f107;
	m_ap[0]				= other.m_ap[0];
	m_ap[1]				= other.m_ap[1];
	m_ap[2]				= other.m_ap[2];
	m_ap[3]				= other.m_ap[3];
	m_ap[4]				= other.m_ap[4];
	m_ap[5]				= other.m_ap[5];
	m_ap[6]				= other.m_ap[6];
	m_isdirty			= other.m_isdirty;
	m_HE				= other.m_HE;
	m_O					= other.m_O;
	m_N2				= other.m_N2;
	m_O2				= other.m_O2;
	m_AR				= other.m_AR;
	m_TOTALMASS			= other.m_TOTALMASS;
	m_H					= other.m_H;
	m_N					= other.m_N;
	m_exosphere_temp	= other.m_exosphere_temp;
	m_altitude_temp		= other.m_altitude_temp;
	return true;
}


/*
C        Neutral Atmosphere Empirical Model from the surface to lower
C          exosphere  MSISE90 (JGR, 96, 1159-1172, 1991)
C         A.E.Hedin 4/24/90;6/3/91(add SAVE)
C         2/11/93 correct switch initialization and mks calculation
C       2/11/97 [AEH] CORRECT ERROR IN GHP6 WHEN USING METER6(.TRUE.)
C           See subroutine GHP6 to specify a pressure rather than
C           altitude.
C     INPUT:
C        IYD - YEAR AND DAY AS YYDDD or DDD (day of year from 1 to 365)
C        SEC - UT(SEC)
C        ALT - ALTITUDE(KM)
C        GLAT - GEODETIC LATITUDE(DEG)
C        GLONG - GEODETIC LONGITUDE(DEG)
C        STL - LOCAL APPARENT SOLAR TIME(HRS)
C        F107A - 3 MONTH AVERAGE OF F10.7 FLUX
C        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
C        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
C           - ARRAY CONTAINING:
C             (1) DAILY AP
C             (2) 3 HR AP INDEX FOR CURRENT TIME
C             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
C             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
C             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
C             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
C                    TO CURRENT TIME
C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 59 HRS PRIOR
C                    TO CURRENT TIME
C        MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
C                 CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
C     Note:  Ut, Local Time, and Longitude are used independently in the
C            model and are not of equal importance for every situation.
C            For the most physically realistic calculation these three
C            variables should be consistent (STL=SEC/3600+GLONG/15).
C            F107, F107A, and AP effects are not large below 80 km
C            and these can be set to 150., 150., and 4. respectively.
C     OUTPUT:
C        D(1) - HE NUMBER DENSITY(CM-3)
C        D(2) - O NUMBER DENSITY(CM-3)
C        D(3) - N2 NUMBER DENSITY(CM-3)
C        D(4) - O2 NUMBER DENSITY(CM-3)
C        D(5) - AR NUMBER DENSITY(CM-3)
C        D(6) - TOTAL MASS DENSITY(GM/CM3)
C        D(7) - H NUMBER DENSITY(CM-3)
C        D(8) - N NUMBER DENSITY(CM-3)
C        T(1) - EXOSPHERIC TEMPERATURE
C        T(2) - TEMPERATURE AT ALT
C
C      TO GET OUTPUT IN M-3 and KG/M3:   CALL METER6(.TRUE.)
C
C      O, H, and N set to zero below 72.5 km
C      Exospheric temperature set to average for altitudes below 120 km.
C
C           The following is for test and special purposes:
C            TO TURN ON AND OFF PARTICULAR VARIATIONS CALL TSELEC(SW)
C               WHERE SW IS A 25 ELEMENT ARRAY CONTAINING 0. FOR OFF, 1.
C               FOR ON, OR 2. FOR MAIN EFFECTS OFF BUT CROSS TERMS ON
C               FOR THE FOLLOWING VARIATIONS
C               1 - F10.7 EFFECT ON MEAN  2 - TIME INDEPENDENT
C               3 - SYMMETRICAL ANNUAL    4 - SYMMETRICAL SEMIANNUAL
C               5 - ASYMMETRICAL ANNUAL   6 - ASYMMETRICAL SEMIANNUAL
C               7 - DIURNAL               8 - SEMIDIURNAL
C               9 - DAILY AP             10 - ALL UT/LONG EFFECTS
C              11 - LONGITUDINAL         12 - UT AND MIXED UT/LONG
C              13 - MIXED AP/UT/LONG     14 - TERDIURNAL
C              15 - DEPARTURES FROM DIFFUSIVE EQUILIBRIUM
C              16 - ALL TINF VAR         17 - ALL TLB VAR
C              18 - ALL TN1 VAR           19 - ALL S VAR
C              20 - ALL TN2 VAR           21 - ALL NLB VAR
C              22 - ALL TN3 VAR           23 - TURBO SCALE HEIGHT VAR
C
C              To get current values of SW: CALL TRETRV(SW)
*/
/////////////////////////////////////////////////////////////////////////////
// nxmsis90

nxBOOL nxmsis90::InvokeMsis90()
{
	nxBOOL	ok =nxTRUE;
	nxLONG	iyd;
	double	sec;
	double	alt;
	double	glat;
	double	glon;
	double	stl;
	double	D[8];
	double	T[2];
	double	utsecs;
	int		day,month,year,hour,mins,secs;
	double	ticks;
	nxLONG   mass = 48;

	if (m_isdirty)
	{
		m_mjd.GetUTC( &day, &month, &year, &hour, &mins, &secs, &ticks);
		iyd    = (nxLONG)( year*1000 + m_mjd.DayOfYear());
		utsecs = (hour*3600.0 + mins*60.0 + secs + ticks);
		sec    = (float)utsecs;
		stl    = (float)( utsecs/3600 + m_longitude/15);
		alt    = (float)m_altitude;
		glat   = (float)m_latitude;
		glon   = (float)m_longitude;

		{
			std::unique_lock<std::mutex>	lock( g_mutex );			// lock the mutex, It will unlock when it goes out of scope
			FORTRAN_GTD6(  &iyd, &sec, &alt, &glat, &glon, &stl, &m_f107a, &m_f107, &m_ap[0], &mass, &D[0], &T[0]);
		}
		m_HE             = (double)D[0];
		m_O              = (double)D[1];
		m_N2             = (double)D[2];
		m_O2             = (double)D[3];
		m_AR             = (double)D[4];
		m_TOTALMASS      = (double)D[5];
		m_H              = (double)D[6];
		m_N              = (double)D[7];
		m_exosphere_temp = (double)T[0];
		m_altitude_temp  = (double)T[1];
		m_isdirty        = nxFALSE;
	}
	return ok;
}

void nxmsis90::UpdateCachedProfiles( double mjd, double latitude, double longitude )
{
	SetMjd( mjd );
	SetLatitude(latitude);
	SetLongitude(longitude);
}

void nxmsis90::SetMjd(double newVal)
{
	m_mjd.MJD(newVal);
	m_isdirty = nxTRUE;
}


void nxmsis90::SetHeight(double newVal)
{
	m_altitude = newVal;
	m_isdirty = nxTRUE;
}


void nxmsis90::SetLatitude(double newVal)
{
	m_latitude = newVal;
	m_isdirty = nxTRUE;
}


void nxmsis90::SetLongitude(double newVal)
{
	m_longitude = newVal;
	m_isdirty = nxTRUE;
}


void nxmsis90::SetF107Avg(double newVal)
{
	m_f107a = (float)newVal;
	m_isdirty = nxTRUE;
}


/*---------------------------------------------------------------------------
 *                        nxmsis90::SetAP                         2019-11-22 */
/** **/
/*---------------------------------------------------------------------------*/

bool nxmsis90::SetAP( const double* ap, int n)
{
	bool ok;

	ok = (n ==7);
	if (ok)
	{
		for (int i=0; i < 7; i++)
		{
			m_ap[i] = ap[i];
		}
	}
	return ok;
}

void nxmsis90::SetF107(double newVal)
{
	m_f107 = (float)newVal;
	m_isdirty = nxTRUE;
}

double nxmsis90::HE()
{
	InvokeMsis90();
	return m_HE;
}

double nxmsis90::H()
{
	InvokeMsis90();
	return m_H;
}

double nxmsis90::N()
{
	InvokeMsis90();
	return m_N;
}

double nxmsis90::AR()
{
	InvokeMsis90();
	return m_AR;
}

double nxmsis90::N2()
{
	InvokeMsis90();
	return m_N2;
}

double nxmsis90::O()
{
	InvokeMsis90();
	return m_O;
}

double nxmsis90::O2()
{
	InvokeMsis90();
	return m_O2;
}

double nxmsis90::O2_O2()
{
	InvokeMsis90();
	return m_O2*m_O2;
}

double nxmsis90::TotalMass()
{
	InvokeMsis90();
	return m_TOTALMASS;
}

double nxmsis90::T_Exosphere()
{
	InvokeMsis90();
	return m_exosphere_temp;
}

double nxmsis90::T()
{
	InvokeMsis90();
	return m_altitude_temp;
}

/*<--------------------------------------------------------------------------
 *'					nxmsis90::MeanMolecularWeight                  2002-10-17
 *	Returns the mean molecular weight at the current location.  Below 80.0
 *	km it simply returns a constant average molecular weight. Above 80 km it
 *	works it out from the abundance of each species.
 *>------------------------------------------------------------------------*/

double nxmsis90::MeanMolecularWeight()
{
	const double M0 = 28.9644;					// Mean molecular weight of air, US. Standard atmosphere 1976
	if (m_altitude  <= 80.0) return M0;			// Mean molecular weight below 80 is fixed.

	InvokeMsis90();
	double 	msum = 0.0;
	double	nsum = 0.0;

/*
	msum += 4.002602*m_HE;
	msum += 15.9994*m_O;
	msum += 28.01348*m_N2;
	msum += 31.9988*m_O2;
	msum += 39.948*m_AR;
	msum += 1.00794*m_H;
	msum += 14.00674*m_N;
*/
	msum =    4.0*m_HE			// Use the same Molecular weights as MSIS Fortran code
		   + 16.0*m_O			// Less accurate but is consistent with the TotalMass calculation
		   + 28.0*m_N2			// and gives a better estimate of pressure
		   + 32.0*m_O2
		   + 40.0*m_AR
		   + 14.0*m_N;

	nsum = m_HE + m_O + m_N2 + m_O2 + m_AR + m_H + m_N;
	return msum/nsum;
}

/*<--------------------------------------------------------------------------
 *'					nxmsis90::P                                    2002-10-17
 *	Calculates the pressure at the current location.  Uses the TotalMass,
 *	Mean Molecular Weight and Ideal gas law to get the pressure.
 *
 *	Returns the pressure in dynes.
 *>------------------------------------------------------------------------*/

double nxmsis90::P()
{
	InvokeMsis90();
    return (TotalMass()/MeanMolecularWeight())*(nxcgs::KBOLTZMAN/nxcgs::AMU)*nxmsis90::T();
}

/*<--------------------------------------------------------------------------
 *'					nxmsis90::MeanNumberDensity                    2002-10-21
 *>------------------------------------------------------------------------*/

double nxmsis90::MeanNumberDensity()
{
	InvokeMsis90();
	return TotalMass()/(MeanMolecularWeight()*nxcgs::AMU);
}



/*---------------------------------------------------------------------------
 *'					nxmsis90cached::nxmsis90cached                 2002-10-23
 *-------------------------------------------------------------------------*/

nxmsis90cached::nxmsis90cached()
{
	m_minH		= 0.0;								// Minimum altitude in KM
	m_maxH      = 120.0;							// maximum altitude in KM
	m_deltah	= 1.0;								// altitude resolution in KMS
	m_cacheisdirty = true;

	AddSpecies(SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3);
	AddSpecies(SKCLIMATOLOGY_TEMPERATURE_K);
	AddSpecies(SKCLIMATOLOGY_PRESSURE_PA);
	AddSpecies(SKCLIMATOLOGY_O2_CM3);
	AddSpecies(SKCLIMATOLOGY_O2_O2_CM6);
}

/*---------------------------------------------------------------------------
 *'					nxmsis90cached::SetCacheAltitudeRange          2002-10-23
 *-------------------------------------------------------------------------*/

void nxmsis90cached::SetCacheAltitudeRange	( double minh, double maxh, double deltah )
{
	m_minH		= minh;							// Minimum altitude in KM
	m_maxH      = maxh;							// maximum altitude in KM
	m_deltah	= deltah;						// altitude resolution in KMS
	SetDirtyCache();
}


/*---------------------------------------------------------------------------
 *                nxmsis90cached::SetMaxHeightKMS                 2019-11-22 */
/** **/
/*---------------------------------------------------------------------------*/

bool nxmsis90cached::SetMaxHeightKMS( double maxh )
{
	bool ok;
	ok = (m_maxH >= 1.0) && (m_maxH < 1001.0);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"nxmsis90cached::SetMaxHeightKMS, The maximum height for MSIS must lie between 1 and 1000 kms. Received value [%e]. Nothing has been changed", (double)maxh);
	}
	else
	{
		m_maxH = maxh;							// maximum altitude in KM
	}
	SetDirtyCache();
	return ok;
}

/*---------------------------------------------------------------------------
 *              nxmsis90cached::SetHeightSpacingKMS               2019-11-22 */
/** **/
/*---------------------------------------------------------------------------*/

bool nxmsis90cached::SetHeightSpacingKMS( double deltah)
{
	bool ok;
	ok = (deltah >= 0.0009) && ( deltah <= 20.0);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"nxmsis90cached::SetHeightSpacingKMS, The height spacing for MSIS must lie between 0.001 and 20.0 kms. Received value [%e]. Nothing has been changed", (double)deltah);
	}
	else
	{
		m_deltah = deltah;					// altitude resolution in KMS
	}
	SetDirtyCache();
	return ok;
}

/*---------------------------------------------------------------------------
 *               nxmsis90cached::AddSpeciesProfile                2019-11-22 */
/** **/
/*---------------------------------------------------------------------------*/

bool nxmsis90cached::AddSpecies( const CLIMATOLOGY_HANDLE& key)
{
	bool	ok;

	
	ok = !(m_profiles.find(key) == m_profiles.end());					// See if we already fetch this species
	if (!ok)															// If we dont the lets add it
	{																	// but it must be one that we support
		ok =     (key == SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3)
			  || (key == SKCLIMATOLOGY_TEMPERATURE_K)
			  || (key == SKCLIMATOLOGY_PRESSURE_PA)
			  || (key == SKCLIMATOLOGY_O2_O2_CM6)
			  || (key == SKCLIMATOLOGY_He_CM3)
			  || (key == SKCLIMATOLOGY_O_CM3)
			  || (key == SKCLIMATOLOGY_N2_CM3)
			  || (key == SKCLIMATOLOGY_O2_CM3)
			  || (key == SKCLIMATOLOGY_Ar_CM3)
			  || (key == SKCLIMATOLOGY_H_CM3)
			  || (key == SKCLIMATOLOGY_N_CM3);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"nxmsis90cached::AddSpecies, The requested species is unsupported by MSIS. It has not been added");
		}
		else
		{
			auto status = m_profiles.insert( value_type( key, nxSpline2()) );
			ok = status.second;
		}
		SetDirtyCache();
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *                    nxmsis90cached::SetF10p7                    2019-11-22 */
/** **/
/*---------------------------------------------------------------------------*/
bool nxmsis90cached::SetF10p7( double f107)
{
	SetF107( f107);
	SetDirtyCache();
	return true;
}

/*---------------------------------------------------------------------------
 *                   nxmsis90cached::SetF10p7a                    2019-11-22 */
/** **/
/*---------------------------------------------------------------------------*/
bool nxmsis90cached::SetF10p7avg( double f107a)
{
	SetF107Avg(f107a);					//!< Set the value of F10.7
	SetDirtyCache();
	return true;
}


/*---------------------------------------------------------------------------
 *                     nxmsis90cached::SetAp                      2019-11-22 */
/** **/
/*---------------------------------------------------------------------------*/

bool nxmsis90cached::SetAp( const double* ap, int numpts)
{
	bool ok;

	ok = SetAP( ap, numpts);
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING,"nxmsis90cached::SetAp, error setting Ap with [%d} points. Make sure you have supplied an array of 7 sensible values", (int)numpts);
	}
	SetDirtyCache();
	return ok;
}

/*---------------------------------------------------------------------------
 *                nxmsis90cached::ConfigureSpline                 2019-11-22 */
/** **/
/*---------------------------------------------------------------------------*/

void nxmsis90cached::ConfigureSpline( const CLIMATOLOGY_HANDLE& key, const nx1dArray<double>&	ht, nx1dArray<double>&	profile)
{
	iterator	iter = m_profiles.find(key);
	if (iter != m_profiles.end())
	{
		nxSpline2*	spline = &(iter->second);
		spline->Configure( ht, profile );
	}
}


/*---------------------------------------------------------------------------
 *               nxmsis90cached::IsSupportedSpecies               2019-11-22 */
/** **/
/*---------------------------------------------------------------------------*/

bool nxmsis90cached::IsSupportedSpecies ( const CLIMATOLOGY_HANDLE& species )
{
	return ( m_profiles.find(species) != m_profiles.end() );
}

/*---------------------------------------------------------------------------
 *              nxmsis90cached::InterpolateToHeight               2019-11-22 */
/** **/
/*---------------------------------------------------------------------------*/

bool nxmsis90cached::InterpolateToHeight( const CLIMATOLOGY_HANDLE& key, double h, double* value )
{
	bool	ok;

	iterator iter = m_profiles.find(key);
	ok = (iter != m_profiles.end());
	if(ok)
	{
		nxSpline2*	spline = &(iter->second);
		*value = spline->Interpolate( h );
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *'					nxmsis90cached::UpdateCachedProfiles           2002-10-23
 *-------------------------------------------------------------------------*/

void nxmsis90cached::UpdateCachedProfiles( double mjd, double latitude, double longitude )
{
	int	npts;
	double				Ht;
	nx1dArray<double>	ht;

	nx1dArray<double>	He;
	nx1dArray<double>	O;
	nx1dArray<double>	N2;
	nx1dArray<double>	O2;
	nx1dArray<double>	O4;
	nx1dArray<double>	Ar;
	nx1dArray<double>	H;
	nx1dArray<double>	N;

	nx1dArray<double>	P;
	nx1dArray<double>	T;
	nx1dArray<double>	n;
	bool				isuptodate;

	isuptodate =    (mjd == Mjd() ) 
		         && (latitude == Latitude() )
				 && (longitude == Longitude() )
				 && (!CacheIsDirty());

	if (!isuptodate)
	{
		SetMjd      ( mjd );
		SetLatitude ( latitude );
		SetLongitude( longitude );

		npts = (int)( (m_maxH- m_minH)/m_deltah ) + 1;

		ht.SetSize (npts);

		He.SetSize(npts);
		O.SetSize(npts);
		N2.SetSize(npts);
		O2.SetSize(npts);
		O4.SetSize(npts);
		Ar.SetSize(npts);
		H.SetSize(npts);
		N.SetSize(npts);

		P.SetSize (npts);
		T.SetSize (npts);
		n.SetSize (npts);

		for (int i=0; i < npts; i++)
		{
			Ht = m_minH + i*m_deltah;
			SetHeight(Ht);
			InvokeMsis90();
			ht[i]  = Ht;

			He[i] = nxmsis90::HE();
			O[i]  = nxmsis90::O();
			N2[i] = nxmsis90::N2();
			O2[i] = nxmsis90::O2();
			O4[i] = nxmsis90::O2_O2();
			Ar[i] = nxmsis90::AR();
			H[i]  = nxmsis90::H();
			N[i]  = nxmsis90::N();

			P[i]  = nxmsis90::P();
			T[i]  = nxmsis90::T();
			n[i]  = nxmsis90::MeanNumberDensity();
		}
		ConfigureSpline(SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3, ht, n  );
		ConfigureSpline(SKCLIMATOLOGY_TEMPERATURE_K,        ht, T );
		ConfigureSpline(SKCLIMATOLOGY_PRESSURE_PA,          ht, P  );
		ConfigureSpline(SKCLIMATOLOGY_O2_O2_CM6,            ht, O4 );
		ConfigureSpline(SKCLIMATOLOGY_He_CM3,				ht, He );
		ConfigureSpline(SKCLIMATOLOGY_O_CM3,				ht, O  );
		ConfigureSpline(SKCLIMATOLOGY_N2_CM3,				ht, N2 );
		ConfigureSpline(SKCLIMATOLOGY_O2_CM3,				ht, O2 );
		ConfigureSpline(SKCLIMATOLOGY_Ar_CM3,				ht, Ar );
		ConfigureSpline(SKCLIMATOLOGY_H_CM3,				ht, H  );
		ConfigureSpline(SKCLIMATOLOGY_N_CM3,				ht, N  );
		m_cacheisdirty = false;
	}
}

/*-----------------------------------------------------------------------------
 *					nxmsis90cached::DeepCopy		2008-5-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxmsis90cached::DeepCopy( const nxmsis90cached& other )
{
	bool ok;

	m_profiles  = other.m_profiles;
	ok          = nxmsis90::DeepCopy( other );
	m_minH	    = other.m_minH;											// Minimum altitude in KM
	m_maxH	    = other.m_maxH;											// maximum altitude in KM
	m_deltah	= other.m_deltah;										// altitude resolution in KMS
	m_cacheisdirty = true;
	return ok;
}




