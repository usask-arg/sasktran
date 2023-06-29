#include <skopticalproperties21.h>
#include "sources/sasktranif_opticalimpl/skbrdf_stubs.h"



/*-----------------------------------------------------------------------------
 *					ISKBrdf_Stub_Base		 2016- 12- 12*/
/** **/
/*---------------------------------------------------------------------------*/

ISKBrdf_Stub_LambertianAlbedo::ISKBrdf_Stub_LambertianAlbedo( SKTRAN_BRDF_Lambertian* lambertian)
	                          : ISKBrdf_Stub_Base( lambertian)
{
	m_lambertian = lambertian;
}


/*-----------------------------------------------------------------------------
 *					ISKBrdf_Stub_LambertianAlbedo::~ISKBrdf_Stub_LambertianAlbedo		 2016- 12- 12*/
/** **/
/*---------------------------------------------------------------------------*/

ISKBrdf_Stub_LambertianAlbedo::~ISKBrdf_Stub_LambertianAlbedo()
{}


/*-----------------------------------------------------------------------------
 *					ISKBrdf_Stub_LambertianAlbedo::SetPropertyScalar		 2016- 12- 12*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKBrdf_Stub_LambertianAlbedo::SetPropertyScalar( const char* propertyname, double value)
{
	nxString					name( propertyname);
	bool						ok = m_lambertian != nullptr;

	if (name == "Albedo")
	{
		ok = ok && m_lambertian->SetAlbedo( value );
		if (!ok)
		{
			nxLog::Record( NXLOG_WARNING,"ISKBrdf_Stub_LambertianAlbedo::SetProperty(Albedo), failed to set albedo to value %e", (value));
		}
	}
	else
	{
		ok = ISKBrdf_Stub_Base::SetPropertyScalar( propertyname,value);
	}
	return ok;

}



/*-----------------------------------------------------------------------------
 *					SKTRAN_BRDF_Lambertian::SKTRAN_BRDF_Lambertian		 2016- 12- 9*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_BRDF_Lambertian::SKTRAN_BRDF_Lambertian()
{
	m_albedo = std::numeric_limits<double>::quiet_NaN();
}

SKTRAN_BRDF_Lambertian::SKTRAN_BRDF_Lambertian(double albedo)
{
	m_albedo = albedo;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_BRDF_Lambertian::SetAlbedo		 2016- 12- 9*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_BRDF_Lambertian::SetAlbedo( double albedo )
{
	m_albedo = albedo;
	return NXFINITE(m_albedo);
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_BRDF_Lambertian::BRDF		 2016- 12- 9*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_BRDF_Lambertian::BRDF( double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double DPHI, double* brdf) const
{
	*brdf = m_albedo/nxmath::Pi;
	return NXFINITE(*brdf);
}




#if 0
/*-----------------------------------------------------------------------------
 *					SKTRAN_BRDF::BRDF		 2016- 11- 24*/
/** **/
/*---------------------------------------------------------------------------*/

bool  SKTRAN_BRDF::BRDF( const SKTRAN_SourceTermQueryObject_Base& incoming, 
						 const SKTRAN_SourceTermQueryObject_Base& outgoing, SKTRAN_Stokes_NC* source   ) const
{
	double costhetain;
	double costhetaout;
	const HELIODETIC_POINT&      point  = incoming.GetPoint();
	const HELIODETIC_UNITVECTOR& upH    = point.LocalZenith();
	HELIODETIC_UNITVECTOR		 localdirs[3];								// Gets local direction unit vectors: north, west, up
	const HELIODETIC_UNITVECTOR& indirH  = incoming.GetLookAway();
	const HELIODETIC_UNITVECTOR& outdirH = outgoing.GetLookAway();

	nxVector					 up		(     upH.X(),     upH.Y(),     upH.Z() );
	nxVector					 indir	(  indirH.X(),  indirH.Y(),  indirH.Z() );
	nxVector					 outdir ( outdirH.X(), outdirH.Y(), outdirH.Z() );	
	nxVector					 north  ( localdirs[0].X(), localdirs[0].Y(), localdirs[0].Z() );	// This not true geographic north but is the north which is towards the sun direction 
	nxVector					 west   ( localdirs[1].X(), localdirs[1].Y(), localdirs[1].Z() );	// NOt true geographic west bu 90 degrees from the sun (ie either dawn or dusk not sure which one)
;

	costhetain =  (up & indir);
	costhetaout = (up & outdir);
	NXASSERT(( costhetain >= 0.0  ));
	NXASSERT(( costhetaout >= 0.0 ));

	point.LocalUnitVectors( localdirs, N_ELEMENTS(localdirs) );
	double inx    = indir & north;
	double iny    = indir & west;
	double inazi  = nxmath::atan2d( iny,inx);
	double outx   = outdir & north;
	double outy   = outdir & west;
	double outazi = nxmath::atan2d( outy, outx);

	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_BRDF_Lambertian::SKTRAN_BRDF_Lambertian		 2016- 11- 24*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_BRDF_Lambertian::SKTRAN_BRDF_Lambertian()
{
	m_albedo = 0.3;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_BRDF_Lambertian::BRDF		 2016- 11- 24*/
/** **/
/*---------------------------------------------------------------------------*/

bool  SKTRAN_BRDF_Lambertian::BRDF( const SKTRAN_SourceTermQueryObject_Base& incoming, 
									const SKTRAN_SourceTermQueryObject_Base& outgoing, 
									SKTRAN_Stokes_NC*						 source) const 
{
	const HELIODETIC_POINT&      point  = incoming.GetPoint();
	const HELIODETIC_UNITVECTOR& up     = point.LocalZenith();
	const HELIODETIC_UNITVECTOR& indir  = incoming.GetLookAway();

	double costhetain =  (up & indir);
	source ->SetTo ( m_albedo*costhetain/nxmath::Pi );
	return true;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_BRDF_ThetaInOutAndDeltaAzi::SKTRAN_BRDF_ThetaInOutAndDeltaAzi		 2016- 11- 24*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_BRDF_ThetaInOutAndDeltaAzi::SKTRAN_BRDF_ThetaInOutAndDeltaAzi()
{
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_BRDF_ThetaInOutAndDeltaAzi::~SKTRAN_BRDF_ThetaInOutAndDeltaAzi		 2016- 11- 24*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_BRDF_ThetaInOutAndDeltaAzi::~SKTRAN_BRDF_ThetaInOutAndDeltaAzi()
{
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_BRDF_ThetaInOutAndDeltaAzi::BRDF		 2016- 11- 24*/
/** **/
/*---------------------------------------------------------------------------*/

bool  SKTRAN_BRDF_ThetaInOutAndDeltaAzi::BRDF( const SKTRAN_SourceTermQueryObject_Base& incoming, 
											   const SKTRAN_SourceTermQueryObject_Base& outgoing, 
											   SKTRAN_Stokes_NC* source   ) const
{
	bool   ok;
	double costhetain;
	double costhetaout;
	const HELIODETIC_POINT&      point  = incoming.GetPoint();
	const HELIODETIC_UNITVECTOR& upH    = point.LocalZenith();
	HELIODETIC_UNITVECTOR		 localdirs[3];								// Gets local direction unit vectors: north, west, up
	const HELIODETIC_UNITVECTOR& indirH  = incoming.GetLookAway();
	const HELIODETIC_UNITVECTOR& outdirH = outgoing.GetLookAway();

	nxVector					 up		(     upH.X(),     upH.Y(),     upH.Z() );
	nxVector					 indir	(  indirH.X(),  indirH.Y(),  indirH.Z() );
	nxVector					 outdir ( outdirH.X(), outdirH.Y(), outdirH.Z() );	
	nxVector					 north  ( localdirs[0].X(), localdirs[0].Y(), localdirs[0].Z() );	// This not true geographic north but is the north which is towards the sun direction 
	nxVector					 west   ( localdirs[1].X(), localdirs[1].Y(), localdirs[1].Z() );	// NOt true geographic west bu 90 degrees from the sun (ie either dawn or dusk not sure which one)


	costhetain =  (up & indir);
	costhetaout = (up & outdir);
	NXASSERT(( costhetain >= 0.0  ));
	NXASSERT(( costhetaout >= 0.0 ));

	point.LocalUnitVectors( localdirs, N_ELEMENTS(localdirs) );
	double inx    = indir & north;
	double iny    = indir & west;
	double inazi  = nxmath::atan2d( iny,inx);
	double outx   = outdir & north;
	double outy   = outdir & west;
	double outazi = nxmath::atan2d( outy, outx);
	double	deltaazi = nxmath::inrange( (outazi-inazi), 360.0 );

	ok = Brdf_impl( costhetain, costhetaout, deltaazi, source);
	return ok;
}

#endif

