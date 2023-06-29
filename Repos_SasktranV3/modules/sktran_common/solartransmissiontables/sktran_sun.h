//#pragma once
//#include <modules/sktran_common/include/sktran_common_internals.h>

/*-----------------------------------------------------------------------------
 *					class SKTRAN_Sun_Base					2013-09-27*/
/**	@ingroup sun
 *	Base class for classes representing the sun. This will allow us to 
 *  treat the sun as a point, circle, etc., and optionally include effects
  * such as limb darkening. These probably aren't significant for the
  * OSIRIS mission, but are important for exoplanet RT simulations. 
 **/
/*---------------------------------------------------------------------------*/
class SKTRAN_Sun_Base 
{
	public:
                                             SKTRAN_Sun_Base   ( );
		virtual                             ~SKTRAN_Sun_Base   ( );
		//virtual bool                         MakeThreadsafeFor ( size_t numthreads=0 )                                     = 0;
		virtual void                         UpdateSun         ( )                                const = 0;
		virtual const HELIODETIC_UNITVECTOR& GetSunUnit        ( )                                const = 0;
		virtual void                         SunUnitVector     (       HELIODETIC_UNITVECTOR* h ) const = 0;
		virtual double                       CosAngleToSun     ( const HELIODETIC_UNITVECTOR& h ) const = 0;
		virtual double                       ComponentToSun    ( const HELIODETIC_VECTOR&     h ) const = 0;
};

/*-----------------------------------------------------------------------------
 *					class SKTRAN_Sun_Point					2013-09-27*/
/**	@ingroup sun
 *	Treat the sun as a point in direction HELIODETIC_UNITVECTOR(0,0,1). 
 **/
/*---------------------------------------------------------------------------*/
class SKTRAN_Sun_Point : public SKTRAN_Sun_Base
{
	protected:
		HELIODETIC_UNITVECTOR                m_sunUnit;

	public:
                                             SKTRAN_Sun_Point  ( );
		virtual                             ~SKTRAN_Sun_Point  ( );
		//virtual bool                         MakeThreadsafeFor ( size_t numthreads=0 );
		virtual void                         UpdateSun         ( )                               const override;
		virtual const HELIODETIC_UNITVECTOR& GetSunUnit        ( )                               const override;
		virtual void                         SunUnitVector     (       HELIODETIC_UNITVECTOR* h ) const override;
		virtual double                       CosAngleToSun     ( const HELIODETIC_UNITVECTOR& h ) const override;
		virtual double                       ComponentToSun    ( const HELIODETIC_VECTOR&     h ) const override;
};
