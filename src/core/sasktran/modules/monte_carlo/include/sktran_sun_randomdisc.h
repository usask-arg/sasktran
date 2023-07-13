#pragma once

#include "sktran_montecarlo_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_Sun_RandomDisc_ThreadEntry		 2014- 11- 17*/
/** **/
/*---------------------------------------------------------------------------*/

class SKTRAN_Sun_RandomDisc_ThreadEntry
{
	public:
		HELIODETIC_UNITVECTOR     m_sun;
		SKTRAN_RNG                m_rng;
};

/*-----------------------------------------------------------------------------
 *					class SKTRAN_Sun_RandomDisc					2013-09-27*/
/**	Treat the sun as a disc with apex angle (i.e. angle between unit vector
 *  pointing to disc center and unit vector pointing to disc edge) m_apexAngle. 
 *  The sun is treated as a point sun initially located at 
 *  HELIODETIC_UNITVECTOR(0,0,1) and moving to a randomly chosen point on 
 *  the disc with each call of UpdateSun(). 
 *  The sun is assumed to be a Lambertian emitter with no limb darkening. 
 **/
/*---------------------------------------------------------------------------*/
class SKTRAN_Sun_RandomDisc : public SKTRAN_Sun_Base
{
	private:
		nxThreadStorageMap< SKTRAN_Sun_RandomDisc_ThreadEntry >		m_currentSun1;
		HELIODETIC_UNITVECTOR*               m_currentSun;
		double                               m_sineApexAngle;
		const std::vector<SKTRAN_RNG>*                    m_rngs;

	public:
                                             SKTRAN_Sun_RandomDisc	( );
		virtual                             ~SKTRAN_Sun_RandomDisc	( );
		bool								 InitializeThreadEntry	( SKTRAN_Sun_RandomDisc_ThreadEntry* threadentry) const;
		bool								 InitializeSunVector	( HELIODETIC_UNITVECTOR* sun ) const;
		bool                                 Initialize				( double sineApexAngle, const std::vector<SKTRAN_RNG>& randGens, const size_t numThreads );
		void                                 ReleaseResources		( );
		double                               GetSineApexAngle		( )                                const;

		//virtual bool                         MakeThreadsafeFor     ( size_t numthreads=0 );
		virtual void                         UpdateSun				( ) const override;
		virtual const HELIODETIC_UNITVECTOR& GetSunUnit				( )                                const override;
		virtual void                         SunUnitVector			(       HELIODETIC_UNITVECTOR* h ) const override;
//		virtual void                         SunUnitVector			(       HELIODETIC_VECTOR&     h ) const override;
		virtual double                       CosAngleToSun			( const HELIODETIC_UNITVECTOR& h ) const override;
//		virtual double                       CosAngleToSun			( const HELIODETIC_VECTOR&     h ) const override;
		virtual double                       ComponentToSun			( const HELIODETIC_VECTOR&     h ) const override;

};
