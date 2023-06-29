//#include "sktran_hr_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_SZA		2014-10-30*/
/** A modification of the HR diffuse table to only allow placement of profiles
 *  at discrete solar zenith angles.  Interpolation between profiles is done
 *  linearly in cosine of solar zenith angle.
 *
 *  This diffuse table is used when the HR engine is operating in one dimensional
 *  mode, the only change is how interpolation is done between the profiles.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_HR_Diffuse_Table_SZA : public SKTRAN_HR_Diffuse_Table_CPU
{
	private:
		std::vector<double>							m_cosszas;

	private:
		bool										SZAWeights					( double sza, SKTRAN_HR_WEIGHT_TYPE* szaweights, size_t* szaindex, size_t& numindex ) const ;

	public:
													SKTRAN_HR_Diffuse_Table_SZA (){}
		virtual									   ~SKTRAN_HR_Diffuse_Table_SZA (){}

		virtual bool								ChooseDiffusePoints			( const HELIODETIC_POINT&				pt, 
																				  size_t*								diffuseindex, 
																				  SKTRAN_HR_WEIGHT_TYPE*				diffuseweights, 
																				  size_t*								numpoints ) const override;

		virtual bool								ChooseGroundPoints			( const HELIODETIC_POINT&				pt,
																				  size_t*								diffuseindex,
																				  SKTRAN_HR_WEIGHT_TYPE*				diffuseweights,
																				  size_t&								numpoints ) const override;

		virtual void								SetNumProfileInterp			( size_t numinterp ) override;

};


