//#include "sktran_hr_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_WF_Ray		2014-10-30*/
/**  Storage for the weighting function calculation for a single ray.
 *   
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_HR_WF_Ray_Base {

};

class SKTRAN_HR_WF_Ray : public SKTRAN_HR_WF_Ray_Base
{
	private:
		std::vector<double> m_weights;			//!< the weighting function values for a list of perturbations
		std::vector<double>	m_addedopticaldepth;
		double m_brdfwf;
	public:
		SKTRAN_HR_WF_Ray() { m_brdfwf = 0.0; };
		~SKTRAN_HR_WF_Ray() { };
		
		bool Allocate			( size_t numperturb );
		size_t  NumPerturb      () const { return m_weights.size(); }
		double& WeightAt		( size_t idx ) { return m_weights[idx]; }
		const double& WeightAt  ( size_t idx ) const { return m_weights[idx]; }
		double& AddedOpticalDepth( size_t idx ) { return m_addedopticaldepth[idx]; }
		double& brdfwf() { return m_brdfwf; }
};

class SKTRAN_HR_WF_Ray_Polarized : public SKTRAN_HR_WF_Ray_Base
{
private:
    std::vector<SKTRAN_Stokes_NC> m_weights;			//!< the weighting function values for a list of perturbations
    std::vector<double>	m_addedopticaldepth;
    double m_brdfwf;
public:
    SKTRAN_HR_WF_Ray_Polarized() { m_brdfwf = 0.0; };
    ~SKTRAN_HR_WF_Ray_Polarized() { };

    bool Allocate			( size_t numperturb );
    size_t  NumPerturb      () const { return m_weights.size(); }
    SKTRAN_Stokes_NC& WeightAt		( size_t idx ) { return m_weights[idx]; }
    const SKTRAN_Stokes_NC& WeightAt  ( size_t idx ) const { return m_weights[idx]; }
    double& AddedOpticalDepth( size_t idx ) { return m_addedopticaldepth[idx]; }
    double& brdfwf() { return m_brdfwf; }
};