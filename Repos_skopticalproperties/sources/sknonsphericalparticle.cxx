#include <skopticalproperties21.h>

using namespace std;


/*-----------------------------------------------------------------------------
 *					sk_NonsphericalParticle::sk_NonsphericalParticle		2009-10-15*/
/** **/
/*---------------------------------------------------------------------------*/

sk_NonsphericalParticle::sk_NonsphericalParticle()
{
	m_delta			= 1e-3;
	m_qext			= 0;
	m_qsca			= 0;
	m_qabs			= 0;
	m_Cext			= 0;
	m_Csca			= 0;
	m_Cabs			= 0;
	m_isdirty		= true;
	m_anyang		= false;
	m_aspectRatio   = 2;
	m_radius        = 0.0;
	m_lambda        = 0.0;
	m_refr          = 0.0;
	m_k             = 0.0;
	m_xx            = 0.0;
	m_anyang        = false;		//!< If TRUE then any angle can be placed in m_xmu. If FALSE then angles are mirror symmetric about 90
}


/*-----------------------------------------------------------------------------
 *					sk_NonsphericalParticle::~sk_NonsphericalParticle		2009-10-15*/
/** **/
/*---------------------------------------------------------------------------*/

sk_NonsphericalParticle::~sk_NonsphericalParticle()
{
	m_xmu.erase();
}

/*-----------------------------------------------------------------------------
 *					sk_NonsphericalParticle::operator=			2009-6-2	*/
/*---------------------------------------------------------------------------*/

bool sk_NonsphericalParticle::DeepCopy( const sk_NonsphericalParticle& other )
{
	bool	ok;

	m_isdirty		= true;
	m_radius		= other.m_radius;
	m_lambda		= other.m_lambda;
	m_refr			= other.m_refr;
	m_k				= other.m_k;
	m_xx			= other.m_xx;
	m_delta			= other.m_delta;
	m_qext			= other.m_qext;
	m_qsca			= other.m_qsca;
	m_qabs			= other.m_qabs;
	m_Cext			= other.m_Cext;
	m_Csca			= other.m_Csca;
	m_Cabs			= other.m_Cabs;
	m_anyang		= other.m_anyang;
	m_aspectRatio   = other.m_aspectRatio;
	ok				= m_xmu.DeepCopy( other.m_xmu );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"sk_NonsphericalParticle::DeepCopy, Error copying contents, this object is in an ill-defined state and may produce strange results");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					sk_NonsphericalParticle::SetCrossSections		2009-10-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool sk_NonsphericalParticle::SetCrossSectionsFromFortran( double cext_microns, double csca_microns )
{
	double		area;

	area   = nxmath::Pi*m_radius*m_radius;
	m_qext = cext_microns/area;					// extinction efficiency
	m_qsca = csca_microns/area;					// scattering efficiency
	m_qabs = m_qext - m_qsca;					// absorption efficiency

	m_Cext = cext_microns*1.0e-8;				// extinction cross-section (in units of cm**2)
	m_Csca = csca_microns*1.0e-8;				// scattering cross-section (in units of cm**2)
	m_Cabs = m_Cext - m_Csca;					// absorption cross-section (in units of cm**2)
	return true;
}
		

/*-----------------------------------------------------------------------------
 *					sk_NonsphericalParticle::Set_ComputationAccuracy		2009-10-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool sk_NonsphericalParticle::Set_ComputationAccuracy( double delta )
{
	SetDirty( m_delta != delta );
	m_delta			= delta;
	return true;
}

/*-----------------------------------------------------------------------------
 *					sk_NonsphericalParticle::Set_Wavelength		2009-10-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool sk_NonsphericalParticle::Set_Wavelength( double lamda )
{
	SetDirty( m_lambda != lamda );
	m_lambda		= lamda;
	UpdateXX();
	return true;
}

/*-----------------------------------------------------------------------------
 *					sk_NonsphericalParticle::Set_Radius		2009-10-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool sk_NonsphericalParticle::Set_Radius( double radius )
{
	SetDirty( m_radius != radius );
	m_radius		= radius;
	UpdateXX();
	return true;
}

/*-----------------------------------------------------------------------------
 *					sk_NonsphericalParticle::Set_RefractiveIndex		2009-10-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool sk_NonsphericalParticle::Set_RefractiveIndex( double ri_real, double ri_imag )
{
	complex<double>	val( ri_real, ri_imag);

	SetDirty( m_refr != val );
	m_refr			= val;
	return true;
}



/*-----------------------------------------------------------------------------
 *					sk_NonsphericalParticle::AllocateArrays		2009-10-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool sk_NonsphericalParticle::AllocateArrays( size_t numangles )
{
	bool	ok;

	ok =    (numangles < 901)
		  && m_xmu.SetSize( numangles )
	      && AllocateScatteringAngleMatrixElements ( numangles );

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "sk_NonsphericalParticle::AllocateArrays, Error allocating arrays for %d angles (maximum number of angles is 901)", (int)numangles);
		m_xmu.erase();
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					sk_NonsphericalParticle::Set_SymmetricScatterAngles		2009-10-15*/
/** Define the scattering angles which must be calculated. The user can provide
 *	a list of angular sections with different resoutions in each section.
 *	This allows the forward scatter peak near 0 degrees to be accurately
 *	modeled at high resolution while the other angles can use a much coarser
 *	resolution.
 *
 *	\param startangle
 *		dblarr[numsteps]. 
 *		Specifies the start angle of eachregion in degrees. Must be between
 *		0 and less than 90.  I implicitly assume that the startangle of one
 *		section is greater than the endangle of the last section
 **/
/*---------------------------------------------------------------------------*/

bool sk_NonsphericalParticle::Set_SymmetricScatterAngles( double *startangle,
														 double *endangle,
													     double *resolution,
													     int	 numrange )
{
	double		lastend = 0.0;
	int			numangles;
	int			nsteps;
	int			totalangles;
	bool		ok = true;
	int			i;
	int			j;
	double		theta;
	int			fidx;
	int			bidx;

	m_anyang = false;
	numangles = 0;
	for (i=0; i < numrange && ok; i++)
	{
		ok = (startangle[i] >= lastend) && (endangle[i] <= 90.0);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "sk_NonsphericalParticle::Set_SymmetricScatterAngles, The startangle (%g) must be less than lastangle (%g) and endangle (%g) must be <= 90 degrees", (double)startangle[i], (double)lastend, (double)endangle[i]);
		}
		else
		{
			nsteps     = (int)((endangle[i] - startangle[i])/resolution[i] + 0.5);
			numangles += nsteps;
			lastend    = startangle[i] + (nsteps-1)*resolution[i];
		}
	}
	if (ok)
	{
		totalangles = numangles *2;								// Account for 90 to 180 degree sector
		totalangles++;											// last angle is not exactly 90 (this is always true)
		ok = AllocateArrays(totalangles);						// Allocate the space for the cosine of scattering angles
		numangles = 0;
		fidx      = 0;
		bidx      = totalangles-1;
		for (i=0; i < numrange && ok; i++)
		{
			nsteps   = (int)((endangle[i] - startangle[i])/resolution[i] + 0.5);
			for (j = 0; j < nsteps; j++)
			{
				theta = nxmath::cosd(startangle[i] + j*resolution[i]);
				m_xmu.At(fidx++) =  theta;
				m_xmu.At(bidx--) = -theta;
			}
		}
		m_xmu.At(fidx) = 0;
	}
	if (!ok) m_xmu.erase();
	SetDirty(true);
	return ok;
}

/*-----------------------------------------------------------------------------
 *					sk_NonsphericalParticle::Set_AnyScatterAngles		2009-10-15*/
/**	Define the scattering angles which must be calculated. The user just specifies
 *	the scattering angles.  They are not necessarily symmetric
 *
 *	\param startangle
 *		dblarr[numsteps]
 *		Specifies the start angle of eachregion in degrees. Must be between
 *		0 and less than 90.  I implicitly assume that the startangle of one
 *		section is greater than the endangle of the last section
**/
/*---------------------------------------------------------------------------*/

bool sk_NonsphericalParticle::Set_AnyScatterAngles( nx1dArray<double>& degrees)
{
	size_t			numangles;
	bool		ok;

	m_anyang  = true;
	numangles = degrees.size();
	ok        = AllocateArrays(numangles);						// Allocate the space for the cosine of scattering angles
	if (ok)
	{
		for (size_t i=0; i < numangles; i++)
		{
			m_xmu.At(i) = nxmath::cosd( degrees[i] );
		}
	}
	else m_xmu.erase();
	SetDirty(true);
	return ok;
}
