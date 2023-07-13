#pragma once



/*-----------------------------------------------------------------------------
 *					class skBRDF		 2016- 12- 9*/
/** Implements the Bidirectional Reflectance Distribution function interface
 *	for the Sasktran model.
 *  
 *	There are some different conventions/assumptions made regarding BRDFs
 *	across the literature - I will list all of the options that I know about
 *	and state which ones are used in Sasktran.
 *
 *	Relative azimuth angle convention:
 *		A relative azimuth angle (raa)(phi) of 0 can be defined as forward or
 *		backward scattering - Sasktran defines it as the former. When 
 *		a model defines the raa as the latter, phi must be replaced with
*       pi - phi (equivalently cos(phi) must be replaced with -cos(phi))
 *		to match conventions.
 *
 *	Reflectance quantities:
 *		Several different quantities are often used to represent reflectance,
 *		and naming conventions are not always consistent. Hapke (2012, section
 *		10.2) provides a helpful naming convention for these quantities:
 *
 *		 - bidirectional reflectance: r = I / J
 *			 - (scattered radiance in a given direction) / (power per unit
 *			   area perpendicular to direction of incidence)
 *		 - bidirectional reflectance distribution function: BRDF = I / mu J
 *			 - (scattered radiance in a given direction) / (power per unit
 *			   area on surface)
 *		 - reflectance factor: REFF = pi I / mu J
 *			 - (reflectance of a surface) / (reflectance of a perfectly
 *			   diffuse function under same conditions)
 *		 - radiance factor: RADF = pi I / J
 *			 - (reflectance of a surface) / (reflectance of a perfectly
 *			   diffuse surface under normal incidence)
 *			    
 *		By these definitions, Sasktran uses BRDFs but most of the literature
 *		on remote sensing BRDFs uses REFFs, so their formulas must be scaled
 *		down by pi to match Sasktran's convention.
 *
 *	References:
 *	Hapke B. Theory of reflectance and emittance spectroscopy. Cambridge:
 *	Cambridge University Press, 1993.
 **/
/*---------------------------------------------------------------------------*/

class skBRDF : public nxUnknown
{
	public:

		/*-----------------------------------------------------------------------------
		 *					BRDF										2016- 12- 9*/
		/** Fetch the scalar BRDF for the given wavelength at the given location for 
		 *	the given ray geometry:
		 *
		 *	\param wavelennm
		 *	The wavelength in nanometers at which the BRDF is required. Many BRDF
		 *	implementations may choose to ignore this parameter if they implement a
		 *	wavelength independent BRDF.
		 *
		 *	\param pt
		 *	The location at which the BRDF is required. It is implicitly
		 *	assumed that this is the ground surface at the given latitude, longitude
		 *	and time, i.e. ignore the height field. Many BRDF implementations may choose
		 *	to ignore this parameter if they implement a location/time independent BRDF's
		 *	
		 *	\param MU_in
		 *	Cosine of the zenith angle of the incoming ray. Note that this is the
		 *	cosine of the angle between the local up direction and the incoming ray 
		 *	vector in the direction AWAY from the surface.
		 *	
		 *	\param MU_in
		 *	Cosine of the zenith angle of the out-going ray. Note that this is the
		 *	cosine of the angle between the local up direction and the out-going ray 
		 *	vector in the direction AWAY from the surface.
		 *
		 *	\param COSDPHI
		 *	Cosine of the azimuth angle between the incoming and outgoing ray. We use
		 *	the cosine as we need a value between 0 and \f$\pi\f$ and the cosine does this
		 *	while avoiding the need to specify degrees or radians. Note that the azimuth 
		 *	should be measured betwee the incoming ray vector TOWARD the surface and the 
		 *	outgoing ray vector AWAY from the surface.
		 *
		 *	\param brdf
		 *	Returns the calculated scalar brdf. Most implementations should return either the
		 *	correct value or NaN if there is a problem.
		 *
		 *	\return
		 *	True if success else false.
		 **/
		/*---------------------------------------------------------------------------*/
		virtual bool	BRDF( double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* brdf)  const= 0; /* MUST BE Thread safe for Sasktran version 2.1 */

		virtual bool	IsLambertian() const {return false;}

	protected:
		bool			CheckCosines(double &MU_in, double &MU_out, double &COSDPHI, nxString functionname) const; // forces MUs to correspond to 0 to 89 deg, and forces COSDPHI to be between -1 and 1 inclusive

};


/*-----------------------------------------------------------------------------
 *					SKTRAN_BRDF_Lambertian		 2016- 11- 24*/
/** **/
/*---------------------------------------------------------------------------*/

class SKTRAN_BRDF_Lambertian : public skBRDF
{
	private:
		double			m_albedo;

	public:
						SKTRAN_BRDF_Lambertian();
						SKTRAN_BRDF_Lambertian(double albedo);
		virtual		   ~SKTRAN_BRDF_Lambertian(){};

	public:
		bool			SetAlbedo			( double albedo);
		virtual bool	BRDF				( double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double DPHI, double* brdf) const override;

		virtual bool	IsLambertian() const override {
			return true;
		}
};



/*-----------------------------------------------------------------------------
 *					class SKTRAN_BRDF_Roujean_entry		 2016- 12- 8*/
/** **/
/*---------------------------------------------------------------------------*/

class SKTRAN_BRDF_Roujean_entry 
{
	public:
		double	k0;
		double	k1;
		double	k2;
		double	lambda_nm;

	public:
				SKTRAN_BRDF_Roujean_entry( double ak0, double ak1, double ak2, double alambda_nm )	{ k0 = ak0; k1 = ak1; k2 = ak2; lambda_nm = alambda_nm;}
				SKTRAN_BRDF_Roujean_entry()															{ k0 = 0;   k1 = 0;   k2 = 0;   lambda_nm = 0.0;}
};

/*-----------------------------------------------------------------------------
 *					class SKTRAN_BRDF_Roujean					2016- 12- 8*/
/** Implements the BRDF function described by Roujean et al. The model is
 *	designed to model vegetation surfaces. It contains three adjustable
 *	parameters and the model considers the surface reflectance to consist
 *	of two main processes  (i) a diffuse reflection component and (ii) a 
 *	volume  scattering by a collection of dispersed facets.
 *	
 *  Equation 10 is implemented with values of k0, k1, and k2 from Table 1, 
 *  but with the following changes:
 *	 - cos(phi) is replaced with -cos(phi)
 *	 - the reflectance is divided by pi
     - the reflectance is also divided by 100 (k0, k1, k2 are percentages)
 *
 *	Reference:
 *	Roujean, J.-L., M. Leroy, and P.-Y. Deschamps (1992), A bidirectional
 *	reflectance model of the Earth's surface for the correction of remote
 *	sensing data, J. Geophys. Res., 97(D18), 20455–20468, doi:10.1029/92JD01411.
 **/
/*---------------------------------------------------------------------------*/


class SKTRAN_BRDF_Roujean : public skBRDF
{
	private:
		double		m_k0;
		double		m_k1;
		double		m_k2;
		static		std::map< std::string, SKTRAN_BRDF_Roujean_entry>				m_predefinedsurfaces;
		typedef		std::map< std::string, SKTRAN_BRDF_Roujean_entry>::value_type	value_type;

	private:
		bool			CheckPredefinedSurfaces		();

	public:
						SKTRAN_BRDF_Roujean			();
		virtual		   ~SKTRAN_BRDF_Roujean			() {}
		bool			LoadPredefinedParameters	( const char* paramname );
		bool			SetBRDFParameters			( double k0, double k1, double k2);
		virtual bool	BRDF						( double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double DPHI, double* brdf) const override;
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_BRDF_Snow_Kokhanovsky2012		 2016- 12- 12*/
/** This is an analytical BRDF model for snow publisshed in 2012 by Kokhanovsky and
 *	Breon. The model has 2 free parameters (L and M) which the user should provide. The
 *	authors suggest deriving the free parameters from fitting data at several wavelengths.
 *	This implementation provides default values for the two parameters based upon snow measured
 *	Greenland in May 2006.
 *	
 *	Equation 1 is implemented with the following changes:
 *	 - the reflectance is divided by pi
 *	
 *	References:
 *	A. A. Kokhanovsky and F. M. Breon, "Validation of an Analytical Snow BRDF Model
 *	Using PARASOL Multi-Angular and Multispectral Observations," 
 *	IEEE Geoscience and Remote Sensing Letters, vol. 9, no. 5, pp. 928-932, Sept. 2012.
 *  doi: 10.1109/LGRS.2012.2185775
 *  URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6166850&isnumber=6205680

 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_BRDF_Snow_Kokhanovsky2012 : public skBRDF
{
	private:												//   The Kokhanovsky and Breon paper of 2012 identifies 2 free parameters for the snow BRDF which should be fitted from observations
		double						m_L;					//!< Parameter (1): L is approximately 13d where d is the average optical diameter of snow grains.
		double						m_M;					//!< Parameter (2): M is directly proportional to mass concentration of pollutants in the snow.
		skRTRefractiveIndex_ICE		m_iceri;				//!< Refractive index of ice.

	private:
		double			K0										( double mu) const;
		double			p										( double thetadegrees) const;
		double			R0										( double mus, double muv, double theta) const;
		double			Chi										( double wavelen_nm ) const;

	public:
						SKTRAN_BRDF_Snow_Kokhanovsky2012		();
		virtual		   ~SKTRAN_BRDF_Snow_Kokhanovsky2012		() {}
		bool			SetBRDFParameters						( double L, double M );
		virtual bool	BRDF									( double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* brdf) const override;
};



/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Roujean_Kernel		 2017-08-02                          */
/** This kernel is derived considering a random arrangement of rectangular
*	blocks on a flat surface. It is one of the kernels used by
*	SKTRAN_BRDF_Roujean.
*
*	This is one of eight BRDF kernels listed by Spurr in
*   his 2002 paper. There is a discrepancy between Spurr Equation A.1 and 
*	Roujean Equation 2: the first term in Spurr A.1 should be divided by 2pi.
*	
*	Equation 2 from Roujean was implemented with the following changes:
*	 - cos(phi) is replaced with -cos(phi)
*	 - the reflectance is divided by pi
*
*	References:
*	Robert J. D. Spurr, "A new approach to the retrieval of surface properties
*	from earthshine measurements," Journal of Quantitative Spectroscopy
*	& Radiative Transfer, vol. 83, pp. 15-46, Oct. 9, 2002.
*	doi:10.1016/S0022-4073(02)00283-2
*	url: http://www.sciencedirect.com/science/article/pii/S0022407302002832
*
*	Roujean, J.-L., M. Leroy, and P.-Y. Deschamps (1992), A bidirectional
*	reflectance model of the Earth's surface for the correction of remote
*	sensing data, J. Geophys. Res., 97(D18), 20455–20468, doi:10.1029/92JD01411.
**/
/*---------------------------------------------------------------------------*/
class SKTRAN_BRDF_Roujean_Kernel : public skBRDF
{
public:
					SKTRAN_BRDF_Roujean_Kernel();
	virtual			~SKTRAN_BRDF_Roujean_Kernel() {};

public:
	virtual bool	BRDF(double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* brdf) const override;

};


/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Li_Kernel		 2017-03-06                          */
/** The Li-sparse and Li-Dense kernels were published by Wanner and Strahler in 1995.
*	They are BRDF kernels that model the reflectance of forest canopies. Both
*	kernels share two non-linear parameters: the crown shape ratio b/r and 
*	the relative height ratio h/b (trees are modelled as oblate spheroids of 
*	height b and width r with their centres at height h off the ground). 
*
*	There is a discrepancy between Spurr's Equation A.2/A.3 and Wanner Equation
*	35/50: the final term in d(x) in Spurr A.2/A.3 should be added, not
*	subtracted.
*
*	A variation of the Li-sparse kernel called the Li-sparse-reciprocal
*	kernel is in common use - a slight change to the formula makes
*	the kernel reciprocal (the incomcing and outgoing angles can be
*	swapped without changing the value).
*
*	These are two of eight BRDF kernels listed by Spurr in
*   his 2002 paper.
*	
*	Equation 32 (sparse) and 47 (dense) from Wanner and Equation 39 from 
*	Strahler (sparse reciprocal) are implemented with the following changes:
*	 - cos(phi) is replaced with -cos(phi)
*	 - the reflectance is divided by pi
*
*	References:
*	W. Wanner and A. H. Strahler, "On the derivation of kernels for kernel-driven 
*	modles of bidirectional reflectance," Journal of Geophysical Research, vol. 
*	100, no. d10, pp. 21077-21089, Oct. 20, 1995.
*	doi: 10.1029/95JD02371
*	URL: http://onlinelibrary.wiley.com/doi/10.1029/95JD02371/abstract
*
*	Robert J. D. Spurr, "A new approach to the retrieval of surface properties
*	from earthshine measurements," Journal of Quantitative Spectroscopy 
*	& Radiative Transfer, vol. 83, pp. 15-46, Oct. 9, 2002.
*	doi:10.1016/S0022-4073(02)00283-2
*	url: http://www.sciencedirect.com/science/article/pii/S0022407302002832
*
*	A. H. Strahler, J. P. Muller, "MODIS BRDF/Albedo Product: Algorithm
*	Theoretical Basis Document Version 5.0," April 1999.
*	URL: https://lpdaac.usgs.gov/sites/default/files/public/product_documentation/atbd_mod09_v5.pdf
**/
/*---------------------------------------------------------------------------*/
class SKTRAN_BRDF_Li_Kernel : public skBRDF
{
	protected:			// kernel describe by two non-linear parameters:
		double	m_br;	// crown shape ratio b/r
		double	m_hb;	// relative height ratio h/b

	protected:
		double	overlap(double mus, double muv, double mup) const;					// eqn 33 and 48
		double	cos_primed_scattering_angle(double s, double v, double cp) const;	// eqn 36 and 51
		double	primed_angle(double theta) const;

	public:
				SKTRAN_BRDF_Li_Kernel();
				virtual	~SKTRAN_BRDF_Li_Kernel() {};
		bool	SetBRDFParameters(double br, double hb);							// eqn 37 and 52
};

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_LiSparse_Kernel		 2017-03-06*/
/** **/
/*---------------------------------------------------------------------------*/
class SKTRAN_BRDF_LiSparse_Kernel : public SKTRAN_BRDF_Li_Kernel
{
	public:
		virtual bool BRDF(double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* brdf) const override;
};

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_LiSparseReciprocal_Kernel		 2017-07-31*/
/** The original Li Sparse kernel is not reciprocal (swapping the incoming and
*	outgoing angles changes the value). The Li Sparse Reciprocal kernel
*	fixes this by changing one term. See the MODIS ATBD referenced above.
/*---------------------------------------------------------------------------*/
class SKTRAN_BRDF_LiSparseReciprocal_Kernel : public SKTRAN_BRDF_Li_Kernel
{
public:
	virtual bool BRDF(double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* brdf) const override;
};


/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_LiDense_Kernel			 2017-03-06*/
/** **/
/*---------------------------------------------------------------------------*/
class SKTRAN_BRDF_LiDense_Kernel : public SKTRAN_BRDF_Li_Kernel
{
	public:
		virtual bool BRDF(double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* brdf) const override;
};



/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Ross_Kernel		                           */
/** 
*	The Ross-thick and Ross-thin kernels are volume-scattering kernels that
*	model the reflection of surface facets (leaves). It is purely geometric
*	and depends on no additional parameters.
*
*	These are two of eight BRDF kernels listed by Spurr in
*   his 2002 paper.
*	
*	Equation 7 (thick) and 13 (thin) from Wanner are implemented with the
*	following changes:
*	 - cos(phi) is replaced with -cos(phi)
*	 - the reflectance is divided by pi
*
*	References:
*	W. Wanner and A. H. Strahler, "On the derivation of kernels for kernel-driven
*	modles of bidirectional reflectance," Journal of Geophysical Research, vol.
*	100, no. d10, pp. 21077-21089, Oct. 20, 1995.
*	doi: 10.1029/95JD02371
*	URL: http://onlinelibrary.wiley.com/doi/10.1029/95JD02371/abstract
*
*	Robert J. D. Spurr, "A new approach to the retrieval of surface properties
*	from earthshine measurements," Journal of Quantitative Spectroscopy
*	& Radiative Transfer, vol. 83, pp. 15-46, Oct. 9, 2002.
*	doi:10.1016/S0022-4073(02)00283-2
*	url: http://www.sciencedirect.com/science/article/pii/S0022407302002832
**/
/*---------------------------------------------------------------------------*/
class SKTRAN_BRDF_Ross_Kernel : public skBRDF
{
	public:
					SKTRAN_BRDF_Ross_Kernel();
		virtual		~SKTRAN_BRDF_Ross_Kernel() {}

	protected:
		double	cos_primed_scattering_angle(double s, double v, double cp) const;
};

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_RossThin_Kernel			 2017-03-06*/
/** **/
/*---------------------------------------------------------------------------*/
class SKTRAN_BRDF_RossThin_Kernel : public SKTRAN_BRDF_Ross_Kernel
{
	public: 
		virtual bool BRDF(double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* brdf) const override;
};

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_RossThick_Kernel			 2017-03-06*/
/** **/
/*---------------------------------------------------------------------------*/
class SKTRAN_BRDF_RossThick_Kernel : public SKTRAN_BRDF_Ross_Kernel
{
public:
	virtual bool BRDF(double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* brdf) const override;
};

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Hapke		          2017-07-28             */
/**
*	The Hapke BRDF depends on three parameters and is usually used
*	independently, not in a linear combination with other kernels. It is
*	similar to the Rahman kernel but the hotspot, a peak in
*	reflectivity back in the direction of the incident light due to the
*	absense of shadowing effects, is treated differently.
*
*	This is one of eight BRDF kernels listed by Spurr in
*   his 2002 paper.
*
*	Equation A.10 of Spurr is implemented with the following changes:
*	 - cos(phi) is replaced with -cos(phi)
*	 - the reflectance is divided by pi
*
*	References:
*	Hapke B. Theory of reflectance and emittance spectroscopy. Cambridge: 
*	Cambridge University Press, 1993.
*
*	Robert J. D. Spurr, "A new approach to the retrieval of surface properties
*	from earthshine measurements," Journal of Quantitative Spectroscopy
*	& Radiative Transfer, vol. 83, pp. 15-46, Oct. 9, 2002.
*	doi:10.1016/S0022-4073(02)00283-2
*	url: http://www.sciencedirect.com/science/article/pii/S0022407302002832

**/
/*---------------------------------------------------------------------------*/

class SKTRAN_BRDF_Hapke : public skBRDF
{
private:
	double			m_omega;
	double			m_delta;
	double			m_b0;

public:
					SKTRAN_BRDF_Hapke();
	virtual		   ~SKTRAN_BRDF_Hapke() {};

public:
	bool			SetBRDFParameters(double omega, double delta, double b0);
	virtual bool	BRDF(double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double DPHI, double* brdf) const override;
};

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_CoxMunk		       2017-07-27            */
/**
*
*	This BRDF is based on the work of Cox and Munk who
*	analyzed the surface of water as a function of wind using aerial
*	photographs in the 1950s. Their work includes dependence on wind
*	direction and a Gram-Charlier correction to the Gaussian distribution
*	of water facet orientations, but this BRDF neglects these corrections
*	and assumes a symmetric distribution. This BRDF depends on two
*	parameters: wind speed and the index of refraction of water.
*
*	This is one of eight BRDF kernels listed by Spurr in
*   his 2002 paper. 
*	
*	Equation A.18 from Spurr is implemented with the following changes:
*	 - cos(phi) is replaced with -cos(phi)
*	 - the reflectance is divided by pi
*
*	References:
*	Cox C, Munk W. The measurement of the roughness of the sea surface 
*	from photographs of the Sun glitter. J Opt Soc Ann 1954; 44:838–50.
*
*	Robert J. D. Spurr, "A new approach to the retrieval of surface properties
*	from earthshine measurements," Journal of Quantitative Spectroscopy
*	& Radiative Transfer, vol. 83, pp. 15-46, Oct. 9, 2002.
*	doi:10.1016/S0022-4073(02)00283-2
*	url: http://www.sciencedirect.com/science/article/pii/S0022407302002832

**/
/*---------------------------------------------------------------------------*/

class SKTRAN_BRDF_CoxMunk : public skBRDF
{
private:
	double			m_wind_speed;
	double			m_index_of_refraction;

public:
					SKTRAN_BRDF_CoxMunk();
	virtual		   ~SKTRAN_BRDF_CoxMunk() {};

public:
	bool			SetBRDFParameters(double wind_speed, double index_of_refraction);
	virtual bool	BRDF(double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double DPHI, double* brdf) const override;
};

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_Rahman		           2017-07-28          */
/**
*
*	This is a semiempirical model with 3 parameters, designed
*	to model the reflectance of arbitrary natural surfaces in the visible
*	and near-infrared bands.	
*
*	This is one of eight BRDF kernels listed by Spurr in
*   his 2002 paper. 
*	
*	Equation 2 from Rahman is implemented with the following changes:
*	 - cos(phi) is replaced with -cos(phi)
*	 - the reflectance is divided by pi
*
*	References:
*	Rahman H, Pinty B, Verstraete M. Coupled surface atmosphere reflectance 
*	(CSAR) model: 2. Semi-empirical surface model usable with NOAA AVHRR
*	data. J Geophys Res 1993; 98:20,791–801.
*
*	Robert J. D. Spurr, "A new approach to the retrieval of surface properties
*	from earthshine measurements," Journal of Quantitative Spectroscopy
*	& Radiative Transfer, vol. 83, pp. 15-46, Oct. 9, 2002.
*	doi:10.1016/S0022-4073(02)00283-2
*	url: http://www.sciencedirect.com/science/article/pii/S0022407302002832

**/
/*---------------------------------------------------------------------------*/

class SKTRAN_BRDF_Rahman : public skBRDF
{
private:
	double			m_rho0;
	double			m_k;
	double			m_theta;

public:
					SKTRAN_BRDF_Rahman();
	virtual		   ~SKTRAN_BRDF_Rahman() {};

public:
	bool			SetBRDFParameters(double rho0, double theta, double k);
	virtual bool	BRDF(double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double DPHI, double* brdf) const override;
};

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_LinearCombinationBase	   2017-08-02          */
/**
*	Derived classes allows an arbitrary number of BRDF kernels to be used in
*	a linear combination.
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_BRDF_LinearCombinationBase : public skBRDF
{
protected:
	std::vector<double>					m_f;
	std::vector<skBRDF*>				m_kernels;

public:
	virtual bool	BRDF(double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double DPHI, double* brdf) const override;
};

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_LinearCombination		           2017-07-28          */
/**
*	This class allows an arbitrary number of BRDF kernels to be used in
*	a linear combination.
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_BRDF_LinearCombination : public SKTRAN_BRDF_LinearCombinationBase
{
public:
					SKTRAN_BRDF_LinearCombination();
	virtual		   ~SKTRAN_BRDF_LinearCombination() {};

public:
	int				NumKernels();
	bool			AddKernel(skBRDF* kernel);
	bool			RemoveKernel(int index);
	bool			SetKernelWeights(const double* f, int numpoints);
};

/*-----------------------------------------------------------------------------
*		SKTRAN_BRDF_MODIS_RossThickLiSparseReciprocal        2017-07-31      */
/**
*	This is the BRDF model used by MODIS for their BRDF/albedo product. It
*	consists of a linear combination of an isotropic (lambertian) term, a
*	volume scattering term (Ross-thick), and a geometric scattering term
*	(Li-Sparse-Reciprocal).
*	
*	References:
*	Robert J. D. Spurr, "A new approach to the retrieval of surface properties
*	from earthshine measurements," Journal of Quantitative Spectroscopy
*	& Radiative Transfer, vol. 83, pp. 15-46, Oct. 9, 2002.
*	doi:10.1016/S0022-4073(02)00283-2
*	url: http://www.sciencedirect.com/science/article/pii/S0022407302002832

**/
/*---------------------------------------------------------------------------*/

class SKTRAN_BRDF_MODIS_RossThickLiSparseReciprocal : public SKTRAN_BRDF_LinearCombinationBase
{
public:
					SKTRAN_BRDF_MODIS_RossThickLiSparseReciprocal();
	virtual		   ~SKTRAN_BRDF_MODIS_RossThickLiSparseReciprocal() {};

public:
	bool			SetBRDFParameters(double f_iso, double f_vol, double f_geo);
};