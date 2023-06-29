#include <complex>

/*-----------------------------------------------------------------------------
 *					class sk_NonsphericalParticle		2009-6-3*/
/**  \ingroup miescat
 *	Abstract base class that represents a scattering and absorbing nonspherical
 *	particle. Derived classes are DLL wrapper classes for the T-matrix
 *	(oriented and randomly-oriented) and discrete dipole approximation (DDA)
 *	FORTRAN codes.
 **/
/*---------------------------------------------------------------------------*/

class sk_NonsphericalParticle : public nxUnknown
{
	public:
		enum ENUM_SHAPE					{ SHAPE_HEX, SHAPE_CYLINDER, SHAPE_SPHEROID };
		enum ENUM_ALGORITHM				{ ALG_TMATRAND, ALG_TMATORIENT, ALG_DDA };

	private:
		bool							m_isdirty;
		double							m_radius;		//!< equivalent-sphere radius: units are microns
		double							m_aspectRatio;	//!< aspect ratio: across-face/length b/w parallel faces
		double							m_lambda;		//!< Wavelength: units are microns
		std::complex<double>			m_refr;			//!< The complex refractive index (real,imaginary)
		double							m_k;			//!< 2Pi/lambda
		double							m_xx;			//!< The size parameter (2*pi*reff/wavelength)
		double							m_delta;		//!< accuracy of computations
		nx1dArray<double>				m_xmu;			//!< The array of cos(scattering angles) (0-Pi radians, always, monotone odd and always includes pi/2)
		double							m_qext;			//!< extinction efficiency
		double							m_qsca;			//!< scattering efficiency
		double							m_qabs;			//!< absorption efficiency
		double							m_Cext;			//!< extinction cross section (in units of cm2)
		double							m_Csca;			//!< scattering cross section (in units of cm2)
		double							m_Cabs;			//!< absorption cross section (in units of cm2)
		bool							m_anyang;		//!< If TRUE then any angle can be placed in m_xmu. If FALSE then angles are mirror symmetric about 90

	private:
		void							UpdateXX						() { m_k = nxmath::TWOPI/m_lambda; m_xx = m_k*m_radius;}
		bool							AllocateArrays					( size_t numangles );

	protected:
		void							SetDirty						( bool val )				{ m_isdirty = m_isdirty || val;}
		void							ClearDirty						( )							{ m_isdirty = false;}
		bool							IsDirty							( )							{ return m_isdirty;}
		bool							SetCrossSectionsFromFortran		( double cext_microns, double csca_microns );

	// ---- Protected pure virtual members

	protected:
		virtual bool					AllocateScatteringAngleMatrixElements	( size_t numangles )			= 0;	//!< Allocates extra storage arrays when the scattering angles are defined

	// ---- Public Pure virtual members

	public:
		virtual ENUM_SHAPE				Shape							() const								 = 0;
		virtual ENUM_ALGORITHM			Algorithm						() const								 = 0;
		virtual nxString				Get_DescriptiveParticleString	()	const								 = 0;
		virtual bool					CreateClone						( sk_NonsphericalParticle** clone) const = 0;
		virtual bool					CalculateScattering						( )								 = 0;	//!< Calls the fortran code and sends the cross-sections back to this object via SetCrossSectionsFromFortran
		virtual bool					Get_ScatteringMatrixCoeffs		( std::vector<skRTPhaseMatrix>*	coeffmatrix )		 = 0;			//!< Get Scattering matrix (related to the phase matrix

	// ---- Public members
	public:
										sk_NonsphericalParticle			();
		virtual						   ~sk_NonsphericalParticle			();
		bool							DeepCopy						( const sk_NonsphericalParticle& other );
		bool							Set_AspectRatio					( double ratio )					{ SetDirty( m_aspectRatio != ratio ); m_aspectRatio = ratio; return true;}
		bool							Set_Wavelength					( double lamda );
		bool							Set_Radius						( double radius );
		bool							Set_Delta						( double delta )					{ SetDirty(delta != m_delta); m_delta = delta; return true;}
		bool							Set_RefractiveIndex				( double ri_real, double ri_imag );
		bool							Set_SymmetricScatterAngles		( double *startangle, double *endangle, double *resolution, int numsteps);
		bool							Set_AnyScatterAngles			( nx1dArray<double>& degrees);
		bool							Set_ComputationAccuracy			( double delta );
		const nx1dArray< double>&		Get_CosAngles					( )	const { return m_xmu;}
		size_t							Get_NumAngles					( ) const { return m_xmu.size();}
		double							Get_Wavelength					( ) const { return m_lambda;}
		double							Get_Radius						( ) const { return m_radius;}
		double							Get_AspectRatio					( ) const { return m_aspectRatio;}
		double							Get_SizeParameter				( )	const { return m_xx;}
		double							Get_Lambda						( ) const { return m_lambda;}
		double							Get_Delta						( ) const { return m_delta;}
		double							Get_k							( ) const { return m_k;}
		std::complex<double>			Get_RefractiveIndex				( ) const { return m_refr;}
		double							Qext							( )		  { CalculateScattering(); return m_qext;}					//!< Get extinction efficiency
		double							Qsca							( )	      { CalculateScattering(); return m_qsca;}					//!< Get scattering efficiency
		double							Qabs							( )	      { CalculateScattering(); return m_qabs;}					//!< Get absoptionefficiency
		double							Cext							( )	      { CalculateScattering(); return m_Cext;}					//!< Get extinction cross-section (area in units of cm2)
		double							Csca							( )	      { CalculateScattering(); return m_Csca;}					//!< Get scattering cross-section (area in units of cm2)
		double							Cabs							( )		  { CalculateScattering(); return m_Cabs;}
};


/*-----------------------------------------------------------------------------
 *					class sk_NonsphericalParticleDistributedSizes			2009-6-8*/
/**  \ingroup miescat
 *	Intermediate class to call the Mishchenko T-Matrix FORTRAN code for
 *	randomly-oriented nonspherical particles. The purpose of this class is to match
 *	the 'extinction per particle' scheme of skOpticalProperties to the inherent
 *	size-distribution calculation done within the program TMATRIXRANDOM.dll.
 *
 *	For more details on the variables, please see the source code for the
 *	randomly-oriented T-matrix code (lmq.lp.f and assoc. files)
**/
/*---------------------------------------------------------------------------*/

class sk_NonsphericalParticleDistributedSizes : public sk_NonsphericalParticle
{
	protected:
		double							m_r1;			//!< minimum: size dist. quadrature
		double							m_r2;			//!< maximum: size dist. quadrature
		double							m_b;			//!< The second size distribution parameter passed to mischenko (his parameter B)
		double							m_gam;			//!< third size distribution parameter
		int								m_npnax;		//!< number of size distributions to use (always use 1)
		int								m_ndistr;		//!< which size distribution to use
		int								m_nkmax;		//!< number of PSD Gaussian quadrature points
		double							m_ssalb;		//!< single-scatter albedo
		double							m_asym;			//!< Asymmetry factor

	protected:
		bool			Set_SizeDistIntegrationLimits	();

	public:
						sk_NonsphericalParticleDistributedSizes( );
		bool			DeepCopy						( const sk_NonsphericalParticleDistributedSizes& other);
		void			Set_SSALB						( double ssalb ) { m_ssalb = ssalb;}
		void			Set_ASYM						( double asym  ) { m_asym  = asym;}
		void			Set_NumPartSizePoints			( int num_pts );
		void			Set_SizeDistLogNormal			( double radius, double sigma );
		void			Set_SizeDistModGamma			( double rc, double alpha,double gamma );
		void			Set_SizeDistPowerLaw			( double radius,double veff );
		void			Set_SizeDistGamma				( double lambda,double mu );
		void			Set_SizeDistMonodisperse		( double radius );
//		bool			Set_UncommonParameters			( double RAT,double r1,double r2, int NDISTR,int NP,int NPNAX, int NKMAX, int NDGS )	= 0;//!< TMatrixRandom only: allow user to set these parameters directly for testing purposes

};

 
/*-----------------------------------------------------------------------------
 *					class sk_TMatrixRandomWrapper						2009-6-3 */
/**  \ingroup miescat
 *	This is a wrapper class to call the Mishchenko T-Matrix FORTRAN code for
 *	randomly-oriented nonspherical particles, which has been compiled into
 *	TMATRIXRANDOM.dll. I have compiled the FORTRAN code largely "as is", only
 *	commenting some lines of screen-output print statements.
 *
 *	The T-matrix Fortran code is compiled under the tmatrix directory
 *	of the skopticalproperties/sources/nonspherescatter/tmatrixrandom subproject
 *	directory.
 *
 *	Reference:
 *	M.I. Mishchenko and L.D. Travis, Capabilities and limitations of a current
 *	FORTRAN implementation of the T-matrix method for randomly oriented,
 *	rotationally symmetric scatterers. J. Quant. Spect. Rad. Transf. 60, 3,
 *	309-324, 1998.
 *
 *	Also see Michael Mishchenko's web site
 *  Dr. Michael I. Mishchenko (crmim@giss.nasa.gov)
 *		http://www.giss.nasa.gov/staff/mmishchenko/t_matrix.html
 *        NASA Goddard Institute for Space Studies
 *
 *	Note the following:
 *	Extinction cross-section, Cext = Qext . pi.a^2		a = radius of volume-equivalent sphere to nonspherical particle
 *  Scattering cross-section, Csca = Qsca x pi.a^2
 *	Absorbtion cross-section = (Cext - Csca)
**/
/*---------------------------------------------------------------------------*/

class sk_TMatrixRandomWrapper : public sk_NonsphericalParticleDistributedSizes
{
	private:
		int								m_np;			//!< particle shape: spheroids NP=-1, cylinders NP=-2
		int								m_ndgs;			//!< no. division points over particle surface
		double							m_rat;			//!< RAT=1: part size specified by equal-volume sphere radius r_v. RAT!=1: equal-surf. area sphere radius r_s
		nx1dArray<double>				m_f11;			//!< F11 phase matrix element array
		nx1dArray<double>				m_f22;			//!< F22 phase matrix element array
		nx1dArray<double>				m_f33;			//!< F33 phase matrix element array
		nx1dArray<double>				m_f44;			//!< F44 phase matrix element array
		nx1dArray<double>				m_f12;			//!< F12 phase matrix element array
		nx1dArray<double>				m_f34;			//!< F34 phase matrix element array

	private:
		sk_TMatrixRandomWrapper&		operator=						( const sk_TMatrixRandomWrapper& other );	// Dummy declaration not actually implemented, stops people doing assignments
		bool							AllocateArrays					( size_t numangles );
		bool							Mishchenko_TMatrix				();
		const char*						PartTypeString					()	const;
		bool							SetDefaultAngles				( );

	protected:
		virtual bool					AllocateScatteringAngleMatrixElements	( size_t numangles);

	public:
										sk_TMatrixRandomWrapper			();
		virtual						   ~sk_TMatrixRandomWrapper			();
		bool							DeepCopy						( const sk_TMatrixRandomWrapper& other );
		bool							Set_ShapeCylinder				();
		bool							Set_ShapeSpheroid				();
		bool							Set_SizeEqualVolumeSphere		()					{ m_rat =  1.0; return true;}
		bool							Set_SizeEqualSurfAreaSphere		()					{ m_rat = -2.0; return true;}
		void							Set_NumRadiusQuadraturePoints	( int numpoints )	{ m_nkmax = numpoints; SetDirty(true);}

		virtual ENUM_SHAPE				Shape							() const		{ return SHAPE_CYLINDER; }
		virtual ENUM_ALGORITHM			Algorithm						() const	{ return ALG_TMATRAND; }
		virtual nxString				Get_DescriptiveParticleString	() const;						//!<
		virtual bool					CreateClone						( sk_NonsphericalParticle** clone ) const;
		virtual bool					CalculateScattering				();
		virtual bool					Set_ParticleDistribution		( skRTParticleDist* newdistribution );
		virtual bool					Get_ScatteringMatrixCoeffs		( std::vector<skRTPhaseMatrix>*	coeffmatrix );
};


 /*---------------------------------------------------------------------------
 *					class sk_TMatrixOrientedWrapper						2009-6-3 */
/**	 \ingroup miescat
 *	This is a wrapper class to call the Mishchenko T-Matrix FORTRAN code,
 *	which has been compiled into TMATRIXORIENTED.dll and TMATRIXORIENTED.lib. I
 *	have compiled the FORTRAN code largely "as is", only commenting some lines
 *	of screen-output print statements.
 *
 *	The T-matrix Fortran code is compiled under the tmatrix directory
 *	of the skopticalproperties/sources/nonspherescatter/tmatrixrandom subproject
 *	directory.
 *
 *	Reference:
 *	M.I. Mishchenko, Calculation of the amplitude matrix for a nonspherical
 *	particle in a fixed orientation, Appl. Opt. 39, 6,1026-1031, 2000.
 *
 *	Also see Michael Mishchenko's web site
 *  Dr. Michael I. Mishchenko (crmim@giss.nasa.gov)
 *		http://www.giss.nasa.gov/staff/mmishchenko/t_matrix.html
 *        NASA Goddard Institute for Space Studies
 *
 *	My intent with this class is to cache the scattering properties for one particle
 *	(as defined by the attributes in abstract base class sk_NonsphericalParticle)
 *	for one orientation of this particle, with one solar zenith angle.
 *	As such, the caching descriptor for objects of this type will be:
 *		<PSD caching descriptor>_<particle type>_<aspect ratio>_<random/oriented>
 *			_<alphaStr>_<betaStr>_<solar zenith angle>_<wavelength>
 *	for example,
 *		gamma2_000001.28000_000001.35000_cyl_2.000_orient_0.000_0.000_89.000_7500.dat
 *	for a 2:1 oblate cylinder with its rotation axis oriented vertically and
 *	an 89 degree solar zenith angle.
 */
 /*-------------------------------------------------------------------------*/

class sk_TMatrixOrientedWrapper : public sk_NonsphericalParticle
{
	private:
		double								m_rat;			//!< RAT=1: part size specified by equal-volume sphere radius r_v. RAT!=1: equal-surf. area sphere radius r_s
		int									m_np;			//!< particle shape: spheroids NP=-1, cylinders NP=-2
		int									m_ndgs;			//!< no. division points over particle surface
		double								m_alpha;		//!< Euler angles (in degrees) specifying the orientation of
		double								m_beta;			//!< scattering particle relative to laboratory reference
		double								m_thetaInc;		//!< zenith angle of the incident beam in degrees
		double								m_thetaSca;		//!< zenith angle of the scattered beam in degrees
		double								m_phiInc;		//!< azimuth angle of the incident beam in degrees
		double								m_phiSca;		//!< azimuth angle of the scattered beam in degrees
		double								m_ssalb;		//!< The single scatter albedo
		nx1dArray<std::complex<double> >	m_S11;			//!< S11 amplitude matrix element array
		nx1dArray<std::complex<double> >	m_S12;			//!< S12 amplitude matrix element array
		nx1dArray<std::complex<double> >	m_S21;			//!< S21 amplitude matrix element array
		nx1dArray<std::complex<double> >	m_S22;			//!< S22 amplitude matrix element array

	private:
		sk_TMatrixOrientedWrapper&			operator=								( const sk_TMatrixOrientedWrapper& other );
		bool								Mishchenko_TMatrix						( );
		const char*							PartTypeString							( ) const;
		bool								SetDefaultAngles						( );

	protected:
		virtual bool						AllocateScatteringAngleMatrixElements	( size_t numangles );

	public:
											sk_TMatrixOrientedWrapper();
		virtual							   ~sk_TMatrixOrientedWrapper();
		bool								DeepCopy						( const sk_TMatrixOrientedWrapper& other );

	public:
		virtual ENUM_SHAPE					Shape							( ) const					{ return SHAPE_CYLINDER; }
		virtual ENUM_ALGORITHM				Algorithm						( ) const					{ return ALG_TMATORIENT; }
		virtual nxString					Get_DescriptiveParticleString	( ) const;
		virtual bool						Get_ScatteringMatrixCoeffs		( std::vector<skRTPhaseMatrix>*	coeffmatrix  );
		virtual bool						CreateClone						( sk_NonsphericalParticle** clone) const;
		virtual bool						CalculateScattering				();

		bool								Set_SizeEqualVolumeSphere		( )	{ m_rat = 1.0; return true;}
		bool								Set_SizeEqualSurfAreaSphere		( )	{ m_rat = -2.0; return true;}
		bool								Set_ShapeCylinder				();
		bool								Set_ShapeSpheroid				();
		bool								Set_ParticleOrientation			( double alpha,double beta );
		bool								Set_SunDirection				( double theta0 );

};

/*-----------------------------------------------------------------------------
 *					class sk_DiscreteDipoleWrapper	2009-6-3*/
/**  \ingroup miescat
 *	This is a wrapper class to call the Draine-Flatau Discrete Dipole Approximation
 *	FORTRAN code, which has been compiled into DISCRETEDIPOLE.dll and
 *	DISCRETEDIPOLElib. I have compiled the FORTRAN code largely "as is", only
 *	commenting some lines of screen-output print statements.
 *
 *	The Discrete Dipole Fortran code is compiled under the DiscreteDipole directory
 *	of the skopticalproperties/sources/nonspherescatter/tmatrixrandom subproject
 *	directory.
 *
 *	Reference:
 *	Draine, B.T., and Flatau, P.J. 2008, "User Guide to the Discrete Dipole
 *	Approximation Code DDSCAT 7.0", http://arxiv.org/abs/0809.0337v5 .
 *	Draine, B.T. "The Discrete-Dipole approximation and its application to
 *	interstellar graphite grains", Astrophys. J., 333, 848-872, 10/1988.
 *
 *	Also see Bruce Draine's web site
 *  Dr. Bruce T. Draine (draine@astro.princeton.edu)
 *		http://www.astro.princeton.edu/~draine/index.html
 *        Princeton University Observatory
 *
 *	For details on target orientation with respect to lab frame and associated angles
 *	and rotations, please see the user guide referenced above.
 *
 *	My intent with this class is to cache the scattering properties for one particle
 *	(as defined by the attributes in abstract base class sk_NonsphericalParticle)
 *	for one orientation of this particle, with one solar zenith angle.
 *	As such, the caching descriptor for objects of this type will be:
 *		<PSD caching descriptor>_<particle type>_<aspect ratio>_<random/oriented>
 *			_<alphaStr>_<betaStr>_<solar zenith angle>_<wavelength>
 *	for example,
 *		gamma2_000001.28000_000001.35000_cyl_2.000_orient_0.000_0.000_89.000_7500.dat
 *	for a 2:1 oblate cylinder with its rotation axis oriented vertically and
 *	an 89 degree solar zenith angle.
**/
/*---------------------------------------------------------------------------*/

class sk_DiscreteDipoleWrapper : public sk_NonsphericalParticle
{
	private:
		double							m_Theta;		//!< second Euler angle (zenith) for target rotation in lab frame
		double							m_Phi;			//!< first Euler angle (z azimuth) for target rotation in lab frame
		double							m_beta;			//!< third Euler angle (z' azimuth) for target rotation in lab frame
		bool							m_labFrame;		//!< scat. dir. defined in LF?
		int								m_dimAll;		//!< init mem allocation NX NY NZ dimensioning allowance
		double							m_eta;			//!<
		double							m_gamma;		//!<
		double							m_minRad;		//!< min radius for computation
		double							m_maxRad;		//!< max radius for computation
		double							m_deltaRad;		//!< delta rad
		int								m_numRad;		//!< number of radii
		int								m_numPlanes;	//!< number of scattering planes
		nx1dArray<double>				m_phiPlane;		//!< azimuth angle of scattering plane(s)
		double							m_scat1;		//!< min scatter angle
		double							m_scat2;		//!< max scatter angle
		double							m_dScat;		//!< scatter angle res
		double							m_theta1;		//!<
		double							m_theta2;		//!< scatter angle res
		int								m_numTheta;		//!<
		double							m_phi1;			//!<
		double							m_phi2;			//!<
		int								m_numPhi;		//!<
		double							m_beta1;		//!<
		double							m_beta2;		//!<
		int								m_numBeta;		//!<
		nxString						m_targetShape;	//!< shape string for target type
		nxString						m_databaseDir;	//!< directory containing ddscat.par, diel.tab, etc.
		nx3dArray<double>				m_f;			//!< [[4x4]x num scattering angles] scattering matrix

	private:
		sk_DiscreteDipoleWrapper&		operator=								( const sk_DiscreteDipoleWrapper& other );
		bool							Draine_DDSCAT							();
		bool							WriteParamFile							();
		bool							Set_UniformScatterAngles				( double deltaTheta );
		bool							Set_SizeDistIntegrationLimits			();


	protected:
		virtual bool					AllocateScatteringAngleMatrixElements	( size_t numangles );

	public:
										sk_DiscreteDipoleWrapper				();
		virtual						   ~sk_DiscreteDipoleWrapper				();
		bool							DeepCopy								( const sk_DiscreteDipoleWrapper& other );
		bool							Set_Angles								( double theta,double phi );

	public:
		virtual ENUM_SHAPE				Shape									( ) const			{ return SHAPE_HEX; }
		virtual ENUM_ALGORITHM			Algorithm								( )	const			{ return ALG_DDA; }
		virtual nxString				Get_DescriptiveParticleString			( ) const;
		virtual bool					CalculateScattering						( );
		virtual bool					Get_ScatteringMatrixCoeffs				( std::vector<skRTPhaseMatrix>*	coeffmatrix );
		virtual bool					CreateClone								( sk_NonsphericalParticle** clone) const;
};
