#include <skopticalproperties21.h>
#include <nxbase_threads.h>

using namespace std;

static std::mutex	g_mutexlock;
//static nxMutex 		g_mutexlock;				// A mutex to avoid calling fortran code from multiple threads


   extern "C" void  DDASCAT(	double*			QEXSUM,			// REAL*8
								double*			QABSUM,			// REAL*8
								double*			QSCSUM,			// REAL*8
								int*			NSCAT,			// INT
								double*			F				// REAL*8
								);		


/*---------------------------------------------------------------------------
 *'			sk_DiscreteDipoleWrapper::sk_DiscreteDipoleWrapper		2009-6-2
 *-------------------------------------------------------------------------*/

sk_DiscreteDipoleWrapper::sk_DiscreteDipoleWrapper()
{
	Set_Radius( 0.5 );			// 
	Set_Delta ( 1.00e-4 );		// override while testing this out
	m_deltaRad          = 0.0;
	m_Theta				= 50;
	m_Phi				= 5;
	m_beta				= 30;
	m_gamma				= 1e-2;
	m_dimAll			= 125;			
	m_eta				= 0.25;
	m_scat1				= 0;		
	m_scat2				= 180;		
	m_dScat				= 0.3333333;
	m_maxRad            = Get_Radius();
	m_minRad			= Get_Radius();
	m_theta1            = m_Theta;
	m_theta2	        = m_Theta;
	m_phi1              = m_Phi;
	m_phi2		        = m_Phi;
	m_beta1             = m_beta;
	m_beta2	            = m_beta;
	m_numRad			= 1;
	m_numBeta           = 1;
	m_numTheta	        = 1;
	m_numPhi            = 1;
	m_labFrame			= true;
	m_numPlanes			= 1;
	m_phiPlane.SetSize( m_numPlanes );
	m_phiPlane.At(0)	= 45;
//	m_phiPlane.At(1)	= 135;
	m_targetShape		= "\'HEX_PRISM\'";
	m_databaseDir		= "e:/software/ddascat/doc";
	Set_UniformScatterAngles( m_dScat );
}

/*---------------------------------------------------------------------------
 *'			sk_DiscreteDipoleWrapper::~sk_DiscreteDipoleWrapper		2009-6-2
 *-------------------------------------------------------------------------*/

sk_DiscreteDipoleWrapper::~sk_DiscreteDipoleWrapper()
{
	m_phiPlane.erase();
	m_f.erase();
}


/*-----------------------------------------------------------------------------
 *					sk_DiscreteDipoleWrapper::operator=			2009-6-2	*/
/*---------------------------------------------------------------------------*/

bool sk_DiscreteDipoleWrapper::DeepCopy( const sk_DiscreteDipoleWrapper& other )
{
	bool	ok;

	ok          = sk_NonsphericalParticle::DeepCopy( other );
	m_Theta		= other.m_Theta;
	m_Phi		= other.m_Phi;
	m_beta		= other.m_beta;
	m_gamma		= other.m_gamma;
	m_labFrame	= other.m_labFrame;		
	m_dimAll	= other.m_dimAll;		
	m_eta		= other.m_eta;			
	m_minRad	= other.m_minRad;		
	m_maxRad	= other.m_maxRad;		
	m_deltaRad	= other.m_deltaRad;		
	m_numRad	= other.m_numRad;		
	m_numPlanes	= other.m_numPlanes;	
	m_scat1		= other.m_scat1;		
	m_scat2		= other.m_scat2;		
	m_dScat		= other.m_dScat;		
	m_theta1	= other.m_theta1;		
	m_theta2	= other.m_theta2;		
	m_numTheta	= other.m_numTheta;		
	m_phi1		= other.m_phi1;			
	m_phi2		= other.m_phi2;			
	m_numPhi	= other.m_numPhi;		
	m_beta1		= other.m_beta1;		
	m_beta2		= other.m_beta2;		
	m_numBeta	= other.m_numBeta;		
	ok			= ok && m_phiPlane.DeepCopy( other.m_phiPlane );		
	ok			= ok && m_f.DeepCopy( other.m_f );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"sk_DiscreteDipoleWrapper::DeepCopy, Error copying contents, this object is in an ill-defined state and may produce strange results");
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *'				sk_DiscreteDipoleWrapper::CreateClone				2009-6-3
 * 
 *-------------------------------------------------------------------------*/
bool sk_DiscreteDipoleWrapper::CreateClone(sk_NonsphericalParticle** userclone) const
{
	sk_DiscreteDipoleWrapper*	clone;
	bool						ok;

	clone = new sk_DiscreteDipoleWrapper;
	ok = (clone != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"sk_DiscreteDipoleWrapper::CreateClone, Error creating clone object");
	}
	else
	{
		clone->AddRef();
		ok = clone->DeepCopy(*this);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"sk_DiscreteDipoleWrapper::CreateClone, Error performing deep Copy, I'm goin to destroy the clone");
			clone->Release();
			clone = NULL;
		}
	}
	*userclone = clone;
	return ok;
}

/*---------------------------------------------------------------------------
 *'		sk_DiscreteDipoleWrapper::Set_Radius				2009-6-25
 * 
 *-------------------------------------------------------------------------*/
bool sk_DiscreteDipoleWrapper::Set_SizeDistIntegrationLimits()
{
	// maybe update this later to run multiple radii within single DLL call
	m_maxRad = Get_Radius();
	m_minRad = Get_Radius();
	return true;
}


/*---------------------------------------------------------------------------
 *'		sk_DiscreteDipoleWrapper::Get_DescriptiveParticleString		2009-6-2
 * Return a string that defines a particle according to its parameters, which
 * is used for caching the computed particle's scattering properties
 *-------------------------------------------------------------------------*/
nxString sk_DiscreteDipoleWrapper::Get_DescriptiveParticleString() const
{
	nxString	basename;

	basename.sprintf("hex_%3.2f_orient_%5.3f_%5.3f", (double)Get_AspectRatio(), (double)m_Theta, (double)m_Phi );
	return ( basename ); 
}

/*---------------------------------------------------------------------------
 *'			sk_DiscreteDipoleWrapper::Set_Angles				2009-6-18
 *	Set the direction of the sun: input solar zenith angle.
 *	
 *-------------------------------------------------------------------------*/
bool sk_DiscreteDipoleWrapper::Set_Angles( double theta,double phi )
{
	SetDirty( m_Theta != theta );
	SetDirty( m_Phi != phi );
	m_Theta		= theta;
	m_Phi		= phi;
	return true;
}


/*---------------------------------------------------------------------------
 *'		sk_DiscreteDipoleWrapper::Set_UniformScatterAngles		2009-6-18
 *
 *-------------------------------------------------------------------------*/
bool sk_DiscreteDipoleWrapper::Set_UniformScatterAngles( double deltaTheta )
{
	double	startangle[1] = {  0.0};
	double  endangle[1]   = { 90.0};
	double  resolution[1] = {deltaTheta};
    bool ok;

	ok =  Set_SymmetricScatterAngles( startangle, endangle, resolution, 1 );
	return ok;
}

/*---------------------------------------------------------------------------
 *'		sk_DiscreteDipoleWrapper::AllocateScatteringMatrixElements		2009-6-16
 *
 *-------------------------------------------------------------------------*/
bool sk_DiscreteDipoleWrapper::AllocateScatteringAngleMatrixElements( size_t numangles )
{
	bool	ok;

	ok = m_f.SetSize( 4,4,numangles*m_numPlanes );		
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "sk_DiscreteDipoleWrapper::AllocateScatteringMatrixElements, could not allocate scattering matrix elements");
		m_f.erase();
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *'					sk_DiscreteDipoleWrapper::CalculateScattering		2009-6-2
 *-------------------------------------------------------------------------*/

bool sk_DiscreteDipoleWrapper::CalculateScattering()
{
	bool	ok;

	ok = !IsDirty();
	if (!ok)
	{
		ok = (Get_NumAngles() > 0);
		ok = ok && Set_SizeDistIntegrationLimits ();
		ok = ok && Draine_DDSCAT();
		if (ok) ClearDirty();
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *'					sk_DiscreteDipoleWrapper::Draine_DDSCAT			2009-6-2
 *
 *-------------------------------------------------------------------------*/

bool sk_DiscreteDipoleWrapper::Draine_DDSCAT()
{
	bool				ok;
	int					numscatang;
	nx1dArray<double>	qexsum(2);
	nx1dArray<double>	qabsum(2);
	nx1dArray<double>	qscsum(2);
	double				qext;
	double				qsca;
	double				area;
	double				cext;
	double				csca;
	double				radius;

	{
	std::unique_lock<std::mutex>	lock(g_mutexlock);
//	g_mutexlock.Lock();																// grab the mutex, Not sure if Wiscombe is thread safe

	numscatang  = (int)Get_NumAngles();
	radius      = Get_Radius();

	ok			= WriteParamFile();			// write the configuration file
	DDASCAT(			qexsum.UnsafeArrayBasePtr(),
						qabsum.UnsafeArrayBasePtr(),
						qscsum.UnsafeArrayBasePtr(),
						&numscatang,
						m_f.UnsafeArrayBasePtr()
						);
     }
//	g_mutexlock.Unlock();														// Release the mutex.

	qext	= 0.5*(qexsum.At(0)+qexsum.At(1));	// mean cross sections over incident polarizations
	qsca	= 0.5*(qscsum.At(0)+qscsum.At(1));

	area    = nxmath::Pi*radius*radius;
	cext    = qext*area;					// extinction cross-section (in units of microns**2)
	csca    = qsca*area;					// scattering cross-section (in units of microns**2)		

	ok		= ok && SetCrossSectionsFromFortran( cext, csca);
	return ok;
}


/*---------------------------------------------------------------------------
 *'			sk_DiscreteDipoleWrapper::Get_ScatteringMatrixCeoffs		2009-6-2
 *	Get the scattering matrix elements (as a function of scattering angle) for 
 *	this nonspherical scattering particle.
 *	Scattering matrix elements (Bohren and Huffman, 1983, eqn 3.16; a.k.a. 
 *	transformation matrix elements: Hansen and Travis, 1974, eqn 2.36; Liou, 
 *	2002, 5.2.105) are set in this routine. These are the following:
 *
 *	 F = 	( F11	F12		F13		F14	)
 *			( F21	F22		F23		F24	)
 *			( F31	F32		F33		F34	)
 *			( F41	F42		F43		F44	)
 *	 (at each value of cos(\Theta) = m_xmu)
 *	Note that the normalization
 *		F(\Theta)*(4*pi/k^2/Csca) = P(\Theta)
 *	to obtain the phase matrix is done within the calling routine of the 
 *	extinction class, skOpticalProperties_Ice::CalculateCrossSections.
 *
 *	The arrays are set up as follows:
 *	coeffmatrix(0,0,*) = F11		*= all angles in m_xmu (see Get_CosAngles())
 *	coeffmatrix(1,0,*) = F21		*= all angles in m_xmu
 *	coeffmatrix(0,1,*) = F12		*= all angles in m_xmu
 *  ...
 *  coeffmatrix(3,2,*) = F43		*= all angles in m_xmu
 *  coeffmatrix(3,3,*) = F44		*= all angles in m_xmu,
 *-------------------------------------------------------------------------*/

bool sk_DiscreteDipoleWrapper::Get_ScatteringMatrixCoeffs( std::vector<skRTPhaseMatrix>*	coeffmatrix )
{
	bool								ok;
	size_t								i;
	size_t								numangles;

	numangles = Get_NumAngles();
	ok = ( numangles > 0);
	if (ok)
	{
		CalculateScattering();								// Update Cache so output parameters reflect latest input parameters
//	c   = (4.0*nxmath::Pi)/(m_k*m_k*m_Csca);			// Get the normalization constant reciprocal of equation 2.40 in Hansen And Travis

		try
		{

			coeffmatrix->resize(numangles);
			NXTRACE_ONCEONLY( firsta,("***** CHECK ***** sk_DiscreteDipoleWrapper::Get_ScatteringMatrixCoeffs, Truitt needs to confirm this is the correct way to unpack the m_f matrix\n"));
			NXTRACE_ONCEONLY( firstb,("***** CHECK ***** sk_DiscreteDipoleWrapper::Get_ScatteringMatrixCoeffs, Still need to ensure the scalar phase matrix integrates to 4 pi over the unit sphere\n"));

			for (i=0; i < numangles; i++)
			{
				(*coeffmatrix)[i].At(1,1) = m_f.At( 0,0,i ); 
				(*coeffmatrix)[i].At(1,2) = m_f.At( 0,1,i ); 
				(*coeffmatrix)[i].At(1,3) = m_f.At( 0,2,i );
				(*coeffmatrix)[i].At(1,4) = m_f.At( 0,3,i );

				(*coeffmatrix)[i].At(2,1) = m_f.At( 1,0,i );
				(*coeffmatrix)[i].At(2,2) = m_f.At( 1,1,i );
				(*coeffmatrix)[i].At(2,3) = m_f.At( 1,2,i );
				(*coeffmatrix)[i].At(2,4) = m_f.At( 1,3,i );

				(*coeffmatrix)[i].At(3,1) = m_f.At( 2,0,i );
				(*coeffmatrix)[i].At(3,2) = m_f.At( 2,1,i );
				(*coeffmatrix)[i].At(3,3) = m_f.At( 2,2,i );
				(*coeffmatrix)[i].At(3,4) = m_f.At( 2,3,i );

				(*coeffmatrix)[i].At(4,1) = m_f.At( 3,0,i );
				(*coeffmatrix)[i].At(4,2) = m_f.At( 3,1,i );
				(*coeffmatrix)[i].At(4,3) = m_f.At( 3,2,i );
				(*coeffmatrix)[i].At(4,4) = m_f.At( 3,3,i );
			}
		}
		catch (...)
		{
			nxLog::Record(NXLOG_WARNING,"sk_DiscreteDipoleWrapper::Get_ScatteringMatrixCoeffs, error allocating memory for phase matrices.");
			ok = false;
		}
	}
	if (!ok) coeffmatrix->clear();
	return ok;
}

/*---------------------------------------------------------------------------
 *'				sk_DiscreteDipoleWrapper::WriteParamFile			2009-6-12
 *	Write the configuration file 'ddscat.par' for driving the DDSCAT algorithm.
 *	This is performed just before the call to Draine_DDSCAT, within the mutex 
 *	lock.
 *
 *  The accuracy of computations is limited by the number of dipoles used to 
 *	represent the scattering target.  For reliable computation of scattering
 *	matrix elements, it is required that
 *			|m|k_0d <= 0.3,	
 *	where d is the interdipole spacing.
 *	I have done a fit of the (aspect ratio 2) column and plate dipole numbers
 *	that will ensure this condition.
 *	As well, the FFT algorithm used works best for dipole numbers that are
 *	of the form 2^i*3^j*5^k.  I match the across-face and along-axes dipole
 *	numbers to these values.
 *  //TODO: The location of the necessary files (ddscat.par,diel.tab) is quite
 *			flaky and error-prone: fix these.
 *-------------------------------------------------------------------------*/

bool sk_DiscreteDipoleWrapper::WriteParamFile()
{
	bool						ok;
	nx1dArray<int>				defaultdipolenumbers;
	std::vector<int>			nums;
	nxString					filename,dielfilename,frameStr,t,dipolenumsfilename;
	nxFile						f;
	int							fundNum,longNum,numAxis,numFace,orientation;
	double						aspectratio;

	filename.sprintf( "%s/ddscat.par",(const char*)m_databaseDir );
	dielfilename.sprintf( "\'%s/diel.tab\'",(const char*)m_databaseDir );
	dipolenumsfilename.sprintf( "%s/defaultdipolenumbers.txt",(const char*)m_databaseDir );
	defaultdipolenumbers.InputColumnMajorText( dipolenumsfilename );

	nums        = defaultdipolenumbers.STLVector();
	aspectratio = Get_AspectRatio();

	if (m_numRad>1)	m_deltaRad = (m_maxRad-m_minRad)/(m_numRad-1);
	else			m_deltaRad = 0;

	if (aspectratio > 1.0)		
	{
		fundNum		= 75 + int((30.0/1.12)*(m_maxRad - 2.24));
		longNum		= int(aspectratio*(double)fundNum);
		numAxis		= nums.at( lower_bound( nums.begin(),nums.end(), (double)fundNum ) - nums.begin() );
		numFace		= nums.at( lower_bound( nums.begin(),nums.end(), (double)longNum ) - nums.begin() );
		orientation = 3;
	}
	else
	{
		fundNum		= 90 + int((30.0/1.12)*(m_maxRad - 2.24));
		longNum		= int(aspectratio*(double)fundNum);
		numAxis		= nums.at( lower_bound( nums.begin(),nums.end(), (double)longNum ) - nums.begin() );
		numFace		= nums.at( lower_bound( nums.begin(),nums.end(), (double)fundNum ) - nums.begin() );
		orientation = 1;
	}
	if (m_labFrame)		frameStr.sprintf( "\'LFRAME\'" );
	else				frameStr.sprintf( "\'TFRAME\'" );

	f.Open( filename, "w");
	ok = f.IsOpen();
	if (ok)
	{
		t.sprintf( "\' ========= Parameter file for v7.0.7 ===================\'\n" );				f.Write( t,t.GetLength(),1 );
		t.sprintf( "\'**** Preliminaries ****\'\n" );												f.Write( t,t.GetLength(),1 );
		t.sprintf( "\'NOTORQ\' = CMTORQ*6 (NOTORQ, DOTORQ) -- either do or skip torque calculations\n" );	f.Write( t,t.GetLength(),1 );
		t.sprintf( "\'PBCGS2\' = CMDSOL*6 (PBCGS2, PBCGST, PETRKP) -- select solution method\n" );	f.Write( t,t.GetLength(),1 );
		t.sprintf( "\'GPFAFT\' = CMDFFT*6 (GPFAFT, FFTMKL)\n" );									f.Write( t,t.GetLength(),1 );
		t.sprintf( "\'GKDLDR\' = CALPHA*6 (GKDLDR, LATTDR)\n" );									f.Write( t,t.GetLength(),1 );
		t.sprintf( "\'NOTBIN\' = CBINFLAG (NOTBIN, ORIBIN, ALLBIN)\n" );							f.Write( t,t.GetLength(),1 );
		t.sprintf( "\'**** Initial Memory Allocation ****\'\n" );									f.Write( t,t.GetLength(),1 );
		t.sprintf( "%d  %d  %d  = dimensioning allowance for target generation\n",				
							(int)m_dimAll,(int)m_dimAll,(int)m_dimAll );							f.Write( t,t.GetLength(),1 );
		t.sprintf( "\'**** Target Geometry and Composition ****\'\n" );								f.Write( t,t.GetLength(),1 );
		t.sprintf( "%s = CSHAPE*9 shape directive\n",
							(const char*)m_targetShape );														f.Write( t,t.GetLength(),1 );
		t.sprintf( "%d  %d  %d  = shape parameters 1 - 3\n",	
							(int)numAxis,(int)numFace,(int)orientation );							f.Write( t,t.GetLength(),1 );
		t.sprintf( "1        = NCOMP = number of dielectric materials\n" );							f.Write( t,t.GetLength(),1 );
		t.sprintf( "%s = file with refractive index 1\n",
							(const char*)dielfilename );															f.Write( t,t.GetLength(),1 );
		t.sprintf( "\'**** Error Tolerance ****\'\n" );												f.Write( t,t.GetLength(),1 );
		t.sprintf( "%e = TOL = MAX ALLOWED (NORM OF |G>=AC|E>-ACA|X>)/(NORM OF AC|E>))\n",
							(double)Get_Delta() );														f.Write( t,t.GetLength(),1 );
		t.sprintf( "\'**** Interaction cutoff parameter for PBC calculations ****\'\n" );			f.Write( t,t.GetLength(),1 );
		t.sprintf( "%e =  GAMMA (1e-2 is normal, 3e-3 for greater accuracy)\n",(double)m_gamma );	f.Write( t,t.GetLength(),1 );
		t.sprintf( "\'**** Angular resolution for calculation of <cos>, etc. ****\'\n" );			f.Write( t,t.GetLength(),1 );
		t.sprintf( "%.3f = ETASCA (number of angles is proportional to [(3+x)/ETASCA]^2 )\n",
							(double)m_eta );														f.Write( t,t.GetLength(),1 );
		t.sprintf( "\'**** Wavelengths (micron) ****\'\n" );										f.Write( t,t.GetLength(),1 );
		t.sprintf( "%.3f %.3f 1 \'LIN\' = wavelengths (first,last,how many,how=LIN,INV,LOG)\n",
							(double)Get_Lambda(),(double)Get_Lambda() );							f.Write( t,t.GetLength(),1 );
		t.sprintf( "\'**** Effective Radii (micron) ****\'\n" );									f.Write( t,t.GetLength(),1 );
		t.sprintf( "%.3f %.3f %d \'LIN\' = aeff (first,last,how many,how=LIN,INV,LOG)\n",
							(double)m_minRad,(double)m_maxRad,(int)m_numRad );						f.Write( t,t.GetLength(),1 );
		t.sprintf( "\'**** Define Incident Polarizations ****\'\n" );								f.Write( t,t.GetLength(),1 );
		t.sprintf( "(0.,0) (1.,0.) (0.,0.) = Polarization state e01 (k along x axis)\n" );			f.Write( t,t.GetLength(),1 );
		t.sprintf( "2 = IORTH  (=1 to do only pol. state e01; =2 to also do orth. pol. state)\n" );	f.Write( t,t.GetLength(),1 );
		t.sprintf( "\'**** Specify which output files to write ****\'\n" );							f.Write( t,t.GetLength(),1 );
		t.sprintf( "0 = IWRKSC (=0 to suppress, =1 to write \".sca\" file for each target orient).\n" );f.Write( t,t.GetLength(),1 );
		t.sprintf( "0 = IWRPOL (=0 to suppress, =1 to write \".pol\" file for each (BETA,THETA))\n" );	f.Write( t,t.GetLength(),1 );
		t.sprintf( "\'**** Prescribe Target Rotations ****\'\n" );									f.Write( t,t.GetLength(),1 );
		t.sprintf( "%.3f %.3f %d = BETAMI, BETAMX, NBETA  (beta=rotation around a1)\n",	
							(double)m_beta1,(double)m_beta2,(int)m_numBeta );						f.Write( t,t.GetLength(),1 );
		t.sprintf( "%.3f %.3f %d = THETMI, THETMX, NTHETA (theta=angle between a1 and k)\n",
							(double)m_theta1,(double)m_theta2,(int)m_numTheta );					f.Write( t,t.GetLength(),1 );
		t.sprintf( "%.3f %.3f %d = PHIMIN, PHIMAX, NPHI (phi=rotation angle of a1 around k)\n",
							(double)m_phi1,(double)m_phi2,(int)m_numPhi );							f.Write( t,t.GetLength(),1 );
		t.sprintf( "\'**** Specify first IWAV, IRAD, IORI (normally 0 0 0) ****\'\n" );				f.Write( t,t.GetLength(),1 );
		t.sprintf( "0  0  0    = first IWAV, first IRAD, first IORI (0 0 0 to begin fresh)\n" );	f.Write( t,t.GetLength(),1 );
		t.sprintf( "\'**** Select Elements of S_ij Matrix to Print ****\'\n" );						f.Write( t,t.GetLength(),1 );
		t.sprintf( "9	= NSMELTS = number of elements of S_ij to print (not more than 9)\n" );		f.Write( t,t.GetLength(),1 );
		t.sprintf( "11 22 33 44 12 13 21 34 43	= indices ij of elements to print\n" );				f.Write( t,t.GetLength(),1 );
		t.sprintf( "\'**** Specify Scattered Directions ****\'\n" );								f.Write( t,t.GetLength(),1 );
		t.sprintf( "%s = CMDFRM (LFRAME, TFRAME for Lab Frame or Target Frame)\n",
							(const char*)frameStr );																f.Write( t,t.GetLength(),1 );
		t.sprintf( "%d = NPLANES = number of scattering planes\n",
							(int)m_numPlanes );														f.Write( t,t.GetLength(),1 );
		for(int i=0;i<m_numPlanes;i++)
		{
			t.sprintf( "%.3f %.3f %.3f %.3f = phi, thetan_min, thetan_max, dtheta (in deg) for plane %d\n",
				(double)m_phiPlane.At(i),(double)m_scat1,(double)m_scat2,(double)m_dScat,(int)i );	f.Write( t,t.GetLength(),1 );
		}
		f.Close();
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"sk_DiscreteDipoleWrapper::WriteParamFile, Error writing ice crystal cache to file <%s>", (const char*)filename );
	}
	return ok;
}



//TODO: make an enumeration-type definition for the (following) target types possible in DDASCAT
//static enum ENUM_DDA_SHAPE{	'FROM_FILE',	// read shape and composition data from file, for dielectric tensors that are diagonal in Target Frame
//							'ANIFRMFIL',	// read shape and composition data from file, for general anisotropic dielectric tensors with orientation angles THETADF,PHIDF,BETAD relative to Target Frame.
//							'ANIELLIPS',	// ellipsoid of anisotropic material
//							'ANIRCTNGL',	// homogeneous anisotropic rectangular target
//							'ANI_ELL_2',	// two touching anisotropic ellipsoids of materials 1-6
//							'ANI_ELL_3',	// three touching anisotropic ellipsoids of materials 1-9
//							'BISLINPBC',	// bilayer slab with periodic grid of lines parallel to z on top, with y-period/d=PYD [or, if PYD=0, a single line]
//							'CONELLIPS',	// two concentric ellipsoids of materials 1,2
//							'CYLINDER1',	// homogeneous finite cylinder
//							'CYLNDRCAP',	// homogeneous cylinder with hemispherical endcaps
//							'CYLNDRPBC',	// 1-d or 2-d array of finite cylinders
//							'DSKBLYPBC',	// 1-d or 2-d array of disk on bilayer rect. slab
//							'DSKRCTNGL',	// single disk on rectangular slab
//							'DSKRCTPBC',	// 1-d or 2-d array of disk on rectangular slab
//							'DW1996TAR',	// 13-cube target used by Draine & Weingartner 1996
//							'ELLIPSOID',	// ellipsoid (homogeneous and isotropic)
//							'ELLIPSPBC',	// 1-d or 2-d array of ellipsoids
//							'ELLIPSO_2',	// two touching isotropic ellipsoids of materials 1 and 2
//							'ELLIPSO_3',	// three touching isotropic ellipsoids of materials 1,2,3
//							'GAUSS_SPH',	// gaussian sphere target
//							'HEXGONPBC',	// 1-d or 2-d array of hexagonal prisms
//							'HEX_PRISM',	// homogeneous hexagonal prism
//							'LYRD_SLAB',	// layered slab target, with up to 4 separate material layers
//							'LYRSLBPBC',	// 1-d or 2-d array of layered rect. slab targets, with up to 4 material layers
//							'MLTBLOCKS',	// collection of cubic blocks defined by data in file 'blocks.par'
//							'RCTGLBLK3',	// isolated target: 3 rectangular blocks with centers on x-axis
//							'RCTGLPRSM',	// homogeneous rectangular prism
//							'RCTGL_PBC',	// 1-d or 2-d array of rectangular targets
//							'SPHERES_N',	// multisphere target union of N spheres
//							'SPHRN_PBC',	// 1-d or 2-d array of multisphere target
//							'SPHROID_2',	// two touching spheroids with symmetry axes at specified angle!
//							'SPH_ANI_N',	// multisphere target, with spheres that can have different, anisotropic, composition
//							'TETRAHDRN',	// regular tetrahedron
//							'TRILYRPBC',	// periodic target: 3 layer rectangular structure)
//							'TRNGLPRSM',	// triangular prism (homogeneous and isotropic)
//							'UNIAXICYL'		// cylinder of unixaxial material
//						};
