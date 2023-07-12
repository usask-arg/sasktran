
// Forward declaration to disable implicit conversion
class SKTRAN_ScatMat_MIMSNC;
class SKTRAN_PhaseMat_MIMSNC; // The trailing MIMSNC could probably be removed from this name
class SKTRAN_ScatMat_Rot;


/*-----------------------------------------------------------------------------
 *					SKTRAN_Stokes_RONC		 2016- 7- 12*/
/** @ingroup stokes 
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_Stokes_NC
{
	private: 
		SKTRAN_StokesScalar m_I, m_Q, m_U;

	public: 
		//SKTRAN_Stokes_RONC(){ SetTo( 0.0 );}
		SKTRAN_StokesScalar I() const { return m_I;}
		SKTRAN_StokesScalar Q() const { return m_Q;}
		SKTRAN_StokesScalar U() const { return m_U;}
		SKTRAN_StokesScalar V() const { return 0.0;}
		void Assign_I ( SKTRAN_StokesScalar val ) { m_I  = val;}
		void Assign_Q ( SKTRAN_StokesScalar val ) { m_Q  = val;}
		void Assign_U ( SKTRAN_StokesScalar val ) { m_U  = val;}
		void Assign_V ( SKTRAN_StokesScalar val ) { ;}
		void Add_I    ( SKTRAN_StokesScalar val ) { m_I += val;}
		void Add_Q    ( SKTRAN_StokesScalar val ) { m_Q += val;}
		void Add_U    ( SKTRAN_StokesScalar val ) { m_U += val;}
		void Add_V    ( SKTRAN_StokesScalar val ) { ;}

		void SetTo ( double d ){ m_I=m_Q=m_U=d;}
		void SetTo ( double I, double Q, double U, double V ) { m_I=I; m_Q=Q; m_U=U;}
		//SKTRAN_Stokes_RONC& operator=  ( const SKTRAN_Stokes_RONC &other ) { m_I =other.I(); m_Q =other.Q(); m_U =other.U(); return *this;}
        SKTRAN_Stokes_NC& operator+= ( const SKTRAN_Stokes_NC &other ) { m_I+=other.m_I; m_Q+=other.m_Q; m_U+=other.m_U; return *this;}
        SKTRAN_Stokes_NC& operator*= ( double val ) { m_I*=val; m_Q*=val; m_U*=val; return *this;}
        SKTRAN_Stokes_NC  operator + ( const SKTRAN_Stokes_NC& other ) const { SKTRAN_Stokes_NC result(*this); result+=other; return result;}
        SKTRAN_Stokes_NC  operator * ( double val )                    const { SKTRAN_Stokes_NC result(*this); result*=val;   return result;}
		void LMultBy (const SKTRAN_ScatMat_MIMSNC& s);
        void Normalize();
		void                RotatePolarPlaneThru( double cosEta, double sinEta )                       //!< Rotate the basis for the electric field components used in this vector through angle eta (Tony's thesis page 44)
        { 
            // inlined for HR polarized mode
            SKRTFLOAT cte, ste, a, b;

            cte = cosEta*cosEta - sinEta*sinEta;
            ste = 2.0*cosEta*sinEta;

            a = cte*m_Q - ste*m_U;
            b = ste*m_Q + cte*m_U;

            m_Q = a;
            m_U = b;
        }
};

SKTRAN_Stokes_NC operator*(const skRTPhaseMatrix&p, const SKTRAN_Stokes_NC& v);


/*-------------------------------------------------------------------------- 
 *					class SKTRAN_ScatMat_MIMSNC					  */
/**	@ingroup stokes
 *	A lighter version of the skRTPhaseMatrix class. Assumes the species
 *  of interest is (M)acroscopically isotropic and (M)irror-(S)ymmetric, 
 *  and that the interaction of  (C)ircular polarization with the particle 
 *  is (N)egligible. That is, it is assumed the scattering matrix can 
 *  be written, 
 *  \code
 *  F(\Theta) = [  p11  p21    0    0 ;
 *                 p21  p22    0    0 ;
 *                   0    0  p33    0 ;
 *                   0    0    0    0  ]
 * \endcode
 *
 *	\par Storage and Layout
 *  The class uses local storage for the 4 unique array elements of #SKRTFLOAT.
 *  See documentation of skRTPhaseMatrix for memory considerations and indexing
 *  details. 
 *
 *	\par Scattering Matrix Normalization
 *	We use the same normalization as Hansen and Travis: the phase
 *	matrix when integrafted over \f$4\pi\f$ solid angle gives \f$4\pi\f$.
 *
 *	\par Iteration and access
 *	This class currently DOES NOT provide the iterators supplied
 *  by skRTPhaseMatrix, nor does it return references to matrix elements. 
 *
 */
/*-------------------------------------------------------------------------*/

class SKTRAN_ScatMat_MIMSNC
{
	private:
		SKRTFLOAT				m_p11;
		SKRTFLOAT               m_p21;
		SKRTFLOAT               m_p22;
		SKRTFLOAT               m_p33;
		
		//SKTRAN_ScatMat_RONC( const SKTRAN_PhaseMat_MIMSNC& other ) /*=delete*/; // Explicitly deleted copy constructors aren't supported in VS2012; uncomment when we move to 2015
		//SKTRAN_ScatMat_RONC( const SKTRAN_ScatMat_Rot& other ) /*=delete*/; // Explicitly deleted copy constructors aren't supported in VS2012; uncomment when we move to 2015
		
	public:
        SKTRAN_ScatMat_MIMSNC          ( );                                         /*!< Create a blank instance */
        explicit SKTRAN_ScatMat_MIMSNC ( const skRTPhaseMatrix& p );
        SKTRAN_ScatMat_MIMSNC          ( const SKTRAN_ScatMat_MIMSNC& other );  /*!< Copy the other phase matrix to this instance */
		
        void                  AssignAt            ( int row, int col, SKRTFLOAT val );              /*!< Assign <em>val</em> to element <em>(row,col)</em> using a 1 based nomenclature*/
		SKRTFLOAT             At                  ( int row, int col ) const;                       /*!< Returns a reference to the element index by <em>(row,col)</em> using a 1 based nomenclature*/
        SKRTFLOAT             p11                 ( ) const { return m_p11; }
		void                  SetTo               ( double val );                                   /*!< Set all of the elements of the phase matix to \e val */
		void                  SetTo               ( float val );                                    /*!< Set all of the elements of the phase matix to \e val */
		void                  SetToIdentity       ( );                                          /*!< This becomes an identity matrix, with (4,4) still assumed to be zero. */
		//SKTRAN_ScatMat_RONC&  RMultBy             ( const SKTRAN_ScatMat_RONC& other );       /*!< This becomes the result of right multiplying itself by #other */
		//SKTRAN_ScatMat_RONC&  LMultBy             ( const SKTRAN_ScatMat_RONC& other );       /*!< This becomes the result of left multiplying itself by #other */
		void                  AddToThis           ( const SKTRAN_ScatMat_MIMSNC& other, double weight );
        SKTRAN_ScatMat_MIMSNC&  operator +=         ( const SKTRAN_ScatMat_MIMSNC& other);        /*!< Add the other phase matrix to this one */
        SKTRAN_ScatMat_MIMSNC&  operator *=         ( double value);                                  /*!< Multiply the phase matrix by \e val */
        SKTRAN_ScatMat_MIMSNC&  operator *=         ( float  value);                                  /*!< Multiply the phase matrix by \e val */
        SKTRAN_ScatMat_MIMSNC&  operator =          ( const SKTRAN_ScatMat_MIMSNC& other);        /*!< Assign (copy) the other phase matrix to this one */
		SKTRAN_Stokes_NC        operator *          ( const SKTRAN_Stokes_NC& stokes )     const;     /*!< Multiply this phase matrix by a Stokes Vector to get a Stokes vector */
        SKTRAN_ScatMat_MIMSNC   operator *          ( double value)                        const;     /*!< Multiply this phase matrix by a constant and return the answer */
        SKTRAN_ScatMat_MIMSNC   operator *          ( float  value)                        const;
        SKTRAN_ScatMat_MIMSNC   operator -          ( const SKTRAN_ScatMat_MIMSNC& other)   const;
        SKTRAN_ScatMat_MIMSNC   operator +          ( const SKTRAN_ScatMat_MIMSNC& other)   const;
		
		void LApplyTo  ( SKTRAN_Stokes_NC* v ) const;
        void Normalize ( );

		friend class SKTRAN_PhaseMat_MIMSNC;
};


/*--------------------------------------------------------------------------
 *					class SKTRAN_ScatMat_Rot					  */
/**	@ingroup stokes
 *  Stores a matrix to rotate the coordinate system in which a stokes 
 *  vector is defined. The user is responsible for making sure inputs
 *  are sensible -- no error checking is performed. The matrix stored 
 *  is equivalent to:
 *
 *	R.SetTo(0.0);
 *  R.At(1,1) =  1.0;
 *  R.At(2,2) =  cosTwoEta;
 *  R.At(2,3) = -sinTwoEta;
 *  R.At(3,2) =  sinTwoEta;
 *  R.At(3,3) =  cosTwoEta;
 *  R.At(4,4) =  1.0;
	(See Mishchenko 2002, ) 
 */
/*-------------------------------------------------------------------------*/
class SKTRAN_ScatMat_Rot
{
	private:
		double m_cos2eta;
		double m_sin2eta;

	public:
		SKTRAN_ScatMat_Rot         ( );
		SKTRAN_ScatMat_Rot         ( double cosEta, double sinEta );
		void             SetAngle  ( double cosEta, double sinEta );
		SKTRAN_Stokes_NC operator* ( const SKTRAN_Stokes_NC& v ) const;
		double           At        ( int row, int col ) const;

		friend class SKTRAN_PhaseMat_MIMSNC;
};

/*--------------------------------------------------------------------------
 *					class SKTRAN_PhaseMat_MIMSNC					  */
/**	@ingroup stokes
 *  Similar to SKTRAN_ScatMat_RONC, but spans the space phase matrices
 *  for which P_X4==P_4X==0. This allows the class to represent any product
 *  of matrices of the form given by SKTRAN_ScatMat_RONC and rotation 
 *  matrices. 
 *
 */
/*-------------------------------------------------------------------------*/
class SKTRAN_PhaseMat_MIMSNC // The trailing MIMSNC could probably be removed from this name
{
	private:
		std::array<SKRTFLOAT,3*3>	m_p; // Storage format is ROW-MAJOR (C++ style; skRTPhaseMatrix is column-major)
		
	public:
                              SKTRAN_PhaseMat_MIMSNC ( );                                         /*!< Create a blank instance */
                              explicit SKTRAN_PhaseMat_MIMSNC ( const skRTPhaseMatrix& p );
                              SKTRAN_PhaseMat_MIMSNC ( const SKTRAN_ScatMat_Rot&  other );  /*!< Copy the other phase matrix to this instance */
                              SKTRAN_PhaseMat_MIMSNC ( const SKTRAN_ScatMat_MIMSNC& other );  /*!< Copy the other phase matri to this instance */
                              SKTRAN_PhaseMat_MIMSNC ( const SKTRAN_PhaseMat_MIMSNC& other );  /*!< Copy the other phase matri to this instance */
		void                  AssignAt            ( int row, int col, SKRTFLOAT val );              /*!< Assign <em>val</em> to element <em>(row,col)</em> using a 1 based nomenclature*/
		SKRTFLOAT             At                  ( int row, int col ) const;                       /*!< Returns a reference to the element index by <em>(row,col)</em> using a 1 based nomenclature*/
		SKRTFLOAT             p11                 ( ) const { return m_p[0]; }                       /*!< Returns the value p11, i.e. the scalar scattering coefficient (phase function) */
		void                  SetTo               ( double val );                                   /*!< Set all of the elements of the phase matix to \e val */
		void                  SetTo               ( float val );                                    /*!< Set all of the elements of the phase matix to \e val */
		void                  SetToIdentity       ( );                                          /*!< This becomes an identity matrix, with (4,4) still assumed to be zero. */
		SKTRAN_PhaseMat_MIMSNC&  RMultBy             ( const SKTRAN_ScatMat_Rot&  other );       /*!< This becomes the result of right multiplying itself by #other */
		SKTRAN_PhaseMat_MIMSNC&  LMultBy             ( const SKTRAN_ScatMat_Rot&  other );       /*!< This becomes the result of left multiplying itself by #other */
		SKTRAN_PhaseMat_MIMSNC&  RMultBy             ( const SKTRAN_ScatMat_MIMSNC& other );       /*!< This becomes the result of right multiplying itself by #other */
		SKTRAN_PhaseMat_MIMSNC&  LMultBy             ( const SKTRAN_ScatMat_MIMSNC& other );       /*!< This becomes the result of left multiplying itself by #other */
		SKTRAN_PhaseMat_MIMSNC&  RMultBy             ( const SKTRAN_PhaseMat_MIMSNC& other );       /*!< This becomes the result of right multiplying itself by #other */
		SKTRAN_PhaseMat_MIMSNC&  LMultBy             ( const SKTRAN_PhaseMat_MIMSNC& other );       /*!< This becomes the result of left multiplying itself by #other */
		SKTRAN_PhaseMat_MIMSNC&  operator +=         ( const SKTRAN_PhaseMat_MIMSNC& other);        /*!< Add the other phase matrix to this one */
		SKTRAN_PhaseMat_MIMSNC&  operator *=         ( double value);                                  /*!< Multiply the phase matrix by \e val */
		SKTRAN_PhaseMat_MIMSNC&  operator *=         ( float  value);                                  /*!< Multiply the phase matrix by \e val */
		SKTRAN_PhaseMat_MIMSNC&  operator =          ( const SKTRAN_ScatMat_Rot&  other);        /*!< Assign (copy) the other phase matrix to this one */
		SKTRAN_PhaseMat_MIMSNC&  operator =          ( const SKTRAN_ScatMat_MIMSNC& other);        /*!< Assign (copy) the other phase matrix to this one */
		SKTRAN_PhaseMat_MIMSNC&  operator =          ( const SKTRAN_PhaseMat_MIMSNC& other);        /*!< Assign (copy) the other phase matrix to this one */
        SKTRAN_Stokes_NC         operator *          ( const SKTRAN_Stokes_NC& stokes )     const;     /*!< Multiply this phase matrix by a Stokes Vector to get a Stokes vector */
		SKTRAN_PhaseMat_MIMSNC   operator *          ( double value)                        const;     /*!< Multiply this phase matrix by a constant and return the answer */
		SKTRAN_PhaseMat_MIMSNC   operator *          ( float  value)                        const;
		SKTRAN_PhaseMat_MIMSNC   operator -          ( const SKTRAN_PhaseMat_MIMSNC& other)   const;
		SKTRAN_PhaseMat_MIMSNC   operator +          ( const SKTRAN_PhaseMat_MIMSNC& other)   const;

        void LApplyTo( SKTRAN_Stokes_NC* v ) const;

		friend class SKTRAN_ScatMat_MIMSNC;
};