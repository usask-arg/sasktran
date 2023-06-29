
typedef double			QUATERNION[4];						/* Define the quaternions used for attitude */
#define NUMQUATERNION	(sizeof(QUATERNION)/sizeof(double))


/*---------------------------------------------------------------------------
 *						class nxQuaternion									*/
/**	\ingroup nxMath_Quaternion
 *	A class for manipulating quaternions.
*/
/*---------------------------------------------------------------------------*/

class nxQuaternion
{
	private:
		double			m_q[4];

	public:
						nxQuaternion			();
						nxQuaternion			( double a, double b, double c, double d);
						nxQuaternion			( double a, const nxVector& qv );
	public:
		nxQuaternion	operator*				( const nxQuaternion& other );
		nxQuaternion&	operator*=				( const nxQuaternion& other );
		double&			operator[]				( int index) { return m_q[index];}
		void			FromScalarAndVector		( double a, const nxVector& qv);
		void			FromScalars				( double a, double b, double c, double d);
		void			FromQUATERNION			( const QUATERNION& Q );
		void			FromEuler				( const nxVector& eulerdegrees);
		void			FromAxisToAxisRotation  ( const nxVector& from, const nxVector& to );
		void			SetToIdentity			();
		double			Norm					() const;
		void			Normalize				();
		double			ScalarComponent			() const ;
		nxVector		VectorComponent			() const ;
		double			DotProduct				( const nxQuaternion& other) const;
		nxQuaternion	Conjugate				();
		nxVector		RotateVector			( const nxVector& original );
		nxQuaternion	UnitSlerp				( const nxQuaternion& endquat, double alpha) const;
		double			X						() const		{ return m_q[1];}
		double			Y						() const		{ return m_q[2];}
		double			Z						() const		{ return m_q[3];}
		double			W						() const		{ return m_q[0];}



};

