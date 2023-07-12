#if !defined (NXBASE_NXVECTOR_H)
#define NXBASE_NXVECTOR_H 1

#include "../science/geodesy/timestmp.h"

/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/

//#include "ndlmath.h"
//#include "timestmp.h"

typedef double NXVECTOR[3];

/*----------------------------------------------------------------------------
 *					class nxVector											*/
/** \ingroup nxMath_Vector
 *	Implements 3-D vectors as X,Y,Z
**/
/*----------------------------------------------------------------------------*/

class nxVector
{
   private:
	   NXVECTOR		m_v;
	public:
      				nxVector();																				//!< Default constructor  Set the vector to zero				
       				nxVector ( double x, double y, double z) {SetCoords(x,y,z);}							//!< Custom Constructor Set the coordinates to xyz
       void			FromSequence( const double fixedarray[3] )	 {SetCoords(fixedarray[0], fixedarray[1], fixedarray[2]);}	//!< Assignment method specifically for Python sequences via SWIG
	   nxVector		AsSequence  () const						 {return *this;}											//!< Returns the nxVector as a sequence to Python via SWIG
       void   		SetCoords( double x, double y, double z) { m_v[0]=x; m_v[1]=y; m_v[2]=z;}				//!< Set the coordinates of the vector
       void   		SetCoords( const NXVECTOR* other) { m_v[0] = (*other)[0]; m_v[1] = (*other)[1]; m_v[2] = (*other)[2];}	//!< Set the coordinates of the vector from a pointer to a 3 element array
	   void			FromLatLong ( double latitude, double longitude, double magnitude = 1.0 );				//!< Set x,y,z vector to geocentric latitude (degrees), longitude (degrees) and magnitude
																				
//     nxBOOL  		isNULL	( ) const { return ( (m_v[0]==0.0) && (m_v[1]==0.0) && (m_v[2]==0.0));}			//!< Return true if the vector is zero
       bool  		IsZero	( ) const { return ( (m_v[0]==0.0) && (m_v[1]==0.0) && (m_v[2]==0.0));}			//!< Return true if the vector is zero
       bool  		IsValid ( ) const { return ( NXFINITE(m_v[0]) && NXFINITE(m_v[1]) && NXFINITE(m_v[2]));}			//!< Return true if the vector is zero
	   void			SetInvalid () { m_v[0]= std::numeric_limits<double>::quiet_NaN(); m_v[1] = m_v[0]; m_v[2] = m_v[0];}
	   double		X(void) const {return m_v[0];}															//!< return the X component
	   double		Y(void) const {return m_v[1];}															//!< return the Y component
	   double		Z(void) const {return m_v[2];}															//!< return the Z component
	   NXVECTOR&	Data()        {return m_v;}

       double 		AngleTo							( const nxVector & V2 )		const;						//!< Get the angle of this vector to another vector (0-180 degrees)
       nxVector 	UnitVector						( void )					const;						//!< Return the unitvector of this vector
       nxVector 	ComponentPerpendicularTo		( const nxVector &Z )		const;						//!< Get the component of this vector (cross product) perpendicular to the passed in vector "Z"
       nxVector 	ComponentParallelTo				( const nxVector &Z )		const;						//!< Get the component of this vector (dot product) parallel to the passed in vector "Z"
       int    		IndexOfMaxComponent				( void );												//!< Return the index o fthe maximum component, x=0, y=1, z=2
       double 		Longitude						( void )					const;						//!< Get the geocentric longitude (degrees) of this vector
       double 		Latitude						( void )					const;						//!< Get the geocentric latitude (degrees) of this vector
       nxVector 	EquatorialToGeographic			( const nxTimeStamp &Tnow )	const;						//!< Convert from equatorial (ECI) coordinates to Geographic coordinates
       nxVector 	GeographicToEquatorial			( const nxTimeStamp &Tnow )	const;						//!< Convert from Geographic coordinate to quatorial (ECI) coordinates
       nxVector 	GeographicToGeomagneticDipole	( void ) const;											//!< Convert from Geographic coordinates to Geomagnetic Dipole coordinates
       void   		TransformToNewPole				( const nxVector &XPrime, const nxVector &ZPrime );		//!< Transform this vector to a new rotated coordinate system, XPrimt and YPrime
       void   		TransformToNewPole				( double theta, double chi );							//!< Transform this vector to a new rotated coordinate system , theta, chi
       void   		RotateAboutXaxis				( double theta );										//!< Rotate this vector an angle theta degrees about the X axis
       void   		RotateAboutZaxis				( double theta );										//!< Rotate this vector an angle theta degrees about the Z axis

       double 		Dot		( const nxVector &v2) const;													//!< Get the dot product of this vector with v2
       double 		Magnitude(void) const;																	//!< Get the magnitude of this vector
       nxVector 	Cross( const nxVector &v2 ) const;     													//!< Get the cross product of this vector with v2 (this^v2).

	   nxVector&		operator=  ( const NXVECTOR* other) { SetCoords(other); return *this;}						//!< Set this element to the 3 element array whose address is given
       nxVector 		operator+  ( const nxVector &v2 ) const;													//!< return (*this + v2) 
	   nxVector&		operator+= ( const nxVector &v2 ) { m_v[0] += v2.m_v[0]; m_v[1] += v2.m_v[1]; m_v[2] += v2.m_v[2]; return *this;}	//!< (*this += v2)
       nxVector&		operator-= ( const nxVector &v2 ) { m_v[0] -= v2.m_v[0]; m_v[1] -= v2.m_v[1]; m_v[2] -= v2.m_v[2]; return *this;}	//!< (*this -= v2)
	   nxVector&		operator*= ( const nxVector &v2 ) { m_v[0] *= v2.m_v[0]; m_v[1] *= v2.m_v[1]; m_v[2] *= v2.m_v[2]; return *this;}	//!< (*this*v2)
	   nxVector&		operator*= ( const double val )   { m_v[0] *= val;     m_v[1] *= val; 	m_v[2] *= val; return *this;}	//!< (*this*val)
	   nxVector&		operator/= ( const double val )   { NXASSERT(0.0 != val); m_v[0] /= val; m_v[1] /= val; m_v[2] /= val; return *this; }
       nxVector 		operator-( const nxVector &v2 ) const;														//  return (*this - v2)
       nxVector			operator-(void)const;
       nxVector 		operator+( const double val ) const;				//!< add a scalar from a nxVector.
       nxVector 		operator-( const double val ) const;			 	//!< subtract a scalar from a nxVector
       nxVector 		operator/( const double val ) const;			 	//!< divide nxVector by a scalar.
       nxVector 		operator^( const nxVector &v2 ) const;				//!< cross product of two vectors.
       double 			operator&( const nxVector &v2 ) const;				//!< Dot product of two vectors.
       double 			operator!(void)const ;								//!< Magnitude of a nxVector.
						operator const NXVECTOR* () const 	{return &m_v;}	//!< Typecast to pointer to an array of 3 elements
friend  nxVector operator*( const nxVector& v1, const double val ); 		//!< multiply nxVector by a scalar.
friend  nxVector operator*( const double val, const nxVector& v1  ); 		//!< multiply nxVector by a scalar.
};

 nxVector 	operator*( const nxVector& v1, const double val ); 				//!< multiply nxVector by a scalar.
 nxVector 	operator*( const double val, const nxVector& v1  );				//!< multiply nxVector by a scalar.

#endif




