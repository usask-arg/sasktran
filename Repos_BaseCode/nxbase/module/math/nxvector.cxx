#include "nxbase_math.h"		// In house nxVector manipulation sofwtare.
/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/
using namespace::nxmath;


//----------------------------------------------------------------------------
//			nxVector::Constructor
//----------------------------------------------------------------------------

nxVector::nxVector()
{
  m_v[0] = 0.0;
  m_v[1] = 0.0;
  m_v[2] = 0.0;
}

//----------------------------------------------------------------------------
//			nxVector::UnitVector
//
//	Returns a unit vector in the same direction as this vector.
//	If the vector is zero, then it returns (0,0,0).
//
//----------------------------------------------------------------------------

nxVector nxVector::UnitVector() const
{
  double mag;
  nxVector unit;

  if (IsZero() || !IsValid() )
  {
     unit.SetCoords(0.0, 0.0, 0.0);
  }
  else
  {
     mag = Magnitude();
     unit = *this/mag;
  }
  return unit;
}
//----------------------------------------------------------------------------
//			nxVector::ComponentPerpendicularTo
//
//	Returns the component of this vector which is perpendicular to a
//	passed vector.
//----------------------------------------------------------------------------

nxVector nxVector::ComponentPerpendicularTo( const nxVector &Z ) const
{
   nxVector zaxis;
   nxVector xaxis;
   double zcomp;

   zaxis  = Z.UnitVector();					// Nominal direction of Freja spin axis in ECI coords.
   zcomp  = Dot( zaxis );					// Dot product of this vector along Z axis.
   xaxis  = *this - (zaxis*zcomp);				// Subtract this component, leaves component perpendicular.
   return xaxis;						// return component perpendicular.
}

//----------------------------------------------------------------------------
//			nxVector::ComponentParallelTo
//
//	Returns the component of this vector which is parallel to a
//	passed vector.
//----------------------------------------------------------------------------

nxVector nxVector::ComponentParallelTo( const nxVector &Z ) const
{
   nxVector zaxis;
   nxVector xaxis;
   double zcomp;

   zaxis  = Z.UnitVector();					// Nominal direction of Freja spin axis in ECI coords.
   zcomp  = Dot( zaxis );					// Dot product of this vector along Z axis.
   xaxis  = (zaxis*zcomp);					// Subtract this component, leaves component perpendicular.
   return xaxis;						// return component perpendicular.
}

//----------------------------------------------------------------------------
//			nxVector::TransformToNewPole
//
//	Converts the current X,Y,Z coordinate system into a new system
//	X',Y',Z'.  The Z' and X' coordinates are given in units of X,Y,Z
//
//	Assume the following,
//
//	Z = "Spin Axis pole" of system.
//	X = Location of zero degrees longitude.
//	Y = Direction pointing EAST at zero degrees longitude.
//
//----------------------------------------------------------------------------

void nxVector::TransformToNewPole( const nxVector &Xprime, const nxVector &Zprime )
{
   nxVector  Xp,Yp,Zp;
   double  x,y,z;

   Xp = Xprime.UnitVector();
   Zp = Zprime.UnitVector();
   Yp = Zp^Xp;

   x  = Dot( Xp );
   y  = Dot( Yp );
   z  = Dot( Zp );

   SetCoords( x, y, z );
}


//----------------------------------------------------------------------------
//			nxVector::TransformToNewPole
//
//	Converts the current X,Y,Z coordinate system into a new system
//	X',Y',Z'.  Notable examples are geocentric to geomagnetic coordinate
//	conversions and equatorial to topocentric conversions.
//
//	Assume the following,
//
//	Z = "Spin Axis" of system.
//	X = Location of zero degrees longitude.
//	Y = Direction pointing EAST at zero degrees longitude.
//
//	Then,
//	Theta = Longitude of X' in (X,Y,Z) coordinate system.
//	Chi   = Latitude  of Z' in (X,Y,Z) coordinate system.
//
//----------------------------------------------------------------------------

void nxVector::TransformToNewPole( double theta, double chi ) 
{
   double 	xprime;
   double	yprime;
   double	zprime;
   double	costheta;
   double	sintheta;
   double       coschi;
   double	sinchi;

   costheta = cosd(theta);
   sintheta = sind(theta);
   coschi   = cosd(chi);
   sinchi   = sind(chi);

   nxVector vrot_x(  sinchi*costheta,  sinchi*sintheta, -coschi );
   nxVector vrot_y( -sintheta,         costheta,          0.0    );
   nxVector vrot_z(  coschi*costheta,  coschi*sintheta, sinchi );

   xprime = Dot( vrot_x );
   yprime = Dot( vrot_y );
   zprime = Dot( vrot_z );

   SetCoords( xprime, yprime, zprime );
}
//----------------------------------------------------------------------------
//			nxVector::AngleTo
//
//	Returns the angle in degrees ( 0 to +180) between this vector and the
//	parameter vector.
//
//----------------------------------------------------------------------------

double nxVector::AngleTo( const nxVector &V2 ) const
{
   nxVector v1;
   nxVector v2;
   double cs;

   v1 = UnitVector();			// Unit vector of this vector.
   v2 = V2.UnitVector();		// Unit vector of second vector.
   cs = v1 & v2;			// Cos (angle) = dot product of two unit vectors.
   if (cs >  1.0) cs =  1.0;
   if (cs < -1.0) cs = -1.0;
   return acosd(cs);
}

//----------------------------------------------------------------------------
//			nxVector::RotateAboutXaxis
//
//	Rotates the Y and Z components about the X axis through an
//	angle theta degrees.  Angle is measured positive from Y axis
//	toward Z axis.
//
//	Converts the internal X,Y,Z components
//
//	Originally written for conversion between Equatorial and
//	ecliptic coordinates which are not efficiently calculated with
//	TransformToNewPole
//
//----------------------------------------------------------------------------

void nxVector::RotateAboutXaxis( double theta )
{
   double    costheta;
   double    sintheta;
   double    yprime;

   costheta = cosd(theta);
   sintheta = sind(theta);

   yprime  =  costheta*m_v[1]  + sintheta*m_v[2];
   m_v[2]      = -sintheta*m_v[1]  + costheta*m_v[2];
   m_v[1]      = yprime;
}

//----------------------------------------------------------------------------
//			nxVector::RotateAboutZaxis
//
//	Rotates the X and Y components about the Z axis through an
//	angle theta degrees.  Angle is measured positive from Y axis
//	toward Z axis.
//
//	Converts the internal X,Y,Z components
//
//	Originally written for conversion between Equatorial and
//	geographic coordinates which are not efficiently calculated with
//	TransformToNewPole
//
//----------------------------------------------------------------------------

void nxVector::RotateAboutZaxis( double theta )
{
   double    costheta;
   double    sintheta;
   double    xprime;

   costheta = cosd(theta);
   sintheta = sind(theta);

   xprime  =  costheta*m_v[0]  + sintheta*m_v[1];
   m_v[1]      = -sintheta*m_v[0]  + costheta*m_v[1];
   m_v[0]      = xprime;
}

//----------------------------------------------------------------------------
//			nxVector::Longitude
//
//	Returns the longitude of a given point in the current X,Y,Z
//	coordinate system.  Returns answer in range 0-360 degrees.
//----------------------------------------------------------------------------

double nxVector::Longitude() const
{
   double x;

   if (( m_v[0] == 0.0) && (m_v[1] == 0.0))
   {
      return 0.0;
   }
   else
   {
      x = atan2d( m_v[1],m_v[0]);
      x = inrange(x,360.0);
      return x;
   }
}

//----------------------------------------------------------------------------
//			nxVector::Latitude
//
//	Returns the latitude of a given point in the current X,Y,Z coordinate
//	system.
//
//----------------------------------------------------------------------------

double nxVector::Latitude() const
{
  double ss;
  double mag;

  mag = Magnitude();
  if (mag > 0.0) ss = asind( m_v[2]/mag);
  else           ss = 0.0;
  return ss;
}

//----------------------------------------------------------------------------
//			nxVector::EquatorialToGeographic
//
//	Converts the current vector from an Equatorial system to a Geographic
//	System, using the mean vernal equinox.
//
//----------------------------------------------------------------------------

nxVector nxVector::EquatorialToGeographic( const nxTimeStamp &Tnow ) const
{
   double  theta;
   nxVector  tmp;

   theta = Tnow.GMST()*360.0;
   tmp   = *this;
   tmp.RotateAboutZaxis(theta);
   return tmp;
}

//----------------------------------------------------------------------------
//			nxVector::GeographicToEquatorial
//
//	Converts the current vector from an Geographic system to a Equatorial
//	System, using the mean vernal equinox.
//
//----------------------------------------------------------------------------

nxVector nxVector::GeographicToEquatorial( const nxTimeStamp &Tnow ) const
{
   double  theta;
   nxVector  tmp;

   theta = Tnow.GMST()*360.0;
   tmp   = *this;
   tmp.RotateAboutZaxis( -theta);
   return tmp;
}


//----------------------------------------------------------------------------
//			nxVector::GeographicToGeomagneticDipole
//
//	Returns the vector in Centred Geomagnetic Dipole
//	coordinates.  I extracted the coordinates of the centred dipole
// 	from the MSIS86 fortran model.
//
//	The original (geographic) coordinates of the vector are unchanged.	
//----------------------------------------------------------------------------

nxVector nxVector::GeographicToGeomagneticDipole() const
{
   static const double GGMLatOffset  =  90.0 - 11.4;
   static const double GGMLongOffset = -69.8;
   nxVector ggm;


   ggm = *this;
   ggm.TransformToNewPole( GGMLongOffset, GGMLatOffset);
   return ggm;
}


//----------------------------------------------------------------------------
//			nxVector::Dot
//----------------------------------------------------------------------------

double nxVector::Dot(   const nxVector &v2) const
{
   return   m_v[0]*v2.m_v[0]
	  + m_v[1]*v2.m_v[1]
	  + m_v[2]*v2.m_v[2];
}

//------------------------------------------------------------------------------
//			nxVector::Magnitude
//------------------------------------------------------------------------------

double nxVector::Magnitude() const
{
  return sqrt( m_v[0]*m_v[0] + m_v[1]*m_v[1] + m_v[2]*m_v[2] );
}

//
// --------------- nxVector::Cross --------------------
//


nxVector nxVector::Cross(  const nxVector &v2) const
{
  nxVector v;

  v.m_v[0] = m_v[1]*v2.m_v[2] - m_v[2]*v2.m_v[1];
  v.m_v[1] = m_v[2]*v2.m_v[0] - m_v[0]*v2.m_v[2];
  v.m_v[2] = m_v[0]*v2.m_v[1] - m_v[1]*v2.m_v[0];
  return v;
}

//-------------------------------------------------------------------------------
//			nxVector::operator+   add two vectors
//-------------------------------------------------------------------------------

nxVector nxVector::operator+(  const nxVector &v2) const
{
 nxVector v;

 v.m_v[0] = m_v[0] + v2.m_v[0];
 v.m_v[1] = m_v[1] + v2.m_v[1];
 v.m_v[2] = m_v[2] + v2.m_v[2];
 return v;
}

//-------------------------------------------------------------------------------
//			vector::operator-
//	Unary negation of the vector
//-------------------------------------------------------------------------------

nxVector	nxVector::operator-() const
{
	nxVector v( -m_v[0], -m_v[1], -m_v[2] );
	return v;
}

//-------------------------------------------------------------------------------
//			nxVector::operator-  subtract nxVector from nxVector.
//-------------------------------------------------------------------------------

nxVector nxVector::operator-(  const nxVector &v2) const
{
 nxVector v;

 v.m_v[0] = m_v[0] - v2.m_v[0];
 v.m_v[1] = m_v[1] - v2.m_v[1];
 v.m_v[2] = m_v[2] - v2.m_v[2];
 return v;
}

//-------------------------------------------------------------------------------
//			nxVector::operator+  add scaler to nxVector.
//-------------------------------------------------------------------------------

nxVector nxVector::operator+(  const double val ) const
{
 nxVector v;

 v.m_v[0] = m_v[0] + val;
 v.m_v[1] = m_v[1] + val;
 v.m_v[2] = m_v[2] + val;
 return v;
}

//-------------------------------------------------------------------------------
//			nxVector::operator- subtract scaler from nxVector
//-------------------------------------------------------------------------------

nxVector nxVector::operator-(  const double val ) const
{
 nxVector v;

 v.m_v[0] = m_v[0] - val;
 v.m_v[1] = m_v[1] - val;
 v.m_v[2] = m_v[2] - val;
 return v;
}

//-------------------------------------------------------------------------------
//			nxVector::operator/ divide nxVector by scalar
//-------------------------------------------------------------------------------

nxVector nxVector::operator/(  const double val ) const
{
	nxVector v( m_v[0]/val, m_v[1]/val, m_v[2]/val );
	return v;
}


 nxVector 	operator*( const nxVector& v1, const double val ) 	// multiply nxVector by a scalar.
{
	nxVector v( v1.m_v[0]*val, v1.m_v[1]*val, v1.m_v[2]*val );
	return v;
}

	
 nxVector 	operator*( const double val, const nxVector& v1  ) 	// multiply nxVector by a scalar.
{
	nxVector v( v1.m_v[0]*val, v1.m_v[1]*val, v1.m_v[2]*val );
	return v;
}


//----------------------------------------------------------------------------
//			nxVector::operator^
//
//	cross product of two vectors
//----------------------------------------------------------------------------

nxVector nxVector::operator^(  const nxVector &v2) const
{
   return Cross( v2);
}

//----------------------------------------------------------------------------
//			nxVector::operator&	Dot product of two vectors.
//
//	Dot product of two vectors.
//
//----------------------------------------------------------------------------

double nxVector::operator&(  const nxVector &v2) const
{
   return Dot( v2 );
}

//----------------------------------------------------------------------------
//			nxVector::operator!
//
//	Magnitude of a nxVector.
//
//----------------------------------------------------------------------------

double nxVector::operator!() const
{
   return Magnitude();
}

//----------------------------------------------------------------------------
//			nxVector::IndexOfMaxComponent
//
//   returns a number 0,1 or 2 corresponding to the index of the maximum
//   component.
//
//   0 means X is largest component
//   1 means Y is largest component
//   2 means Z is largest component
//----------------------------------------------------------------------------

   
 int nxVector::IndexOfMaxComponent()
 {
   double maxval;
   int    maxidx;
    
   if ( fabs(m_v[0]) > fabs(m_v[1]) )
   {
      maxval = fabs(m_v[0]);
      maxidx = 0;
   }
   else
   {
      maxval = fabs(m_v[1]);
      maxidx = 1;
   }
   if (maxval < fabs(m_v[2]) )
   {
 //     maxval = fabs(m_v[2]);
      maxidx = 2;
   }
   return maxidx;
}

void nxVector::FromLatLong( double latitude, double longitude, double magnitude )
{
	double h;

	m_v[2] = sind(latitude)*magnitude;
	h = cosd(latitude)*magnitude;
	m_v[0] = cosd(longitude)*h;
	m_v[1] = sind(longitude)*h;
}

   
