#include "nxbase_geodesy.h"

size_t SKTRAN_CoordinateTransform_V2::m_numinstances = 0;
size_t SKTRAN_GridDefBase_V2::m_numinstances         = 0;


/*-----------------------------------------------------------------------------
 *					HELIODETIC_UNITVECTOR::FromVector		2013-10-28*/
/** **/
/*---------------------------------------------------------------------------*/

HELIODETIC_UNITVECTOR& HELIODETIC_UNITVECTOR::FromVector	(const HELIODETIC_VECTOR &v)
{
	double	f = 1.0/v.Magnitude();;

	SetCoords( v.X()*f, v.Y()*f, v.Z()*f );
	return *this;
}

/*-----------------------------------------------------------------------------
 *					HELIODETIC_VECTOR::SetCoords		2009-1-26*/
/** **/
/*---------------------------------------------------------------------------*/

void HELIODETIC_VECTOR::SetCoords	( double x, double y, double z )
{
	m_data[0] = x;
	m_data[1] = y;
	m_data[2] = z;
}


/*-----------------------------------------------------------------------------
 *					HELIODETIC_VECTOR::SetCoords		2009-1-26*/
/** **/
/*---------------------------------------------------------------------------*/

void HELIODETIC_VECTOR::SetCoords( const HELIODETIC_UNITVECTOR& unit, double magnitude )
{
	m_data[0] = unit.X()*magnitude;
	m_data[1] = unit.Y()*magnitude;
	m_data[2] = unit.Z()*magnitude;
}



/*-----------------------------------------------------------------------------
 *					HELIODETIC_VECTOR::SetCoords		2009-1-26*/
/** **/
/*---------------------------------------------------------------------------*/

void HELIODETIC_VECTOR::SetCoords	( const HELIODETIC_UNITVECTOR& unit)
{
	m_data[0] = unit.X();
	m_data[1] = unit.Y();
	m_data[2] = unit.Z();
}


/*-----------------------------------------------------------------------------
 *					HELIODETIC_VECTOR::UnitVector		2009-1-26*/
/** **/
/*---------------------------------------------------------------------------*/

HELIODETIC_UNITVECTOR HELIODETIC_VECTOR::UnitVector() const
{
	HELIODETIC_UNITVECTOR	unit;
	double					f;

	f = Magnitude();
	if (f > 0.0) f = 1.0/f;
	unit.SetCoords( m_data[0]*f, m_data[1]*f, m_data[2]*f );
	return unit;
}



/*-----------------------------------------------------------------------------
 *					HELIODETIC_POINT::Initialize		2009-1-26*/
/** **/
/*---------------------------------------------------------------------------*/

void HELIODETIC_POINT::Initialize( const HELIODETIC_UNITVECTOR& unit, double magnitude, const SKTRAN_CoordinateTransform_V2* coords )
{
	m_direction = unit;
	m_radius    = magnitude;
	m_heightm   = (coords != nullptr) ? coords->RadiusToAltitude(m_radius) : std::numeric_limits<double>::quiet_NaN();

#if defined(NXDEBUG)
	HELIODETIC_VECTOR	v;
	v.SetCoords( m_direction, m_radius );
	if (coords != nullptr)
	{
		m_geo = coords->HelioVectorToGeographic( v);
	}
#endif
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayBaseGeometry_V21::LocationAsPoint		2013-10-28*/
/** **/
/*---------------------------------------------------------------------------*/

void HELIODETIC_POINT::FromVector( const HELIODETIC_VECTOR&	v, const SKTRAN_CoordinateTransform_V2* coords )
{
	HELIODETIC_UNITVECTOR	unit;
	double					r;
	double					f;

	r = v.Magnitude();
	f = (1.0/r);
	unit.SetCoords( v.X()*f, v.Y()*f, v.Z()*f );
	Initialize( unit, r, coords );
}


/*-----------------------------------------------------------------------------
 *					HELIODETIC_POINT::Clear		2009-1-26*/
/** **/
/*---------------------------------------------------------------------------*/

void HELIODETIC_POINT::Clear()
{
	m_direction.Clear();
	m_radius =-9999999.0;
	m_heightm = -9999999.0;
}


/*-----------------------------------------------------------------------------
 *					HELIODETIC_POINT::CosZenithAngle		2009-1-26*/
/** **/
/*---------------------------------------------------------------------------*/

double HELIODETIC_POINT::CosZenithAngle( const HELIODETIC_UNITVECTOR& look ) const
{
	double	coszen;

	coszen = look & m_direction;
	if (coszen > 1.0)
	{
		NXASSERT(( coszen < 1.000001 ));
		coszen = 1.0;
	}
	if (coszen < -1.0)
	{
		NXASSERT(( coszen > -1.000001 ));
		coszen = -1.0;
	}
	return coszen;
}


/*-----------------------------------------------------------------------------
 *					HELIODETIC_POINT::LocalUnitVectors		2009-1-26*/
/** Gets the local unit vectors**/
/*---------------------------------------------------------------------------*/

bool HELIODETIC_POINT::LocalUnitVectors( HELIODETIC_UNITVECTOR* localunitarray, size_t numv)  const
{
	double					coszen = m_direction.Z();
	double					sinzen = sqrt( 1.0 - coszen*coszen );
	double					sinzenrecip;
	HELIODETIC_UNITVECTOR*	north = localunitarray;
	HELIODETIC_UNITVECTOR*	west  = localunitarray + 1;
	HELIODETIC_UNITVECTOR*	up    = localunitarray + 2;
	bool					ok;

	ok = (numv == 3);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"HELIODETIC_POINT::LocalUnitVectors, You must pass in exactly 3 unit vetor array. You passed in %d", (int)numv);
	}
	else
	{
		if (sinzen < 0.005 )				// Sun and location are with 0.28 degrees of each other
		{
			north->SetCoords( -1,  0, 0 );		// North
			west-> SetCoords(  0, -1, 0 );		// West
			up->   SetCoords(  0,  0, 1 );		// Up
			if( coszen < 0.0 ){
				north->Negate();
				west-> Negate();
				up->   Negate();
			}
		}
		else
		{
			sinzenrecip = 1.0/sinzen;
			west->SetCoords(     m_direction.Y()*sinzenrecip,		// west = direction ^ sun
								-m_direction.X()*sinzenrecip,		// But sun is along Z axis = (0,0,1)
								 0 );								// Therefore cross product is really quite simple.

			north->SetCoords(	west->Y()*m_direction.Z(),			//north = west ^ up
							   -west->X()*m_direction.Z(),
								west->X()*m_direction.Y() - west->Y()*m_direction.X() );

			*up = m_direction;																	// Up is same as this unit vector
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					HELIODETIC_POINT::TransformToLocalZenithCoords		2009-1-26*/
/** **/
/*---------------------------------------------------------------------------*/

HELIODETIC_UNITVECTOR HELIODETIC_POINT::TransformToLocalZenithCoords( const HELIODETIC_UNITVECTOR& v, const HELIODETIC_UNITVECTOR* localunitarray ) const
{
	const HELIODETIC_UNITVECTOR&	north = localunitarray[0];
	const HELIODETIC_UNITVECTOR&	west  = localunitarray[1];
	const HELIODETIC_UNITVECTOR&	up    = localunitarray[2];
	HELIODETIC_UNITVECTOR			newv;

	double	x,y,z;

	x = v & north;
	y = v & west;
	z = v & up;
	newv.SetCoords(x,y,z);
	return newv;
}


/*-----------------------------------------------------------------------------
 *					HELIODETIC_POINT::Vector		2009-1-26*/
/** **/
/*---------------------------------------------------------------------------*/

HELIODETIC_VECTOR HELIODETIC_POINT::Vector() const
{
	HELIODETIC_VECTOR	v;

	v.SetCoords( m_direction, m_radius );
	return v;
}

/*-----------------------------------------------------------------------------
 *					HELIODETIC_POINT::Vector		2014-06-02*/
/** **/
/*---------------------------------------------------------------------------*/

HELIODETIC_UNITVECTOR HELIODETIC_POINT::UnitVector() const
{
	return m_direction;
}


/*-----------------------------------------------------------------------------
 *					HELIODETIC_BASIS::ProduceBasis		 2014- 12- 10*/
/** **/
/*---------------------------------------------------------------------------*/

bool HELIODETIC_BASIS::ProduceBasis ( const HELIODETIC_POINT& pt, const HELIODETIC_UNITVECTOR& look ) 
{
	HELIODETIC_UNITVECTOR	loc = pt.UnitVector();			// Unit vector to the point

	return ProduceBasis( loc, look );
}
bool HELIODETIC_BASIS::ProduceBasis ( const HELIODETIC_VECTOR& pt, const HELIODETIC_UNITVECTOR& look ) 
{
	HELIODETIC_UNITVECTOR	loc = pt.UnitVector();			// Unit vector to the point

	return ProduceBasis( loc, look );
}



bool HELIODETIC_BASIS::ProduceBasis ( const HELIODETIC_UNITVECTOR&  , const HELIODETIC_UNITVECTOR& look )
{
	x.SetCoords( -look.X(), -look.Y(), -look.Z() );
	HELIODETIC_VECTOR temp;
	temp.SetCoords( -look.Z()*look.X(), -look.Y()*look.Z(), 1.0 - look.Z()*look.Z() ); // = sunDir - look*dot(sunDir,propDir)
	if(1e-6<temp.Magnitude()){
		y = temp.UnitVector();
		temp.SetCoords(x.Y()*y.Z() - x.Z()*y.Y(), x.Z()*y.X() - x.X()*y.Z(), x.X()*y.Y() - x.Y()*y.X() );
		z = temp.UnitVector();
	} else{
		HELIODETIC_UNITVECTOR unit_x;
		unit_x.SetCoords( 1.0, 0.0, 0.0 ); // Sun direction is z, terminator direction is y, so choose x as perpendicular to both
		double norm;
		norm = x.X()*unit_x.X() + x.Y()*unit_x.Y() + x.Z()*unit_x.Z();
		double ax, ay, az;
		if( 1e-6 < 1 - norm*norm ){
			// Can take component of x perpendicular to look
			ax = unit_x.X() - norm*x.X();
			ay = unit_x.Y() - norm*x.Y();
			az = unit_x.Z() - norm*x.Z();
			norm = sqrt( ax*ax + ay*ay + az*az );
			ax = ax / norm;
			ay = ay / norm;
			az = az / norm;
			y.SetCoords( ax, ay, az );
		} else{
			// Use component of y perpendicular to look instead
			HELIODETIC_UNITVECTOR unit_y;
			unit_y.SetCoords( 0.0, 1.0, 0.0 );
			norm = x.X()*unit_y.X() + x.Y()*unit_y.Y() + x.Z()*unit_y.Z();
			ax = unit_y.X() - norm*x.X();
			ay = unit_y.Y() - norm*x.Y();
			az = unit_y.Z() - norm*x.Z();
			norm = sqrt( ax*ax + ay*ay + az*az );
			ax = ax / norm;
			ay = ay / norm;
			az = az / norm;
			y.SetCoords( ax, ay, az );
	//		nxLog::Record(NXLOG_WARNING, "HELIODETIC_BASIS::ProduceBasis, Took unit_y path ().\n");
		}

		ax = x.Y()*y.Z() - x.Z()*y.Y();
		ay = x.Z()*y.X() - x.X()*y.Z();
		az = x.X()*y.Y() - x.Y()*y.X();
		norm = sqrt( ax*ax + ay*ay + az*az );
		ax = ax / norm;
		ay = ay / norm;
		az = az / norm;

		z.SetCoords(ax,ay,az);
	}

	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefBase_V2::SKTRAN_GridDefBase_V2		2007-11-21*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_GridDefBase_V2::SKTRAN_GridDefBase_V2()
{
	++m_numinstances;
//	m_gridvalues    = NULL;
//	m_numgridpoints = 0;
	m_gridsearchmode = GRIDSEARCH_NONUNIFORM;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefBase_V2::~SKTRAN_GridDefBase_V2		2007-11-21*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_GridDefBase_V2::~SKTRAN_GridDefBase_V2()
{
	ReleaseResources();
	--m_numinstances;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefBase_V2::ReleaseResources		2007-11-21*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_GridDefBase_V2::ReleaseResources()
{
	m_gridvalues.clear();
//	if (m_gridvalues != NULL) delete [] m_gridvalues;
//	m_gridvalues = NULL;
//	m_numgridpoints = 0;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefBase_V2::AllocateGridArray		2007-11-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_GridDefBase_V2::AllocateGridArray( size_t numpoints )
{
	bool	ok;

	ReleaseResources();
	ok = (numpoints == 0);
	if (!ok)
	{
		m_gridvalues.resize(numpoints);
		ok = (m_gridvalues.size() == numpoints);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "SKTRAN_GridDefBase_V2::AllocateGridArray, Error allocating space for %Iu elements", (size_t)numpoints );
			ReleaseResources();
		}
	}

	SetGridSearchMode_NonUniform( );

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefBase_V2::CopyGridArray		2007-11-21*/
/** **/
/*---------------------------------------------------------------------------*/


bool SKTRAN_GridDefBase_V2::CopyGridArray(  const double* source, size_t numpoints )
{
	m_gridvalues.assign( source, source+numpoints);
	SetGridSearchMode_NonUniform( );
	return true;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefBase_V2::CopyGridArray		2014-1-31*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_GridDefBase_V2::CopyGridArray(  const std::vector<double>& source)
{
	m_gridvalues = source;
	SetGridSearchMode_NonUniform( );
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefBase_V2::DeepCopy		2009-6-17*/
/** Copies the Grid Array from other into this instance. Note there is no
 *	checking of lock counts or anything else. The User must make sure no one
 *	else is reliant upon the underlying pointers not being changed.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_GridDefBase_V2::DeepCopy( const SKTRAN_GridDefBase_V2& other )
{
	bool		ok = true;

	m_gridvalues = other.m_gridvalues; 
	//ok = CopyGridArray( other.m_gridvalues, other.m_numgridpoints );
	SetGridSearchMode( other.GetGridSearchMode( ) );

	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefBase_V2::IndexOfPointEqualOrAbove		2007-11-21*/
/** Retrieves the index of the grid point whose value is equal or just greater than
 *	the value passed in.  Returns false if no element is greater than value
 *	and sets the index to the last element of the array.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_GridDefBase_V2::IndexOfPointEqualOrAbove( double value, SKTRAN_GridIndex* index ) const
{
	const_iterator start  = m_gridvalues.begin();
	const_iterator	finish = m_gridvalues.end();
	const_iterator	iter;
	bool			ok;

	iter   = LowerBound( start, finish, value );
	ok     = !(iter == finish );
	if (ok) *index = iter - start;
	else    *index = m_gridvalues.size()-1;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefBase_V2::IndexOfPointBelow		2007-11-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_GridDefBase_V2::IndexOfPointBelow( double value, SKTRAN_GridIndex* index ) const
{
	const_iterator	start  = m_gridvalues.begin();
	const_iterator	finish = m_gridvalues.end();
	const_iterator	iter;
	bool			ok;

	iter   = LowerBound( start, finish, value );		// find point greater than or equal
	ok     = !(iter == start);								// Make sure we are not at the start
	if (ok) *index = (--iter - start);						// go down 1 to get point definitely below
	else    *index = 0;
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefBase_V2::IndexOfPointBelowOrEqual		2008-1-31*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_GridDefBase_V2::IndexOfPointBelowOrEqual( double value, SKTRAN_GridIndex* index ) const
{
	const_iterator	start  = m_gridvalues.begin();
	const_iterator	finish = m_gridvalues.end();
	const_iterator	iter;
	bool			ok;

	iter   = UpperBound( start, finish, value );				// Find the point definitely greater than this value
	ok     = !(iter == start);
	if (ok) *index = (--iter - start);
	else    *index = 0;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefBase_V2::FindBoundingIndices		2008-1-10*/
/** Finds the indices in the grid that bound a given value x and returns
 *	the bounding indices. It also returns the weights that should be applied
 *	to do a linear interpolation of the indexed values.
 *
 *	It always fails if the grid has less than 2 points.  It never fails for more than
 *	two points if allowoutofrange is true.  If allowoutofrange is false then the
 *	code will return fail if x is outside the bounds of the grid but the indices
 *	and interpolation indices are still set to the best values.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_GridDefBase_V2::FindBoundingIndices( double x, ENUM_INTERPOLATIONMODE outrange, SKTRAN_GridIndex* lowercell, double* lowerweight, SKTRAN_GridIndex* uppercell, double* upperweight) const
{
	const_iterator		x1;							// The value just above our value of x
	const_iterator		x0;							// The value just below our value of x
	const_iterator		start;
	const_iterator		finish;
	double				xv;
	double				dx;
	bool				ok;
	bool				outofbounds;

//	NXASSERT(( x < 250000.0 ));		// Look out for old Earth radii values.
	ok = (NumGridPoints() > 1);
	if (!ok)
	{
		ok = (NumGridPoints() == 1) && (outrange == OUTOFBOUND_TRUNCATE);
		if (ok)
		{
			*lowerweight = 1;
			*lowercell   = 0;
			*upperweight = 0;
			*uppercell   = 0;
		}
		else
		{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_GridDefBase_V2::FindBoundingIndices, the grid is too small to define upper and lower bounds");
			*uppercell   = 0;							// Get the index of the upper cell
			*lowercell   = 0;							// Get the index of the lower cell
			*upperweight = 0;
			*lowerweight = 0;
		}
	}
	else
	{
		xv = x;
		outofbounds = (xv <  m_gridvalues.front()) || (xv >= m_gridvalues.back());			// this is OK as long as it is not out of bounds
		if (outofbounds)																// if it is out of bounds
		{																				// then
			if (xv < m_gridvalues.front() ) xv =  m_gridvalues.front();
			if (xv >= m_gridvalues.back() )
			{
				xv = m_gridvalues.at(NumGridPoints()-2);
			}
		}																	// we now havethe upper value as sensible;
		start  = m_gridvalues.begin();										// get the start of the altitude grid
		finish = m_gridvalues.end();										// Get the end of the altitude grid
		x1     = UpperBound( start, finish, xv );							// **** upper_bound is important as there are often 2 ground (diffuse) points with the same altitude (we want the second one) . Dont use lower_bound.  Find the pointer to the value greater than x

		x0           = x1 - 1;									// get the lower value
		*uppercell   = (x1 - start);							// Get the index of the upper cell
		*lowercell   = (x0 - start);							// Get the index of the lower cell
		dx           = *x1 - *x0;								// Get the distance between the two points in the grid
		NXASSERT((dx > 0 ));									// It should always be positive if the grid is in ascending order.

		*upperweight   =  (x - *x0)/dx;							// Get the linear interpolation weight to be applied to the upper index
		*lowerweight   = 1.0 - *upperweight;					// Get the linear interpolation weight to be applied to the lower index

		ok = (!outofbounds) || (outrange == OUTOFBOUND_EXTRAPOLATE);
		if (!ok)													// If we are out of bound but not in interpolation mode (by default, ok is false so we have implemented OUTOFBOUND_ERROR)
		{															// then
			if (outrange == OUTOFBOUND_ZERO)					// if we dont want to go outside the array
			{													// then
				*lowerweight = 0;								// clear all the weights
				*upperweight = 0;								// so there is no contribution
				ok           = true;							// And this interpolation is ok
			}
			else if	(outrange == OUTOFBOUND_TRUNCATE)			// other option is OUTOFBOUND_TRUNCATE
			{													// so
				ok = true;										// this is ok
				if (x0 == start)								// if we are at the beginning
				{												// then
					*lowerweight = 1;							// all the weight is on the first point in the array
					*upperweight = 0;							// no weight on the 2nd point
				}												// case of OUTOFBOUND_ERROR is alread done with ok returning as FALSE
				else if (x1 == (finish-1))						// otherwise if we are beyond the last point
				{												// then
					*lowerweight = 0;							// set the first weight to zero
					*upperweight = 1;							// set the last weight to 1
				}												// case of OUTOFBOUND_ERROR is already handled as ok is already false
			}
			else if ( outrange == OUTOFBOUND_ERROR)				// if we are out of bounds then error only if requested
			{													// and if we are out of bounds
				ok = ( x <= back() ) || ( x >= front());		// we occassionaly get an exact floating point compare ( so check and correct)
			}
			else
			{
				NXASSERT((false));								// new OUTOFBOUND otpion has not been implemented in this code
			}
		}														// done handling all of the out of bound errors
	}															// done handling the linear interpolation
	return ok;													// return ok.
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefBase_V2::FindingBoundingIndices		2008-1-12*/
/** **/
/*---------------------------------------------------------------------------*/

size_t SKTRAN_GridDefBase_V2::FindingBoundingIndices	( double radius, ENUM_INTERPOLATIONMODE outrange, SKTRAN_GridIndex* index, double* weight, size_t maxvertices ) const
{
	double				lowerweight = 0.0;
	double				upperweight = 0.0;
	size_t				numvertices = 0;
	SKTRAN_GridIndex	lowercell = 0;
	SKTRAN_GridIndex	uppercell = 0;
	bool				ok;
	size_t				idx;

	ok = ( maxvertices >= 2 );
	ok = ok && FindBoundingIndices( radius, outrange, &lowercell, &lowerweight, &uppercell, &upperweight);
	if (!ok)								// If we failed to find the interpolant
	{										// then
		numvertices = 0;					// set the number of vertices to 0
	}										// and that is that
	else									// otherwise interpolation worked
	{										// so
		idx = 0;							// get the  number of non-zero weights
		if (lowerweight != 0.0)				// if the lower value is non zero
		{									// then
			index [idx] = lowercell;		// add the index to our array
			weight[idx] = lowerweight;		// and the lower weight to our weight array
			idx++;							// increment number of weights
		}									// and
		if (upperweight != 0.0)				// do the same for the upper weight
		{
			index[idx]  = uppercell;
			weight[idx] = upperweight;
			idx++;
		}
		numvertices = idx;					// count number of non-zero weights
	}										// and that is that
	return numvertices;						// return number of vertices found
}

/*-----------------------------------------------------------------------------
 *			SKTRAN_CoordinateTransform_V2::SetGridSearchMode_Uniform		2013-10-24*/
/** **/
/*---------------------------------------------------------------------------*/
bool SKTRAN_GridDefBase_V2::SetGridSearchMode_Uniform(  )
{
	bool ok = true;
	
	double tolerance = 1e-7;	// Allow the grid to be non-uniform by up to a tenth of a millimeter
	double delta     = 1.0;

	if(1 < m_gridvalues.size()){
		delta = At(1) - At(0);
		ok = ok && 1e-10 < delta;
		for(const_iterator it = begin()+1; it < end(); ++it)
		{
			ok = ok && fabs((*it-*(it-1)-delta)) < tolerance;
		}
	}

	if(ok){
		m_gridsearchmode = GRIDSEARCH_UNIFORM;
		m_InvUniformDelta   = 1.0 / delta;
	} else{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_GridDefBase_V2::SetGridSearchMode_Uniform, Grid is not uniform, cannot switch mode.");
	}

	return ok;
}

/*-----------------------------------------------------------------------------
 *			SKTRAN_CoordinateTransform_V2::SetGridSearchMode_NonUniform		2013-10-24*/
/** **/
/*---------------------------------------------------------------------------*/
bool SKTRAN_GridDefBase_V2::SetGridSearchMode_NonUniform(  )
{
	m_gridsearchmode = GRIDSEARCH_NONUNIFORM;	// Can always treat any grid like it's non-uniform -- no checking necessary
	return true;
}

/*-----------------------------------------------------------------------------
 *			SKTRAN_CoordinateTransform_V2::SetGridSearchMode		2013-10-24*/
/** Search can be optimized if the grid is known to be uniform. To attempt
  * optimization call this method with argument ENUM_GRIDSEARCHMODE::GRIDSEARCH_UNIFORM.
  * Method returns true if grid is uniform and optimization can be done, false if 
  * values in grid vary by more than 1e-7. **/
/*---------------------------------------------------------------------------*/
bool SKTRAN_GridDefBase_V2::SetGridSearchMode( SKTRAN_GridDefBase_V2::ENUM_GRIDSEARCHMODE mode )
{
	bool ok = true;

	switch(mode){
	case GRIDSEARCH_NONUNIFORM:
		ok = ok && SetGridSearchMode_NonUniform();
		break;
	case GRIDSEARCH_UNIFORM:
		ok = ok && SetGridSearchMode_Uniform();		// Check that the grid is uniform before configuring uniform search
		break;
	default:
		ok = false;
		break;
	}

	if(!ok) nxLog::Record(NXLOG_WARNING, "SKTRAN_GridDefBase_V2::SetGridSearchMode, Could not change mode.");
	return ok;
}

/*-----------------------------------------------------------------------------
 *			SKTRAN_CoordinateTransform_V2::LowerBound		2013-10-24*/
/** Return the first element e such that value<=e **/
/*---------------------------------------------------------------------------*/
SKTRAN_GridDefBase_V2::const_iterator SKTRAN_GridDefBase_V2::LowerBound ( const_iterator start, const_iterator finish, double value ) const
{
	const_iterator iter;
	switch(m_gridsearchmode){
	case GRIDSEARCH_NONUNIFORM:
		iter = std::lower_bound( start, finish, value );
		break;
	case GRIDSEARCH_UNIFORM:
		iter = LowerBound_Uniform( start, finish, value );
		break;
	default:
		nxLog::Record( NXLOG_ERROR, "SKTRAN_GridDefBase_V2::LowerBound, Grid search mode is invalid.");
		break;
	}
	return iter;
}

/*-----------------------------------------------------------------------------
 *			SKTRAN_CoordinateTransform_V2::UpperBound		2013-10-24*/
/** Return the first element e such that value<e **/
/*---------------------------------------------------------------------------*/
SKTRAN_GridDefBase_V2::const_iterator SKTRAN_GridDefBase_V2::UpperBound ( const_iterator start, const_iterator finish, double value ) const
{
	const_iterator iter;
	switch(m_gridsearchmode){
	case GRIDSEARCH_NONUNIFORM:
		iter = std::upper_bound( start, finish, value );
		break;
	case GRIDSEARCH_UNIFORM:
		iter = UpperBound_Uniform( start, finish, value );
		break;
	default:
		nxLog::Record( NXLOG_ERROR, "SKTRAN_GridDefBase_V2::LowerBound, Grid search mode is invalid.");
		break;
	}
	return iter;
}

/*-----------------------------------------------------------------------------
 *			SKTRAN_CoordinateTransform_V2::LowerBound_Uniform		2013-10-24*/
/** **/
/*---------------------------------------------------------------------------*/
SKTRAN_GridDefBase_V2::const_iterator SKTRAN_GridDefBase_V2::LowerBound_Uniform ( const_iterator start, const_iterator /*finish*/, double value ) const
{

	const_iterator it1, it2;
	it1 = start + (SKTRAN_GridIndex)((value-*start)*m_InvUniformDelta);
	if(*it1 < value){
		it2 = it1 + 1;
	} else{
		if(start!=it1 && value<=*(it1-1)){
			it2 = it1 - 1;
		} else{
			it2 = it1;
		}
	}

	return it2;
}

/*-----------------------------------------------------------------------------
 *			SKTRAN_CoordinateTransform_V2::UpperBound_Uniform		2013-10-24*/
/** **/
/*---------------------------------------------------------------------------*/
SKTRAN_GridDefBase_V2::const_iterator SKTRAN_GridDefBase_V2::UpperBound_Uniform ( const_iterator start, const_iterator finish, double value ) const
{

	const_iterator it1, it2;
	it1 = start + (SKTRAN_GridIndex)((value-*start)*m_InvUniformDelta);
	if( it1 >= finish )
	{
		return finish;
	}
	if(*it1 <= value)
	{
		it2 = it1 + 1;
	}
	else 
	{
		if(start!=it1 && value < *(it1-1))
		{
			it2 = it1 - 1;
		}
		else
		{
			it2 = it1;
		}
	}

	return it2;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_CoordinateTransform_V2::SKTRAN_CoordinateTransform_V2		2008-1-9*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_CoordinateTransform_V2::SKTRAN_CoordinateTransform_V2()
{
	double nan = std::numeric_limits<double>::quiet_NaN();

	m_earthRadius = nan;
	m_latitude    = nan;
	m_longitude   = nan;
	m_mjd         = nan;
	m_groundaltitude_meters = 0.0;				// Use sensible defaults for the moment
	m_toaaltitude_meters    = 100000.0;			// Use NaN in the future 

//	m_numinstances++;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_CoordinateTransform_V2::~SKTRAN_CoordinateTransform_V2		2008-1-9*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_CoordinateTransform_V2::~SKTRAN_CoordinateTransform_V2()
{
//	--m_numinstances;
//	NXTRACE(("SKTRAN_CoordinateTransform_V2::~SKTRAN_CoordinateTransform_V2, destroying object, # of instances left = %Iu\n", (size_t)m_numinstances));
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_CoordinateTransform_V2::ConfigureCoordinates		2008-1-9*/
/** This is the primary method to configure the coordinate system.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_CoordinateTransform_V2::ConfigureCoordinates( double latitude, double longitude, double mjd, const nxVector& sun )
{
	double		sza;
	nxVector	location;
	bool		ok;


	do																	//	Keep looping untl we are well away from the Sun
	{																	// for each loop
		m_truegeoid.FromGeodetic( latitude, longitude );				// Get the true location of this latitude and longitude
		location = m_truegeoid.Location();								// Get the geocentric location
		sza      = location.AngleTo(sun);								// Get the solar zenith angle
		ok       = (sza >= 0.04);										// make sure its ok
		if (!ok) latitude += 0.05;										// if not step in small steps until it is
	} while (!ok);														// Repeat until good.

	m_truegeoid.GetOsculatingSpheroid( &m_earthRadius, &m_centre );		// Get the center and radius of the fitted osculating sphere
	m_oscgeoid.SetTrueSphere(m_earthRadius);							// Define the osculating spheroid
	m_oscgeoid.FromGeodetic( latitude,longitude, 0.0 );					// Set the location of the point
	m_referencepoint_unit = m_oscgeoid.Location().UnitVector();			// Cache the unit vector towardslocation of the reference point in osculating frame
	m_latitude       = latitude;											// Cache the latitude
	m_longitude      = longitude;											// Cache the longitude
	m_mjd            = mjd;
	NXASSERT(( (sun.Magnitude() < 1.0000000001) && (sun.Magnitude() > 0.9999999999)));								// Throw exception if user throws in a position vector for Sun rather than a directional unit vector
	m_sun         = sun;												// Copy the sun over, its a directional unit vector so dont do any translation to osculating sphere center

	ok = ok && ConfigureGlobalTransform();
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_CoordinateTransform_V2::ManuallySetOsculatingSphereRadius		2014-1-31*/
/** A method that allows the advanced user to manually set the osculating
 *	radius. This is typically used to compare outputs with other codes where
 *	we need to get identical output**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_CoordinateTransform_V2::ManuallySetOsculatingSphereRadius( double radius )
{
	double deltar = radius - m_earthRadius;						// Get the change in radius
	m_oscgeoid.SetTrueSphere(radius);							// Define the osculating spheroid to a true spher of the new radius
	m_earthRadius = radius;										// Save the new radius	
	m_centre      = m_centre - deltar*m_referencepoint_unit;	// adjust offset so surface of new sphere is still at the surface of Earth 
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_CoordinateTransform_V2::ConfigureCoordinates		2008-1-11*/
/** This is an alternative method to configure the coordinate system.
 *	It simply calls the primary configuration method. The code takes a small
 *	precaution to avoid being directly underneath the Sun.  This condition
 *	causes trigonometry issues when defining the profile tables in Sasktran and
 *	is best avoided by displacing the reference point.  This will have minor
 *	impact on the SASKTRAN model, especially if multiple diffuse profiles are
 *	used.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_CoordinateTransform_V2::ConfigureCoordinates( const nxVector& observer, const nxVector& look, double mjd, const nxVector& sun )
{
	nxGeodetic		geoid;
	double			latitude;
	double			longitude;
	nxVector		location;
	double			sza;

	geoid.FromTangentPointLocation( observer, look );
	location = geoid.Location();								// Get the tangent point location
	sza      = location.AngleTo( sun );							// get the solar zenith angle
	while ( sza < 0.05 )										// While we are pretty close to the sub-solar point
	{															// we are going to move the reference point
		location = location + 5000.0*look;						// by adding 5 km in the look direction
		sza      = location.AngleTo( sun );						// update the solar zenith angle
	}															// and try again.
	latitude  = geoid.GeodeticLatitude();
	longitude = geoid.GeodeticLongitude();
	return ConfigureCoordinates( latitude, longitude, mjd, sun.UnitVector() );
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_CoordinateTransform_V2::ConfigureGlobalTransform		2007-12-20*/
/** The global transform defines a spherical system wherethe sun is on the
 *	z axis, zero degrees longitude is the great circle joining the
 *	sun and the reference point.  The global transform is the set of three
 *	vectors that define the global X, Y,Z as measured at the equator of this
 *	global system. Points in the diffuse grid are assumed to lie along lines of zero degree
 *	longitude in this grid, which maps to the great cirle joining the
 *	reference point and the sun.
 *
 *	\par Quality Control
 *	2008-02-13, ndl303: I have stepped through this code with a known sun
 *	position (1,0,0) and the reference point at 0 degrees latitude and 70 degrees
 *	longitude. It generated the correct vectors.  Note that the code will fail
 *	if the reference point is directly below the sun.  Under these conditions the
 *	code should choose a reference point offset from the sub-solar point.
 *
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_CoordinateTransform_V2::ConfigureGlobalTransform()
{
	bool		ok;
	nxVector	sun;
	nxVector	refpt;
	nxVector	zprime;
	nxVector	xprime;
	nxVector	yprime;

	zprime    = GetSunUnitVector();								// The Z axis is the direction to the sun, expressed in geographic coords
	refpt     = m_referencepoint_unit;							// Use the reference point in the spatial grid as the definition of 0 degrees longitude.
	xprime    = refpt.ComponentPerpendicularTo( zprime );		// Get the component perpendicular to Z
	ok        = xprime.Magnitude() > 1.0E-4;					// Are observer and sun parallel
	if (!ok)													// Yes,
	{															// Add code later to support this option
		nxLog::Record(NXLOG_WARNING, "SKTRAN_TableDiffusePoints_V21::ConfigureTransform, observer and sun are parallel, This is not properly supported and may cause errors");
	}
	xprime = xprime.UnitVector();								// Normalize the x prime direction
	m_heliounitvector[0] = xprime;								// The Helio X axis is perpendicular to sun in plane of reference point
	m_heliounitvector[1] = zprime ^ xprime;						// The Y axis is a third axis of right hand system
	m_heliounitvector[2] = zprime;								// The Helio Z axis is towards the sun
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_CoordinateTransform_V2::GeographicToHelio		2009-1-26*/
/** **/
/*---------------------------------------------------------------------------*/

HELIODETIC_VECTOR SKTRAN_CoordinateTransform_V2::GeographicToHelio( const nxVector& geographic ) const
{
	double					x, y, z;
	nxVector				geo;
	double					r;
	HELIODETIC_VECTOR		helio;

	r = geographic.Magnitude();

	geo.SetCoords( geographic.X()/r, geographic.Y()/r, geographic.Z()/r );

	x = m_heliounitvector[0] & geo;
	y = m_heliounitvector[1] & geo;
	z = m_heliounitvector[2] & geo;

	helio.SetCoords( x*r, y*r, z*r );
	return helio;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_CoordinateTransform_V2::GeographicToHelio		2009-1-26*/
/** **/
/*---------------------------------------------------------------------------*/

HELIODETIC_UNITVECTOR SKTRAN_CoordinateTransform_V2::GeographicToHelioUnitVector( const nxVector& geo ) const
{
	double					x, y, z;
	double					m;
	HELIODETIC_UNITVECTOR	helio;

	m = geo.Magnitude();
	m = (m > 0.0) ? 1.0/m : 0.0;
	x = m_heliounitvector[0] & geo;
	y = m_heliounitvector[1] & geo;
	z = m_heliounitvector[2] & geo;

	helio.SetCoords( m*x, m*y, m*z );
	return helio;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_CoordinateTransform_V2::GeodeticPointFromSolarCoords		2008-12-19*/
/** Get the geographic latitude and longitude of a point in our grid at the
 *	specified cosine of solar zenith angle.
 **/
/*---------------------------------------------------------------------------*/

GEODETIC_INSTANT SKTRAN_CoordinateTransform_V2::PointToGeodetic( const HELIODETIC_POINT& helio, double mjd ) const
{
	GEODETIC_INSTANT	pt;
	nxVector			v;

	v =      m_heliounitvector[0]*helio.LocalZenith().X()
		   + m_heliounitvector[1]*helio.LocalZenith().Y()
		   + m_heliounitvector[2]*helio.LocalZenith().Z();

	pt.mjd       = mjd;
	pt.heightm   = helio.Altitude();
	pt.latitude  = v.Latitude();
	pt.longitude = v.Longitude();
	return pt;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_CoordinateTransform_V2::HelioVectorToGeographic		2009-1-27*/
/** **/
/*---------------------------------------------------------------------------*/

nxVector  SKTRAN_CoordinateTransform_V2::HelioVectorToGeographic( const HELIODETIC_VECTOR& helio) const
{
	return ( m_heliounitvector[0]*helio.X()
		   + m_heliounitvector[1]*helio.Y()
		   + m_heliounitvector[2]*helio.Z());
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_CoordinateTransform_V2::HelioUnitVectorToGeographic		2009-1-27*/
/** **/
/*---------------------------------------------------------------------------*/

nxVector SKTRAN_CoordinateTransform_V2::HelioUnitVectorToGeographic( const HELIODETIC_UNITVECTOR& helio) const
{
	return ( m_heliounitvector[0]*helio.X()
		   + m_heliounitvector[1]*helio.Y()
		   + m_heliounitvector[2]*helio.Z());
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_CoordinateTransform_V2::GeographicToHeliodetic		2009-1-22*/
/** Converts a geographic position vector to Heliodetic Cooodinates**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_CoordinateTransform_V2::HelioVectorToHelioPoint( const HELIODETIC_VECTOR& geo, HELIODETIC_POINT* point) const
{
	double				r;

	r               = geo.Magnitude();
	point->Initialize( geo.UnitVector(), r, this );
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_CoordinateTransform_V2::ReferencePoint		2010-2-24*/
/** **/
/*---------------------------------------------------------------------------*/

HELIODETIC_POINT SKTRAN_CoordinateTransform_V2::ReferencePoint( double altitude_meters ) const
{
	nxVector				refgeo;
	HELIODETIC_UNITVECTOR	refv;
	HELIODETIC_POINT		location;
	double					r;

	refgeo = ReferencePointUnitVector();
	refv   = GeographicToHelio( refgeo ).UnitVector();
	r      = AltitudeToRadius(altitude_meters);
	location.Initialize( refv, r, this );
	return location;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_CoordinateTransform_V2::SetAtmosphereAltitudeBounds		 2014- 11- 28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_CoordinateTransform_V2::SetAtmosphereAltitudeBounds( double groundalt_meters, double toaalt_meters)
{
	bool	ok;

	m_groundaltitude_meters = groundalt_meters;		// Altitude of the ground in meters
	m_toaaltitude_meters    = toaalt_meters;		// Altitude of the top of the atmosphere in meters.
	
	ok = (m_toaaltitude_meters > m_groundaltitude_meters) && (m_toaaltitude_meters > 5000.0);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_CoordinateTransform_V2::SetAtmosphereAltitudeBounds, the Top of the Atmosphere altitude (%f) is less than the ground altitude or it is too small (did you use kms by accident) or is undefined", (double)m_toaaltitude_meters);
	}
	return ok;
}

