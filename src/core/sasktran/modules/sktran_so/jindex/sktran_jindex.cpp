#include "../sasktranv21_internals.h"
#include <limits>
#include <cmath>
#include <float.h>

/*-----------------------------------------------------------------------------
 *					SKTRANSO_JIndex::ConfigureDiffuseTableIndex		2007-12-21*/
/** Configure the JIndex so it indexes a radiance in the Diffuse Table.
 *	The diffuse table uses 3 indices to index a radiance,
 *		- sza,    indexes position, typically used for sza, ranges from 0-65535.
 *		- radius, indexes position, typically altitude, ranges from at least 0 to 255, but usually 0-65535.
 *		- vertex, indexes radiance on unit sphere, ranges from at least 0-255, but usually 0-65535.
 *	This code sets all of the elements available in the structure and is used by other helper functions
 *	to set their entries.
 **/
/*---------------------------------------------------------------------------*/

void SKTRANSO_JIndex::ConfigureDiffuseTableIndex( SKTRAN_GridIndex szaindex, SKTRAN_GridIndex radiiindex, double szaweight, double radiiweight, size_t vertexindex, double vertexweight /*, double weight*/ )
{
	double	factor;

	factor            = szaweight*radiiweight; // *weight;
	m_szagridindex    = (SKTRAN_GridIndexByte)szaindex;
	m_radiigridindex  = (SKTRAN_GridIndexShort)radiiindex;
	m_vertexindex     = (SKTRAN_GridIndexByte)vertexindex;
	m_weight          = SKTRAN_DBL_TO_GridIndexFloat( vertexweight*factor );

	NXASSERT(( szaindex       < SKTRAN_MAX_GridIndexByte ));		// Make sure we can truncate to unsigned short without difficulty
	NXASSERT(( radiiindex     < SKTRAN_MAX_GridIndexShort ));
	NXASSERT(( vertexindex    < SKTRAN_MAX_GridIndexByte ));

}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_JIndex::ConfigureEmissionTableIndex		 2015- 3- 5*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRANSO_JIndex::ConfigureEmissionTableIndex( SKTRAN_GridIndex positionindex,   SKTRAN_GridIndex heightindex, double szaweight, double altweight, bool isground)
{
	ConfigureDiffuseTableIndex( positionindex, heightindex, szaweight, altweight, isground ? 1 : 0, 1.0);
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_JIndex::IsEmissionGround		 2015- 3- 5*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_JIndex::IsEmissionGround()  const
{
	return (m_vertexindex == 1);
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_JIndex::ConfigureSolarTransmissionTableIndex		2010-2-18*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRANSO_JIndex::ConfigureSolarTransmissionTableIndex	( SKTRAN_GridIndex positionindex,   SKTRAN_GridIndex heightindex, double szaweight, double radiiweight, size_t vertexindex, double vertexweight )
{
	ConfigureDiffuseTableIndex( positionindex, heightindex, szaweight, radiiweight, vertexindex, vertexweight );
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_JIndex::ConfigureGroundPointTableIndex		2010-2-18*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRANSO_JIndex::ConfigureGroundPointTableIndex (SKTRAN_GridIndex pointindex, double pointweight, size_t vertexindex, double vertexweight )					// get the point indexing, weights and the vertex index points to the incoming r
{
	ConfigureDiffuseTableIndex( pointindex, 0, pointweight, 1.0, vertexindex, vertexweight );
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_JIndex::ConfigureScatterMatrixTableIndex		2010-2-18*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRANSO_JIndex::ConfigureScatterMatrixTableIndex( SKTRAN_GridIndex positionindex, SKTRAN_GridIndex heightindex, double szaweight, double radiiweight, size_t vertexindex, double vertexweight )
{
	ConfigureDiffuseTableIndex( positionindex, heightindex, szaweight, radiiweight, vertexindex, vertexweight );
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_JIndex::ConfigureLOSSolarTransmissionIndex		2010-6-17*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRANSO_JIndex::ConfigureLOSSolarTransmissionIndex( SKTRAN_GridIndex quadpointindex, double weight )
{
	ConfigureDiffuseTableIndex( 0, quadpointindex, 1.0, weight, 0, 1.0);
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_JIndexArray::SKTRANSO_JIndexArray		2008-2-2*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_JIndexArray::SKTRANSO_JIndexArray()
{
	m_numcells     = 0;
	m_numQ         = 0;
	m_numindices   = 0;
	m_Qindex       = NULL;
	m_firstQinCell = NULL;
	m_Jarray       = NULL;

	#ifdef NXDEBUG
		m_maxcells   = 0;
		m_maxQ       = 0;
		m_maxindices = 0;
	#endif

}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_JIndexArray::~SKTRANSO_JIndexArray		2008-2-2*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_JIndexArray::~SKTRANSO_JIndexArray()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_JIndexArray::ReleaseResources		2008-2-2*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRANSO_JIndexArray::ReleaseResources()
{
	if ( m_Qindex       != NULL) delete [] m_Qindex;
	if ( m_firstQinCell != NULL) delete [] m_firstQinCell;
	if ( m_Jarray       != NULL) delete [] m_Jarray;
	m_numcells     = 0;
	m_numQ         = 0;
	m_numindices   = 0;
	m_Qindex       = NULL;
	m_firstQinCell = NULL;
	m_Jarray       = NULL;

	#ifdef NXDEBUG
		m_maxcells   = 0;
		m_maxQ       = 0;
		m_maxindices = 0;
	#endif

}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_JIndexArray::QuadraturePointsInCell		2008-2-2*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_JIndexArray::QuadraturePointsInCell( size_t cellidx, size_t* startptindex, size_t* numquadraturepoints) const
{
	bool	ok;

	ok = (cellidx < m_numcells);
	if (ok)
	{
		*startptindex         = m_firstQinCell[cellidx++];
		*numquadraturepoints  = m_firstQinCell[cellidx] - *startptindex ;
	}
	else
	{
		*startptindex = 0;
		*numquadraturepoints = 0;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_JIndexArray::JIndicesAtQuadraturePoint		2008-2-9*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_JIndexArray::JIndicesAtQuadraturePoint( size_t quadraturepointidx, size_t* Jindex, size_t* numentries) const
{
	bool	ok;
	size_t	jindex;

	ok = (quadraturepointidx < m_numQ);
	if (ok)
	{
		jindex       = m_Qindex[ quadraturepointidx++ ];			// Get the start location of this quadrature point in the list 
		*numentries  = m_Qindex[ quadraturepointidx] - jindex;		// Get the number of radiances that contribute to evaluating the source function at this point
		*Jindex      = jindex;							// Get the starting entry for this quadrature point
	}
	else
	{
		*Jindex     = 0;
		*numentries = 0;
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_JIndexArray::JEntriesAtQuadraturePoint		2008-2-2*/
/** Fetch the array of radiance index entries that contribute to evaluating the
 *	quadrature point specified by quadraturepointidx.  It is expected that this
 *	code is called during the optical initialization phase when Sasktran must
 *	convert the source function indexes into actual pointers to source function
 *	radiances.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_JIndexArray::JEntriesAtQuadraturePoint( size_t quadraturepointidx, const SKTRANSO_JIndex** Jstart, size_t* numentries) const
{
	bool	ok;
	size_t	jindex;

	ok = JIndicesAtQuadraturePoint( quadraturepointidx, &jindex, numentries );
	if (ok) *Jstart = m_Jarray + jindex ;
	else    *Jstart = NULL;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_JIndexArray::AllocateMaximumStorage		2008-2-4*/
/** This is a class that allocates the storage specified by the caller.
 *	This is expected to be called either from the SKTRAN_JindexTableFactory
 *	object initialization which is allocating a big storage area reponsible
 *	for buffering copies of SKTRAN_JindexTable or from the DeepCopy method.
 *
 *	In eitehr circumstance the caller is fully aware of the buffer sizes and is
 *	careful not to overflow storage areas.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_JIndexArray::AllocateMaximumStorage( size_t numcells, size_t numQ, size_t numindices)
{
	bool	ok;

	ReleaseResources();

	m_numcells     = 0;
	m_numQ         = 0;
	m_numindices   = 0;
	m_Qindex       = new uint32_t        [ numQ       + 1];
	m_firstQinCell = new uint32_t        [ numcells   + 1];
	m_Jarray       = new SKTRANSO_JIndex [ numindices ];

	ok = ( m_Qindex != NULL ) && ( m_firstQinCell != NULL) && ( m_Jarray != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRANSO_JIndexArray::AllocateMaximumStorage, Error allocating space for maximum storage.");
		ReleaseResources();
	}
	else
	{
		m_Qindex[0]       = 0;				// The first quadrature point starts at index 0 in m_JArray;
		m_firstQinCell[0] = 0;				// The first cell starts at index 0 in m_Qindex;

#ifdef NXDEBUG
		size_t s = sizeof(SKTRANSO_JIndex);
		NXASSERT(( s==16  ));
		m_maxcells = numcells;
		m_maxQ     = numQ;
		m_maxindices = numindices;
#endif

	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_JIndexArray::ResetCounters		2008-2-4*/
/** This is a method provided to support SKTRAN_JIndexTableFactory.  The
 *	factory object allocates one instance of SKTRAN_JIndexTable with a large
 *	storage area. This is used to process all the rays in the geometry.
 *	This resets the counters so a new set of entries can be properly inserted 
 **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_JIndexArray::ResetCounters()
{
	m_numcells     = 0;
	m_numQ         = 0;
	m_numindices   = 0;
	return			true;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_JIndexArray::InsertQuadraturePointEntries		2008-2-4*/
/** Inserts a copy of all of the SKTRANSO_JIndex entries needed to represent
 *	the source function(J) at one quadrature point.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_JIndexArray::InsertQuadraturePointEntries( const SKTRANSO_JIndex* Jlinearsum,  size_t numentries)
{
	size_t	idx;

	for (idx = 0; idx < numentries; idx++ )						// Copy all of the entries
	{															// specified by the user
		NXASSERT(( m_numindices < m_maxindices ));
		m_Jarray[m_numindices++] = Jlinearsum[idx];				// over to our array, increment the output location
	}															// do all of the elements for this soucre function representation
	NXASSERT((m_numQ < m_maxQ));
	m_Qindex[++m_numQ] = m_numindices;							// Indicate where 

	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_JIndexArray::FinishCellEntries		2008-2-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_JIndexArray::FinishCellEntries()
{
	NXASSERT(( m_numcells < m_maxcells ));
	m_firstQinCell[ ++m_numcells] = m_numQ;						// The first Q in the next cell is 
	if (m_numQ == 0 ) ResetCounters();							// If we have no cell entries then reset the counters;
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_JIndexArray::DeepCopy		2008-2-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool  SKTRANSO_JIndexArray::DeepCopy( const SKTRANSO_JIndexArray& other )
{
	bool	ok;
	size_t	idx;

	ok = AllocateMaximumStorage( other.m_numcells, other.m_numQ, other.m_numindices );
	if (ok)
	{
		m_numcells   = other.m_numcells;
		m_numQ       = other.m_numQ;
		m_numindices = other.m_numindices;

		for (idx = 0;  idx < m_numcells + 1; idx++) m_firstQinCell[idx] = other.m_firstQinCell[idx];
		for (idx = 0;  idx < m_numQ + 1;     idx++) m_Qindex[idx]       = other.m_Qindex[idx];
		for (idx = 0;  idx < m_numindices;   idx++) m_Jarray[idx]       = other.m_Jarray[idx];
	}
	else
	{
		nxLog::Record( NXLOG_WARNING,"SKTRANSO_JIndexArray::DeepCopy, Error making a deep copy ofthe JIndex table. That is going to cause a problem or poor accuracy");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_JIndexArray::MultiplyAllWeightsBy		2008-7-18*/
/** **/
/*---------------------------------------------------------------------------*/

void  SKTRANSO_JIndexArray::MultiplyAllWeightsBy( double value )
{
	size_t	idx;

	for (idx = 0;  idx < m_numindices;   idx++)
	{
		m_Jarray[idx].MultiplyWeightBy( value);
	}
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_JIndexArray::DumpTable		2008-4-23*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRANSO_JIndexArray::DumpTable() const
{
	size_t	idx;

	for (idx = 0;  idx < m_numcells + 1; idx++) printf("Cell      %3u  first = %3u\n", (unsigned int )idx, (unsigned int)m_firstQinCell[idx]);
	for (idx = 0;  idx < m_numQ + 1;     idx++) printf("Quadpoint %3u  first = %3u\n", (unsigned int )idx, (unsigned int)m_Qindex[idx]);
	for (idx = 0;  idx < m_numindices;   idx++) printf("Index     %3u        = %3u, %3u, %3u, %24.16e\n", (unsigned int)idx, (unsigned int)m_Jarray[idx].HeightIndex(),(unsigned int)m_Jarray[idx].PositionIndex(),(unsigned int)m_Jarray[idx].VertexIndex(),(double)m_Jarray[idx].VertexWeight()  );
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_JValueTable_V21::SKTRAN_JValueTable_V21		2008-2-7*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_JValueTable_V21::SKTRAN_JValueTable_V21()
{
	m_jindex    = NULL;
//	m_radiances = NULL;
//	m_weights   = NULL;
	m_reservesize = 0;
	m_numnonzero  = 0;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_JValueTable_V21::ReleaseResources		2008-2-7*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_JValueTable_V21::ReleaseResources()
{
	m_radiances.clear();
	m_weights.clear();
	//if (m_radiances != NULL ) delete [] m_radiances;
	//if ( m_weights  != NULL ) delete [] m_weights;
	//m_radiances   = NULL;
	//m_weights     = NULL;
	m_reservesize = 0;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_JValueTable_V21::~SKTRAN_JValueTable_V21		2008-2-7*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_JValueTable_V21::~SKTRAN_JValueTable_V21()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_JValueTable_V21::ReserveStorage		2008-3-26*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_JValueTable_V21::ReserveStorage( size_t numpoints )
{
	bool ok;

	ok = (numpoints <= m_reservesize);
	if (!ok)
	{
		ReleaseResources();
		ok = (numpoints == 0 );
		if (!ok)
		{
			if      (numpoints < 10 ) numpoints += 2;	// Add a few extra points to
			else if (numpoints < 50 ) numpoints += 6;	// accomodate future requests that are
			else                      numpoints += 12;	// similar but not identical. Helps with fragmentation. 
			m_radiances.resize(numpoints); // = new const SKTRAN_StokesScalar* [numpoints];
			m_weights.resize(numpoints); //   = new SKTRAN_StokesScalar        [numpoints];
			//ok = (m_radiances != NULL) && (m_weights != NULL);
			ok = (m_weights.size() == numpoints) && (m_radiances.size() == numpoints);
			if (ok)
			{
				m_reservesize = numpoints;
#if defined(NXDEBUG)
				for (size_t i = 0; i < numpoints; i++)
				{
					m_radiances[i] = NULL;
					m_weights[i]  = std::numeric_limits<SKTRAN_StokesScalar>::quiet_NaN();
				}
#endif
			}
			else
			{
				nxLog::Record(NXLOG_WARNING,"SKTRAN_JValueTable_V21::AllocateStorage, Error allocating storage for %Iu table elements",(size_t)numpoints);
			}
		}
	}
	if (!ok) ReleaseResources();
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_JValueTable_V21::AttachToGeometry		2008-2-7*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_JValueTable_V21::AttachToGeometry( const SKTRANSO_JIndexArray& jindex )
{
	size_t						numentries;
	bool						ok;

	m_jindex   = &jindex;
	ok = (m_jindex != NULL);
	if (ok)
	{
		numentries = m_jindex->NumIndices();					// Get the number of source J that contribute to this table
		ok         = ReserveStorage( numentries );
		if (ok)
		{
			m_numnonzero = numentries;
		}
	}
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING, "SKTRAN_JValueTable_V21::AttachToGeometry, error connecting to JindexTable and allocating memory");
		ReleaseResources();
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_JValueTable_V21::Clear		2009-2-4*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_JValueTable_V21::Clear()
{
	m_numnonzero = 0;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_JValueTable_V21::JWeightsAtQuadraturePoint		2008-2-9*/
/** Get the array of weighst associated with a quadrature point in the
 *	JValueTable.
 *	This function does not seem to be used anymore ndl303, 2012-09-14 so I have commented it out
 **/
/*---------------------------------------------------------------------------*/

/*
bool  SKTRAN_JValueTable_V21::JWeightsAtQuadraturePoint	( size_t quadraturepointidx, SKTRAN_StokesScalar** weights, size_t* numentries)
{
	bool	ok;
	size_t	idx;

	ok = m_jindex->JIndicesAtQuadraturePoint( quadraturepointidx, &idx, numentries );
	if (ok) *weights = &m_weights.at(idx);
	else    *weights = NULL;
	return ok;
}
*/

/*-----------------------------------------------------------------------------
 *					SKTRAN_JValueTable_V21::AdjustWeightsAndTrim		2008-2-9*/
/** Multiply the weights of each quadrature point in this array
 *	by the value stored in the array factor. In addition the code trims out any
 *	points which are zero and make no meaningful contribution to the integrals.
 *	There should only be one call to AdjustWeightsAndTrim after each call to
 *	AttachToGeometry ...unless ConvertJIndexTable (or SetWeightsAndRadiancePtrs) is called to reset the weights
 *  and radiances -Tony Bathgate. The relationship bewteen indexes and quadrature points
 *	is messed up by this call and the user should only call Evaluate after this function
 *
 *	\param factor
 *	An input array [numquadraturepoints] of adjustment factors.  These factor
 *	are applied (multiplied) to the weights of each quadrature point in this table.
 ***/
/*---------------------------------------------------------------------------*/

bool SKTRAN_JValueTable_V21::AdjustWeightsAndTrim( const double* factor, size_t numquadraturepoints )
{
	bool		ok;
	double		f;
	size_t		qidx;
	size_t		i;
	size_t		startindex;
	size_t		numindices;
	size_t		outindex;
	double		v;		

	NXASSERT(( m_numnonzero == m_jindex->NumIndices() ));								// Make sure the user is not calling this function twice as it will be all buggered up.
	ok = (numquadraturepoints == m_jindex->NumQ());										// Make sure we have the correct number of quadrature points
	outindex = 0;
	if (ok)																				// if we do
	{																					// then 
		for (qidx =0; qidx < numquadraturepoints; qidx++)								// iterate over the quadrature points
		{																				// for each quadrature point
			f  = factor[qidx];															// get the multiplication factor for this quadrature point
			if ( f != 0.0)																// If its non zero then continue else skip
			{															
				m_jindex->JIndicesAtQuadraturePoint( qidx, &startindex, &numindices );	// Get the start of this quadrtaure point in our table
				for( i= 0; i < numindices; i++)							// and for each of the
				{														// weights in this quadrature point
					v = m_weights[startindex];							// Get the value of this weight
					if ( v != 0.0)
					{
						m_weights   [outindex] = SKTRAN_DBL_TO_GridIndexFloat( v*f );	// apply the same factor
						//Following line added 2009/08/25- Tony Bathgate
						m_radiances [outindex] = m_radiances[startindex];				// Copy the radiance field, must match weights
						outindex++;
					}
					startindex++;
				}														// do the entire linear sum of J's for this quadrature point.
			}
		}																// do all of the quadrature points
		m_numnonzero = outindex;										// store the number of non-zero elements in this table.
	}																	// and that is that
	else																// otherwise do we have a problem
	{																	// if we do then letthe user know
		nxLog::Record( NXLOG_WARNING,"SKTRAN_JValueTable_V21::AdjustWeights, Quadrature point mismatch, there has been a logic error, oops");
	}																// as things are screwed up. Time to hit the debugger.
	return ok;														// return the status
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_JValueTable_V21::AdjustWeightsByConstantFactor		2009-1-22*/
/** Multiply the weights of each every point in this array
 *	by the constant value. The relationship bewteen indexes and quadrature points
 *	is messed up by this call and the user should only call Evaluate after this function
 *  ...unless ConvertJIndexTable (or SetWeightsAndRadiancePtrs) is called to reset the weights and radiances -Tony Bathgate
**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_JValueTable_V21::AdjustWeightsByConstantFactorAndTrim( double factor )
{
	bool		ok;
	size_t		idx;
	size_t		outindex;
	double		v;		
	double		w;

	//Following line added 2009/08/25 - Tony Bathgate.
	NXASSERT(( m_numnonzero == m_jindex->NumIndices() ));						// Make sure the user is not calling this function twice as it will be all buggered up.
	ok       = (m_weights.size() > 0); // != NULL);
	outindex = 0;
	if (ok)
	{
		for (idx = 0; idx < m_numnonzero; idx++)
		{
			v              = m_weights[idx];									// Get the value of this weight as a double
			w              = v*factor;											// Get the new weight
			if (w != 0.0)														// if the new weight is not zero
			{																	// then 
				m_weights   [outindex] = SKTRAN_DBL_TO_GridIndexFloat( w );		// move the Adjust it by the factor
				m_radiances [outindex] = m_radiances[idx];					// Copy the radiance field 
				outindex++;
			}
		}
		m_numnonzero = outindex;										// store the number of non-zero elements in this table.
	}
	else
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_JValueTable_V21::AdjustWeightsByConstantFactor: Error adjusting factors by constant as weight array is  NULL");
	}
	return ok;														// return the status
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_JValueTable_V21::Evaluate		2008-2-9*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_StokesScalar	SKTRAN_JValueTable_V21::Evaluate	() const
{
	SKTRAN_StokesScalar						sum;
	size_t									idx;
	size_t									numindices;
	size_t									numloops;
	size_t									numrem;
	size_t									l;
	const SKTRAN_StokesScalar* const *		radiances;
	const SKTRAN_StokesScalar*				weights;
	const SKTRAN_StokesScalar* const *		r;
	const SKTRAN_StokesScalar*				w;


	sum        = 0;									// We want to evaluate the sum
	numindices = m_numnonzero;						// So get the number of indices that are non-zero.
	NXASSERT(( (numindices <= m_radiances.size()) && (numindices <= m_weights.size()) ));
	if (numindices > 0)
	{
		if ( (numindices < 12) )		// IF we have just a few elements
		{												// then we just
			radiances  = &m_radiances.front();
			weights    = &m_weights.front();
			for (idx = 0; idx < numindices; idx++ )				// do the summation
			{													// as
				sum +=  *(radiances[idx])*weights[idx];			// that is the quickest
			}
		}
		else
		{
			numloops   = numindices/12;				// break it into chunks of 12
			numrem     = numindices%12;				// but dont forget the remainder
			radiances  = &m_radiances.front();
			weights    = &m_weights.front();
			idx = 0;								// set the index into the radiances and weights
			for (l = 0; l < numloops; l++ )			// for the twelve point loop
			{										// loop ove rthe 12 points
				r    = radiances + idx;				// get the radiance pointer
				w    = weights   + idx;				// get the weight points
				sum +=   *(r[ 0])*w[ 0]				// add the product of the radiance times the weight
					+ *(r[ 1])*w[ 1]				// Do all twelve. This is essentially loop unrolling
					+ *(r[ 2])*w[ 2]				// and helps the compiler speed up the processor
					+ *(r[ 3])*w[ 3]				// by avoiding excessive loops.
					+ *(r[ 4])*w[ 4]				// We might even be able to use SSE/MMX instructions in this region
					+ *(r[ 5])*w[ 5]
					+ *(r[ 6])*w[ 6]
					+ *(r[ 7])*w[ 7]
					+ *(r[ 8])*w[ 8]
					+ *(r[ 9])*w[ 9]
					+ *(r[10])*w[10]
					+ *(r[11])*w[11];			// do all twelve
				idx +=  12;
			}

			for (l = 0; l < numrem; l++ )						// For the last few
			{													// just
				sum +=  *(m_radiances[idx])*m_weights[idx];		// do the calculation
				idx++;											// and that is that
			}
		}
	}
	NXASSERT(( NXFINITE(sum) ));
	return sum;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_JValueTable_V21::DumpTable		2008-4-23*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_JValueTable_V21::DumpTable() const
{
	double r;
	double w;

	for (size_t idx = 0; idx < m_numnonzero; idx++ )			// for the twelve point loop
	{										// loop ove rthe 12 points
		r    = *(m_radiances[idx]);			// get the radiance pointer
		w    = m_weights[idx];				// get the weight points
		printf("%5u %24.15e  %24.15e  \n", (unsigned int)idx, (double)r, (double)w);
	}
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_JValueTable_V21::SetWeightsAndRadiancePtrs		2009/08/25 -Tony Bathgate */
/** 
* This was added to replace part of the function ConvertJIndexTable.  It had 
* to reset m_numnonzero (a private member) but I didn't want to give the converter
* access to it directly.  This shouldn't change functionality in any way except
* it ensures that AdjustWeightsAndTrim won't fail its assertion when processing
* multiple wavelengths.
* 
**/
/*---------------------------------------------------------------------------*/
bool SKTRAN_JValueTable_V21::SetWeightsAndRadiancePtrs( const SKTRANSO_JindexTableBase* table,  ENUM_SKTRAN_JSOURCE jsource  )
{
	bool								ok;
	size_t								idx;
	const SKTRAN_StokesScalar*			radianceptr;
	const SKTRANSO_JIndex*				entry;

	NXASSERT( m_jindex!=NULL );
	ok =  table!=NULL;

	for (idx = 0; idx < m_jindex->NumIndices(); idx++)
	{
		entry				=  m_jindex->At(idx);
		radianceptr			= table->ConvertJIndexToRadiancePtr( entry, jsource );
		ok					= ok && radianceptr!=NULL;
		m_radiances[idx]	= radianceptr;
		m_weights[idx]		=  entry->VertexWeight();
	}

	NXASSERT( ok );
	if( ok )
	{
		m_numnonzero		= m_jindex->NumIndices();
	}

	return ok;
}

