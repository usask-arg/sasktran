#include "include/sktran_hr_internals.h"



/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_LinesOfSightTable::~SKTRAN_HR_LinesOfSightTable		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_HR_LinesOfSightTable::~SKTRAN_HR_LinesOfSightTable()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_LinesOfSightTable::ReleaseResources		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_LinesOfSightTable::ReleaseResources()
{

	//for( size_t idx = 0; idx < m_geometryrays.size(); idx++ )
	//{
	//	delete m_geometryrays[idx];
	//	delete m_opticalrays[idx];
	//}
	//m_geometryrays.clear();
	m_opticalrays.clear();
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_LinesOfSightTable::SetLinesOfSight		2013-06-12*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_LinesOfSightTable::SetLinesOfSight( const SKTRAN_LineOfSightArray_V21& linesofsight, const SKTRAN_CoordinateTransform_V2& coords )
{
	bool						ok = true;

	ok = ok && m_observerlinesofsight.DeepCopy( linesofsight );
	if( !ok )
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_HR_LinesOfSightTable, There were errors copying the lines of sight");
		return ok;
	}

	// translate observer to osculating sphere
	SKTRAN_LineOfSightEntry_V2* entry;
	nxVector offsetobserver;
	for( size_t idx = 0; idx < m_observerlinesofsight.NumRays(); idx++ )
	{
		m_observerlinesofsight.GetRayVar( idx,  &entry );
		offsetobserver = coords.TranslateGeoidToOsculatingSphere( entry->Observer() );
		entry->Configure( offsetobserver, entry->Look(), entry->Mjd() );
	}


	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_LinesOfSightTable::CreateRays		2013-06-12*/
/** This takes the internal linesofsight, set by SetLinesOfSight, and 
 *  converts them to rays which can be used inside the model
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_LinesOfSightTable::CreateRays( const SKTRAN_RayFactory_Base* rayfactory )
{
	bool	ok = true;
	bool	ok1;

	const SKTRAN_LineOfSightEntry_V2*	entry;
	size_t								idx;
	size_t								numrays;
	HELIODETIC_VECTOR					observer;	
	HELIODETIC_UNITVECTOR				look;

	numrays = m_observerlinesofsight.NumRays();
	m_opticalrays.resize( numrays );
	for( idx = 0; idx < numrays; idx++ )
	{
		ok1 = rayfactory->CreateRayObject( &m_opticalrays[idx] );
		ok1 = ok1 && m_observerlinesofsight.GetRay( idx, &entry );
		if (ok1)
		{
			observer = rayfactory->CoordsPtr()->GeographicToHelio(entry->Observer());
			look     = rayfactory->CoordsPtr()->GeographicToHelio(entry->Look()).UnitVector();
			ok1      = m_opticalrays[idx]->MoveObserver( observer, look);
		}
		ok = ok && ok1;
	}

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_LinesOfSightTable::RayAt		2013-06-12*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_RayOptical_Base* SKTRAN_HR_LinesOfSightTable::RayAt( size_t idx )
{
	return m_opticalrays.at(idx).get();
}

const SKTRAN_RayOptical_Base* SKTRAN_HR_LinesOfSightTable::RayAt( size_t idx ) const
{
	return m_opticalrays.at(idx).get();
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_LinesOfSightTable::RayEntryAt		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

std::unique_ptr<SKTRAN_RayOptical_Base>&  SKTRAN_HR_LinesOfSightTable::RayEntryAt( size_t idx )
{
	return m_opticalrays.at(idx);
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_LinesOfSightTable::MeanMJD		2013-06-27*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_HR_LinesOfSightTable::MeanMJD() const
{
	double	mjd = 0;
	size_t	numrays;

	numrays = m_observerlinesofsight.NumRays();
	for (size_t idx = 0; idx < numrays; idx++)
	{
		mjd += m_observerlinesofsight.Entry(idx)->Mjd();
	}
	mjd /= numrays;
	return mjd;
}
