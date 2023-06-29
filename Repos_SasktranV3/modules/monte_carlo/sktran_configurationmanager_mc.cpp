#include "include/sktran_montecarlo_internals.h"


SKTRAN_ConfigurationManager_MC::SKTRAN_ConfigurationManager_MC()
{
}

SKTRAN_ConfigurationManager_MC::~SKTRAN_ConfigurationManager_MC()
{
	ReleaseResources();
}

void SKTRAN_ConfigurationManager_MC::ReleaseResources()
{
}

bool SKTRAN_ConfigurationManager_MC::ConfigureCoordinateTransform(const nxVector& sun, const SKTRAN_LineOfSightArray_V21& linesofsight, double surfaceHeight, double toaHeight, const std::vector<double>& reference_point, bool nadirreferencepointonground )
{
	bool ok = true;

	m_raymanager.Clear();
	m_coordinatesystem.reset( (const SKTRAN_CoordinateTransform_V2*)nullptr);
	if( sun.IsValid() && !sun.IsZero() )
	{
		ok = ok && m_raymanager.SetSun( sun );
	}
	if (reference_point.size() == 4)
	{
		ok = ok && m_raymanager.SetReferencePoint(reference_point.at(0), reference_point.at(1), reference_point.at(2), reference_point.at(3));
	}

	ok = ok && m_raymanager.SetNadirReferencePointOnGround(nadirreferencepointonground);
	ok = ok && m_raymanager.UpdateUndefinedParametersFromLinesOfSight( linesofsight );
//	nxLog::Record(NXLOG_INFO,"**** TODO **** SKTRAN_ConfigurationManager_MC::ConfigureCoordinateTransform, need to set proper mina and max altitudes for MC engine, currently set to 0 and 100000");
	ok = ok && m_raymanager.MakeCoordinateSystem( &m_coordinatesystem, surfaceHeight, toaHeight );
    //printf("\n\n  p: %17.14e %17.14e, %17.14e\n\n", m_coordinatesystem->ReferencePtLatitude(), m_coordinatesystem->ReferencePtLongitude(), m_coordinatesystem->ReferencePointMJD() );

    return ok;
}

