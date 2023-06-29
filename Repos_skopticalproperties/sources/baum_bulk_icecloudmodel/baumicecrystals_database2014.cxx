#include <skopticalproperties21.h>

#include "nxnetcdfio.h"

/**	Set the default truncation angles for all De effective sizes.  Selections are 
 *	based on minimizing the absolute geometric mean (over inc. directions) of the 
 *	photon-conservation multipliers.
 *	For reference on these selections, see Wiensz et al, JQSRT, 2012. 
 */
static double g_cutoff_DE   [] = { 0.0, 40.0, 50.0, 10000.0};			// The De at which we change the automatic forward scatter cuttoff, follows Truitts paper/thesis
static double g_cutoff_angle[] = { 5.0,  5.0,  2.0,     2.0};			// The corresponding cutoff angles


/*-----------------------------------------------------------------------------
 *					skBaumIceCrystals_DataBase::CurrentIndex::CurrentIndex		2014-4-14*/
/** **/
/*---------------------------------------------------------------------------*/

skBaumIceCrystals_DataBase::CurrentIndex::CurrentIndex(const nx1dArray<double>&	values)
							: m_valuearray(values)
{
	m_index_0      = 0;
	m_index_1      = 0;
	m_value_0      = std::numeric_limits<double>::quiet_NaN();
	m_value_1      = std::numeric_limits<double>::quiet_NaN();
	m_currentvalue = std::numeric_limits<double>::quiet_NaN(); 
}

/*-----------------------------------------------------------------------------
 *					skBaumIceCrystals_DataBase::CurrentIndex::UpdateIndices		2014-4-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool skBaumIceCrystals_DataBase::CurrentIndex::UpdateIndices(double v)
{
	bool	ok;

	ok = (v == m_currentvalue);
	if (!ok)
	{
		ok = ( v >= m_value_0 ) && ( v <= m_value_1 );						// See if the current range is still valid
		if (!ok)															// If it is not
		{																	// Then 
			ok = (v >= m_valuearray.front() && v <= m_valuearray.back());	// make sure we are in raneg of the array, otehrwise lookup the bounding indices
			ok = ok && nxLinearInterpolate::FindBoundingIndicesAscending(m_valuearray.begin(), m_valuearray.end(), v, &m_index_0, &m_index_1, &m_value_0, &m_value_1);
		}																	// and we are done
		m_currentvalue = v;													// This is the new current value
		if (!ok)															// If we had problems then we are out of range.
		{																	// so do some diagnostics and fix it up.
			if ( v < m_valuearray.front())
			{
				m_index_0 = 0;
				m_index_1 = 0;
			}
			if ( v > m_valuearray.back() )
			{
				m_index_0 = m_valuearray.size() -1;
				m_index_1 = m_valuearray.size() -1;
			}
			nxLog::Record(NXLOG_WARNING,"skBaumIceCrystals_DataBase::CurrentIndex::UpdateIndices, Value %g was out of range for the array, %g to %g", (double)v, (double)m_valuearray.front(), (double)m_valuearray.back() );
			m_currentvalue = std::numeric_limits<double>::quiet_NaN();
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skBaumIceCrystals_DataBase::skBaumIceCrystals_DataBase		2014-4-14*/
/** **/
/*---------------------------------------------------------------------------*/

skBaumIceCrystals_DataBase::skBaumIceCrystals_DataBase()
	                       : m_config( "USask-ARG","SkOpticalProperties/Baum_IceCrystals/Storage/",nxRegistryConfiguration::GLOBAL_INI, true),
						     m_Wavelenindex( m_wavelen ),
							 m_Deindex     ( m_De ),
							 m_PhaseAngleindex( m_phaseangles)
{
	m_thetaCutoff     = 5.0;
	m_thetacutoff_auto = true;												// Use an automatic theta cut off angle
	m_icecrystalshape = BAUM_GENERAL_HABIT_MIXTURE_SEVERELYROUGH;			// Baum 2014 report. Figure 9, General Habit gives good agreement in the Solar band.
	m_isdirty         = true;
}


/*-----------------------------------------------------------------------------
 *					skBaumIceCrystals_DataBase::~skBaumIceCrystals_DataBase		2014-4-14*/
/** **/
/*---------------------------------------------------------------------------*/

skBaumIceCrystals_DataBase::~skBaumIceCrystals_DataBase()
{
}


/*-----------------------------------------------------------------------------
 *					skBaumIceCrystals_DataBase::FetchBaseDirectory		2014-4-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool skBaumIceCrystals_DataBase::FetchBaseDirectory( )
{
	bool	ok;

	ok = !m_basedir.IsEmpty();
	if (!ok)
	{
		ok = m_config.LocateDirectoryFromKey("2014Database", &m_basedir, false, false, "Browse for the directory containing the BAUM 2014 database files.");
		m_basedir.EnsureLastCharIsDirectoryChar();
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skBaumIceCrystals_DataBase::FetchFilename		2014-4-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool skBaumIceCrystals_DataBase::FetchFilename( nxString* filename )
{
	nxString	name;
	bool		ok = true;

	ok = FetchBaseDirectory();
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skBaumIceCrystals_DataBase::FetchFilename, Cannot load the base directory from the registry");
		*filename = "bad registry";
	}
	else
	{
		switch (m_icecrystalshape)
		{
		case BAUM_AGGREGATE_SOLID_COLUMNS_SEVERELYROUGH :	name = "AggregateSolidColumns_SeverelyRough_AllWavelengths_FullPhaseMatrix"; break;
		case BAUM_GENERAL_HABIT_MIXTURE_SEVERELYROUGH :		name = "GeneralHabitMixture_SeverelyRough_AllWavelengths_FullPhaseMatrix"; break;
		case BAUM_SOLID_COLUMNS_SEVERELYROUGH :				name = "SolidColumns_SeverelyRough_AllWavelengths_FullPhaseMatrix"; break;
		default :
			ok = false;
			nxLog::Record(NXLOG_WARNING, "skBaumIceCrystals_DataBase::FetchFilename, Unknown ICE CRYSTAL shape");
			*filename = "unknonwn ice crystal shape";
			break;
		};
	}
	if (ok)
	{
		// Check to see if the Legendre cached file exists
		if (nxDirectory::FileExists(m_basedir + name + "_LegendreAdded.nc"))
		{
			m_legendrecached = true;
			*filename = m_basedir + name + "_LegendreAdded.nc";
		}
		else
		{
			m_legendrecached = false;
			*filename = m_basedir + name + ".nc";
			ok = nxDirectory::FileExists(*filename);
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "skBaumIceCrystals_DataBase::FetchFilename, file does not exist <%s>", (const char*)(*filename));
			}
		}


	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skBaumIceCrystals_DataBase::SetForwardScatterCutoffAngle		2014-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skBaumIceCrystals_DataBase::SetForwardScatterCutoffAngle( double cutoff_degrees )
{
	m_thetaCutoff      = cutoff_degrees;
	m_thetacutoff_auto = ( 0.0 > cutoff_degrees );

	return true;
}

/*-----------------------------------------------------------------------------
 *					skBaumIceCrystals_DataBase::LoadDatabase		2014-4-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool skBaumIceCrystals_DataBase::LoadDatabase( bool useDeltaEddington )
{
	static std::mutex				g_mutex_xsectionslock;

	nxString						filename;
	nxNetcdfFile					file;
	bool							ok;
	nx2dArray<double>				totalarea;
	nx2dArray<double>				singlescatalbedo;
	nx2dArray<double>				ext_eff;
	nx2dArray<double>				ext_over_iwc;
	double							nan = std::numeric_limits<double>::quiet_NaN();
//	nx2dArray<double>				iwc;
//	nx2dArray<double>				asymparam;

	ok = !m_isdirty;
	if (!ok)
	{
		std::unique_lock<std::mutex>	lock(g_mutex_xsectionslock);			// Lock this thread
		ok = !m_isdirty;															// Make sure we are still dirty once we get the thread back
		if (!ok)																	// if it is dirty then lets do the netcdf file I/O stuff.
		{
			ok = FetchFilename( &filename);
			ok = ok && file.OpenRead(filename);
			ok = ok && file.VarAt( "wavelengths"                     )->LoadData( &m_wavelen,        nan);
			ok = ok && file.VarAt( "effective_diameter"              )->LoadData( &m_De,             nan);
		//	ok = ok && file.VarAt( "ice_water_content"               )->LoadData( &iwc,              nan);
			ok = ok && file.VarAt( "total_area"                      )->LoadData( &totalarea,        nan);
		//	ok = ok && file.VarAt( "asymmetry_parameter"             )->LoadData( &asymparam,        nan);
			ok = ok && file.VarAt( "single_scattering_albedo"        )->LoadData( &singlescatalbedo, nan);
			ok = ok && file.VarAt( "extinction_efficiency"           )->LoadData( &ext_eff,          nan);
			ok = ok && file.VarAt( "extinction_coefficient_over_iwc" )->LoadData( &ext_over_iwc,     nan);
			// p21/p22/p33/p43/p44 in the Baum database file are actually relative to p11.
			// We scale them in InterpolateP11 because it would be more accurate to interpolate these
			// relative values
			ok = ok && file.VarAt( "p11_phase_function"              )->LoadData( &m_p11,            nan);
			ok = ok && file.VarAt( "p21_phase_function"              )->LoadData( &m_p21,            nan);
			ok = ok && file.VarAt( "p22_phase_function"              )->LoadData( &m_p22,            nan);
			ok = ok && file.VarAt( "p33_phase_function"              )->LoadData( &m_p33,            nan);
			ok = ok && file.VarAt( "p43_phase_function"              )->LoadData( &m_p43,            nan); 
			ok = ok && file.VarAt( "p44_phase_function"              )->LoadData( &m_p44,            nan);
			ok = ok && file.VarAt( "phase_angles"                    )->LoadData( &m_phaseangles,    nan);
			if (m_legendrecached)
			{
				ok = ok && file.VarAt("phase_moments")->LoadData(&m_phasemoments, nan);
				ok = ok && file.VarAt("p11_phase_moments")->LoadData(&m_p11_legendre, nan);
			}
			file.Close();
			if (ok)
			{
				m_wavelen              *= 1000.0;														// Convert microns to nanometers
				m_ext_crosssection      = totalarea * ext_eff;											// Extinction cross-section = total area * extinction efficiency (microns**2)
				m_ext_crosssection     *=  1.0E-08;														// Convert squalre microns to square cms.
				m_scatt_crosssection    =  (m_ext_crosssection * singlescatalbedo);						// The scattering cross-section = extinction cross-section times single scattering albedo
				m_abs_crosssection      =  (m_ext_crosssection - m_scatt_crosssection);					// Get the absorption cross-section
				if( useDeltaEddington ){
					ok = TruncateAndComputeDeltaFraction( m_p11, &m_deltafraction);	// Compute forward scatter for each setting of wavelength (w) and effective size De.
				} else{
					ok = TruncateAndComputeDeltaFraction_NoTruncation( m_p11, &m_deltafraction ); // Set delta fraction to zero
				}

			}
			m_isdirty = !ok;
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "skBaumIceCrystals_DataBase::LoadDatabase, There were errors loading the Baum database file <%s>", (const char*) filename);
			}
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skBaumIceCrystals_DataBase::InterpolateCrossSectionArray		2014-4-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool skBaumIceCrystals_DataBase::InterpolateCrossSectionArray( const nx2dArray<double>& xs, double* value )
{
	double	v[4];

	v[0] = xs.At( m_Wavelenindex.Idx0(),  m_Deindex.Idx0() ); 		//-# v[0] = v(x0, y0)
	v[1] = xs.At( m_Wavelenindex.Idx0(),  m_Deindex.Idx1() );		//-# v[1] = v(x0, y1)
	v[2] = xs.At( m_Wavelenindex.Idx1(),  m_Deindex.Idx1() );		//-# v[2] = v(x1, y1)
	v[3] = xs.At( m_Wavelenindex.Idx1(),  m_Deindex.Idx0() );		//-# v[3] = v(x1, y0)

	*value = nxLinearInterpolate::FromSquare( m_Wavelenindex.V(), m_Deindex.V(),
											m_Wavelenindex.V0(), m_Wavelenindex.V1(),
											m_Deindex.V0(),      m_Deindex.V1(), 
											v);
	return true;
}

bool skBaumIceCrystals_DataBase::InterpolatePhaseArrayAtIndex(const nx3dArray<double>& xs, size_t phaseidx, double* value)
{
	double	v[4];

	v[0] = xs.At(m_Wavelenindex.Idx0(), m_Deindex.Idx0(), phaseidx); 		//-# v[0] = v(x0, y0)
	v[1] = xs.At(m_Wavelenindex.Idx0(), m_Deindex.Idx1(), phaseidx);		//-# v[1] = v(x0, y1)
	v[2] = xs.At(m_Wavelenindex.Idx1(), m_Deindex.Idx1(), phaseidx);		//-# v[2] = v(x1, y1)
	v[3] = xs.At(m_Wavelenindex.Idx1(), m_Deindex.Idx0(), phaseidx);		//-# v[3] = v(x1, y0)

	*value = nxLinearInterpolate::FromSquare(m_Wavelenindex.V(), m_Deindex.V(),
		m_Wavelenindex.V0(), m_Wavelenindex.V1(),
		m_Deindex.V0(), m_Deindex.V1(),
		v);
	return true;

}

/*-----------------------------------------------------------------------------
 *					skBaumIceCrystals_DataBase::InterpolateCrossSection		2014-4-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool skBaumIceCrystals_DataBase::InterpolateCrossSections( double w, double de, double* absxs, double* extxs, double* scattxs)
{
	bool ok;

	ok  = (de == 0.0);
	if (ok)
	{
		*absxs = 0.0;
		*extxs = 0.0;
		*scattxs = 0.0;
	}
	else
	{
		ok  =       m_Wavelenindex.UpdateIndices( w );
		ok  = ok && m_Deindex.UpdateIndices( de );
		ok  = ok && InterpolateCrossSectionArray( m_abs_crosssection,   absxs);
		ok  = ok && InterpolateCrossSectionArray( m_ext_crosssection,   extxs);
		ok  = ok && InterpolateCrossSectionArray( m_scatt_crosssection, scattxs);
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skBaumIceCrystals_DataBase::InterpolateCrossSections, errors interpolating Baum 2014 Ice Crystal cross-sections for wavelength %g nm, Effective diameter (De) = %g microns", (double)w, (double)de);
		*absxs = 0.0;
		*extxs = 0.0;
		*scattxs = 0.0;
	}
	return ok;	
}


/*-----------------------------------------------------------------------------
 *					skBaumIceCrystals_DataBase::InterpolateForwardScatter		2014-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skBaumIceCrystals_DataBase::InterpolateForwardScatter( double w, double de, double* forwardscatter)
{
	bool ok;

	ok  = (de == 0.0);
	if (ok)
	{
		*forwardscatter= 0.0;
	}
	else
	{
		ok  =       m_Wavelenindex.UpdateIndices( w );
		ok  = ok && m_Deindex.UpdateIndices( de );
		ok  = ok && InterpolateCrossSectionArray( m_deltafraction,   forwardscatter);
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skBaumIceCrystals_DataBase::InterpolateForwardScatter, errors interpolating Baum 2014 Ice Crystal forward scatter for wavelength %g nm, Effective diameter (De) = %g microns", (double)w, (double)de);
		*forwardscatter = 0.0;
	}
	return ok;	
}


/*-----------------------------------------------------------------------------
 *					skBaumIceCrystals_DataBase::InterpolateP11		2014-4-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool skBaumIceCrystals_DataBase::InterpolateP11( double w, double de, double angle, nx3dArray<double>& p, double* p11value)
{
	bool					ok;
	nx2dArray<double>		p11lower;
	nx2dArray<double>		p11upper;
	double					p11[2];

	ok  =       m_Wavelenindex.UpdateIndices( w );							
	ok  = ok && m_Deindex.UpdateIndices( de );
	ok  = ok && m_PhaseAngleindex.UpdateIndices( angle );

	ok  = ok && InterpolatePhaseArrayAtIndex( p, m_PhaseAngleindex.Idx0(), &p11[0]);
	ok  = ok && InterpolatePhaseArrayAtIndex( p, m_PhaseAngleindex.Idx1(), &p11[1] );
	if (ok)
	{
		*p11value = nxLinearInterpolate::FromTwoPoints( angle, m_PhaseAngleindex.V0(), m_PhaseAngleindex.V1(), p11);
	}
	else
	{
		*p11value = 0.0;
		nxLog::Record(NXLOG_WARNING,"skBaumIceCrystals_DataBase::InterpolateP11, Error interpolating Baum 2014 Ice Crystal phase matrix arrays for wavelength %g nm, Effective Diameter (De) = %g microns, scatter angle = %g degrees.", (double)w, (double)de, (double)angle);
	}
	return ok;
}

bool skBaumIceCrystals_DataBase::InterpolateLegendre(double w, double de, std::vector<double>& legendremoments) 
{
	bool ok;

	ok = m_Wavelenindex.UpdateIndices(w);
	ok = ok && m_Deindex.UpdateIndices(de);

	legendremoments.resize(m_phasemoments.size());
	for (size_t idx = 0; idx < legendremoments.size(); idx++)
	{
		double	v[4];

		v[0] = m_p11_legendre.At(m_Wavelenindex.Idx0(), m_Deindex.Idx0(), idx);
		v[1] = m_p11_legendre.At(m_Wavelenindex.Idx0(), m_Deindex.Idx1(), idx);
		v[2] = m_p11_legendre.At(m_Wavelenindex.Idx1(), m_Deindex.Idx1(), idx);
		v[3] = m_p11_legendre.At(m_Wavelenindex.Idx1(), m_Deindex.Idx0(), idx);

		legendremoments[idx] = nxLinearInterpolate::FromSquare(m_Wavelenindex.V(), m_Deindex.V(),
			m_Wavelenindex.V0(), m_Wavelenindex.V1(),
			m_Deindex.V0(), m_Deindex.V1(),
			v);
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skBaumIceCrystals_DataBase::InterpolatePhaseMatrix		2014-4-15*/
/** Interpolates all elements of the phase matrix to the requested wavelength
 *	effective size and scattering angle. Its a fairly large computations
 **/
/*---------------------------------------------------------------------------*/

bool skBaumIceCrystals_DataBase::InterpolatePhaseMatrix( double w, double de, double angle, skRTPhaseMatrix* phasematrix)
{
	double	p11 = 0.0;
	double	p21 = 0.0;
	double	p22 = 0.0;
	double	p33 = 0.0;
	double	p43 = 0.0;
	double	p44 = 0.0;
	bool	ok;

	ok = (de == 0.0);						// An effective size of 0.0 means "ignore"
	if (!ok)									// so everything is set to 0.0
	{																// so
		ok =       InterpolateP11( w, de, angle, m_p11, &p11);		// interpolate all elements of the phasematrix
		ok = ok && InterpolateP11( w, de, angle, m_p21, &p21);		// that are in the Baum database
		ok = ok && InterpolateP11( w, de, angle, m_p22, &p22);
		ok = ok && InterpolateP11( w, de, angle, m_p33, &p33);
		ok = ok && InterpolateP11( w, de, angle, m_p43, &p43);
		ok = ok && InterpolateP11( w, de, angle, m_p44, &p44);
	}

	if (!ok)														// and if we had problems
	{																// then set the elements to zero
		p11 = 0.0;													// and record an error message
		p21 = 0.0;
		p22 = 0.0;
		p33 = 0.0;
		p43 = 0.0;
		p44 = 0.0;
		nxLog::Record(NXLOG_WARNING,"skBaumIceCrystals_DataBase::InterpolatePhaseMatrix, error interpolating phase matrix to wavelength = %g, De = %g, angle = %g", (double)w, (double)de, (double) angle);
	}

	// Every element except p11 in the baum database is relative to p11
	// The signs on p21 and p43 are adjusted to match the Mie code
	phasematrix->At(1,1) = p11;
	phasematrix->At(2,1) = -1.0 * p21 * p11;
	phasematrix->At(3,1) = 0.0;
	phasematrix->At(4,1) = 0.0;
	phasematrix->At(1,2) = -1.0 *p21 * p11;
	phasematrix->At(2,2) = p22 * p11;
	phasematrix->At(3,2) = 0.0;
	phasematrix->At(4,2) = 0.0;
	phasematrix->At(1,3) = 0.0;
	phasematrix->At(2,3) = 0.0;
	phasematrix->At(3,3) = p33 * p11;
	phasematrix->At(4,3) = p43 * p11;
	phasematrix->At(1,4) = 0.0;
	phasematrix->At(2,4) = 0.0;
	phasematrix->At(3,4) = -1.0 * p43 * p11;
	phasematrix->At(4,4) = p44 * p11;
	return ok;
}


bool skBaumIceCrystals_DataBase::InterpolatePhaseScalar(double w, double de, double angle, double& p11)
{
	bool	ok;

	ok = (de == 0.0);						// An effective size of 0.0 means "ignore"
	if (!ok)									// so everything is set to 0.0
	{																// so
		ok = InterpolateP11(w, de, angle, m_p11, &p11);		// interpolate all elements of the phasematrix
	}

	if (!ok)														// and if we had problems
	{																// then set the elements to zero
		p11 = 0.0;													// and record an error message
		nxLog::Record(NXLOG_WARNING, "skBaumIceCrystals_DataBase::InterpolatePhaseScalar, error interpolating phase matrix to wavelength = %g, De = %g, angle = %g", (double)w, (double)de, (double)angle);
	}

	return ok;
}

/*-----------------------------------------------------------------------------
 *					skBaumIceCrystals_DataBase::ForwardScatterCutoffAngle		2014-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

double skBaumIceCrystals_DataBase::ForwardScatterCutoffAngle( double de )
{
	double cutoff = 0.0;

	if (!m_thetacutoff_auto) cutoff = m_thetaCutoff;
	else                     cutoff = nxLinearInterpolate::EvaluateYatX( de, g_cutoff_DE, g_cutoff_angle, N_ELEMENTS(g_cutoff_DE), nxLinearInterpolate::ENUM_TRUNCATE, 0.0);
	return cutoff;
}


/*-----------------------------------------------------------------------------
 *		skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::TruncateAndComputeDeltaFraction		2011-7-20
 *	Truncate phase function using a 'cutoff angle' criterion and return the directly
 *	forward-scattered fraction, 'f'.
 *	NOTE: this function overwrites the original phase function with its truncated
 *	version.
 *---------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------
 *					skBaumIceCrystals_DataBase::TruncateAndComputeDeltaFraction		2014-4-15*/
/** Computes the **/
/*---------------------------------------------------------------------------*/

bool skBaumIceCrystals_DataBase::TruncateAndComputeDeltaFraction( nx3dArray<double>& p11, nx2dArray<double>* deltafraction   )
{
	size_t				idx;
	size_t				idxr;
	size_t				idxcut;
	double				pcut;
	nx1dArray<double>	xmu;
	nx1dArray<double>	pdelta;
	nx1dArray<double>	phasefunct;
	nxSpline			spline;
	bool				ok;
	bool				ok1;
	double				cutoff;


	ok = deltafraction->SetSize( p11.XSize(), p11.YSize() );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skBaumIceCrystals_DataBase::TruncateAndComputeDeltaFraction, Error allocating space for deltafraction array (%z by %z)", (size_t)p11.XSize(), (size_t)p11.YSize() );
	}
	else
	{
		for (size_t ide =0; ide < m_p11.YSize(); ide++)
		{
			for (size_t iw = 0; iw < p11.XSize(); iw++)
			{
				ok1 = p11.Slice( iw, iw, ide, ide, 0, NXARRAY_STARSELECT, &phasefunct);
				if (ok1)
				{
					cutoff = ForwardScatterCutoffAngle( m_De.At(ide) );
					idxcut =  std::lower_bound( m_phaseangles.begin(), m_phaseangles.end(), cutoff ) - m_phaseangles.begin();	// Get the index of the cutoff angle
					if ( idxcut > 1)																							// Make sure we have at least two points for the spline
					{																											// If we do
						xmu.SetSize( idxcut );																					// then allocate the storage for the cos(angles)
						pdelta.SetSize( idxcut );																				// and for the low angle phase matrix terms
						pcut = phasefunct.At( idxcut );																			// Get the truncated phase matrix value
						for(idx=0; idx < idxcut; idx++)																			// get cos(theta) and delta-funct for angles in fwd peak
						{																										// put in order of incr. cos(theta)
							idxr				= idxcut - idx - 1;
							xmu.At(idx)			= nxmath::cosd( m_phaseangles.At( idxr ) );
							pdelta.At(idx)		= phasefunct.At( idxr ) - pcut;													// delta-funct component
							phasefunct.At(idxr) = pcut;																			// replace orig. with truncated version
						}
						ok1 = spline.Configure( xmu,pdelta );
						if (ok1 ) deltafraction->At( iw, ide ) = 0.5*spline.Integrate( NULL );
						else      deltafraction->At( iw, ide ) = 0.0;
					}
					else
					{
						deltafraction->At( iw, ide ) = 0.0;
					}
				}
				ok = ok && ok1;
			}
		}
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skBaumIceCrystals_DataBase::TruncateAndComputeDeltaFraction, Error computing forward scatter fraction, spline integral may be badly setup");
		}
	}
	return ok;
}


bool skBaumIceCrystals_DataBase::TruncateAndComputeDeltaFraction_NoTruncation( nx3dArray<double>& p11, nx2dArray<double>* deltafraction   )
{
	bool ok;
	ok = deltafraction->SetSize( p11.XSize(), p11.YSize() );
	deltafraction->SetTo( 0.0 );

	return ok;
}



