#include "include/sktran_montecarlo_internals.h"
//#include <sasktranv3.h>

SKTRAN_TableOpticalProperties_MCBase::SKTRAN_TableOpticalProperties_MCBase(){
	//m_scatterCdfMatrices = NULL;
	//m_cdfLookupSpace     = NULL;
	//m_lookupSpaceIndices = NULL;
	m_numAngles          = 0;
}

SKTRAN_TableOpticalProperties_MCBase::~SKTRAN_TableOpticalProperties_MCBase(){
	releaseTable();
}

void SKTRAN_TableOpticalProperties_MCBase::releaseTable(){
	
	//delete[]  m_scatterCdfMatrices;
	//m_scatterCdfMatrices = NULL;

};

bool SKTRAN_TableOpticalProperties_MCBase::allocateCdf(const SKTRAN_Specifications_Base *specs){
	bool ok = true;

	//const SKTRAN_SpecsInternal_V21* specscast = dynamic_cast<const SKTRAN_SpecsInternal_V21*>(specs);
	//ok = ok && NULL!=specscast;
	//ok = ok && allocateCdf( *specscast->DiffuseSpecs()->ScatterAngleGrid() , *specscast->OpticalTableSpecs()->OpticalPropertiesGrid() );
	ok = false;

	return ok;
}

bool SKTRAN_TableOpticalProperties_MCBase::allocateCdf( const SKTRAN_GridDefScatterAngle_V21& scatteranglegrid, const SKTRAN_GridDefOpticalPropertiesRadii_V21& altitudegrid )
{
	bool ok = true;

	size_t numAlts   = altitudegrid.NumAltitudes();
	//size_t numAngles = scatteranglegrid.NumAngles();
	m_numAngles = scatteranglegrid.NumAngles();
	m_numalt = numAlts;
	m_numloc = 1;

	ok = ok && 0<numAlts;
	ok = ok && 0<m_numAngles;

	releaseTable();
	m_scatterCdfMatrices.resize( numAlts* m_numAngles );
	//ok = ok && NULL!=m_scatterCdfMatrices;
	if(!ok) releaseTable();

	return ok;
}

bool SKTRAN_TableOpticalProperties_MCBase::allocateCdf( const SKTRAN_GridDefScatterAngle_V21& scatteranglegrid, const SKTRAN_GridDefOpticalPropertiesRadii_V21& altitudegrid, const SKTRAN_UnitSphere_V2& unitsphere )
{
	bool ok = true;

	size_t numAlts   = altitudegrid.NumAltitudes();
	size_t numloc    = unitsphere.NumUnitVectors();
	//size_t numAngles = scatteranglegrid.NumAngles();
	m_numAngles = scatteranglegrid.NumAngles();
	m_numalt = numAlts;
	m_numloc = unitsphere.NumUnitVectors();

	ok = ok && 0<numAlts;
	ok = ok && 0<m_numAngles;

	releaseTable();
	m_scatterCdfMatrices.resize( numloc * numAlts * m_numAngles );
	//ok = ok && NULL!=m_scatterCdfMatrices;
	if(!ok) releaseTable();

	return ok;
}

bool SKTRAN_TableOpticalProperties_MCBase::allocateCdf(const SKTRAN_GridDefScatterAngle_V21& scatteranglegrid, const SKTRAN_GridDefOpticalPropertiesRadii_V21& altitudegrid, const SKTRAN_UnitSphere_V2& unitsphere, const SKTRAN_GridDefWavelength& wavelengthgrid )
{
	bool ok = true;
	
	if (wavelengthgrid.NumWavelengths() > 0)
	{


		size_t numAlts = altitudegrid.NumAltitudes();
		size_t numloc = unitsphere.NumUnitVectors();
		size_t numwav = wavelengthgrid.NumWavelengths();
		//size_t numAngles = scatteranglegrid.NumAngles();
		m_numAngles = scatteranglegrid.NumAngles();
		m_numalt = numAlts;
		m_numloc = numloc;

		ok = ok && 0 < numAlts;
		ok = ok && 0 < m_numAngles;
		ok = ok && 0 < numwav;

		releaseTable();
		m_scatterCdfMatrices.resize(numloc * numAlts * m_numAngles * numwav);
		//ok = ok && NULL!=m_scatterCdfMatrices;
		if (!ok) releaseTable();

		return ok;
	}
	else
	{
		ok = ok && allocateCdf(scatteranglegrid, altitudegrid, unitsphere);
	}
	return ok;
}

inline bool SKTRAN_TableOpticalProperties_MCBase::makeScatterCdfs( std::unique_ptr< SKTRAN_PolarizationProperties_Base >& scatProps, size_t numangles ){
	return makeScatterCdfs_trapz( scatProps, numangles );
}

bool SKTRAN_TableOpticalProperties_MCBase::makeScatterCdfs_trapz( std::unique_ptr< SKTRAN_PolarizationProperties_Base >& scatProps, size_t numangles ){
	
	bool ok = true;
	double oldVal, newVal, temp;
    size_t pidx = 0;

	std::vector<SKTRAN_PhaseMatrixScalar>::iterator cdfit = m_scatterCdfMatrices.begin(); 

	size_t numCases = scatProps->NumElements() / numangles;
	ok = scatProps->NumElements() == numCases*numangles;
	if(ok){
		for( size_t cidx=0; cidx<numCases; ++cidx ){
			*cdfit = 0.0;
			newVal = scatProps->PhaseMatrixAccess( pidx ); 
			++pidx;
			for(size_t angleid=1; angleid<numangles; angleid++){
				oldVal = newVal;
				newVal = scatProps->PhaseMatrixAccess( pidx );
				++pidx;
				temp = *cdfit + (newVal+oldVal)/2.0;
				*(++cdfit) = temp;
			}
			++cdfit;
		}
		ok = ok && scatProps->NumElements() == pidx; // Sanity check
		if(!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_TableOpticalProperties_MCBase::makeScatterCdfs_trapz, Failed to fully traverse all phase functions!");
	} else{
		nxLog::Record(NXLOG_ERROR, "SKTRAN_TableOpticalProperties_MCBase::makeScatterCdfs_trapz, Configuration is bad, numAngles doesn't divide size of singleScatt.");
	}

	return ok;
}


bool SKTRAN_TableOpticalProperties_MCBase::FindScatterAngleBoundingIndices( std::vector< SKTRAN_PhaseMatrixScalar >::const_iterator cdf, double target, size_t& lowIndex, size_t& uppIndex, double& lowerw, double& upperw) const{

	size_t low, mid, high;

	low = 0;
	high = m_numAngles - 1;

    target *= cdf[high];
	
	while(low < high){
		mid = (low + high) / 2;
		if( target <= cdf[mid] ){
			high = mid - 1;
		} else{
			low = mid + 1;
		}
	}

	lowIndex = high;
	if((target < cdf[lowIndex]) && (lowIndex != 0)) lowIndex--;

	if( (m_numAngles-1) == lowIndex ){
		// upper end; no room for uppIndex>lowIndex
		uppIndex = lowIndex;
		lowerw = 1.0;
		upperw = 0.0;
	} else{
		// there is room in pdf for (uppIndex = lowIndex+1)
		uppIndex = lowIndex+1;
		lowerw = (cdf[uppIndex] - target) / (cdf[uppIndex] - cdf[lowIndex]);
		upperw = (target - cdf[lowIndex]) / (cdf[uppIndex] - cdf[lowIndex]);
	}

	return true;

}


bool SKTRAN_TableOpticalProperties_MCBase::GetCosScattAngleCdf(const SKTRAN_GridIndex &lowercell, const double &lowerw, const SKTRAN_GridIndex &uppercell, const double &upperw, std::vector<SKTRAN_PhaseMatrixScalar>::iterator cdfLookupSpace ) const {

	bool				ok = true;
	
	std::vector<SKTRAN_PhaseMatrixScalar>::const_iterator  lookupSpaceEnd   = cdfLookupSpace + m_numAngles;

	//if(NULL == cdfLookupSpace){
	//	nxLog::Record( NXLOG_ERROR, "SKTRAN_TableOpticalProperties_MCBase::GetCosScattAngleCdf, user must provide an array in which to store the cos(scattAngle) cdf");	// This isn't really the correct message to give
	//	ok = false;
	//}

	auto lowerit = m_scatterCdfMatrices.cbegin() + lowercell*m_numAngles;
	auto upperit = m_scatterCdfMatrices.cbegin() + uppercell*m_numAngles;

	*cdfLookupSpace = 0.0;
	++cdfLookupSpace;
	for(; cdfLookupSpace < lookupSpaceEnd; cdfLookupSpace++){
		*cdfLookupSpace = *(++lowerit)*lowerw + *(++upperit)*upperw;
	}

	return ok;

}

bool SKTRAN_TableOpticalProperties_MCBase::GetCosScattAngleCdf( size_t* altindex, double* altweight, size_t numalt, size_t* locindex, double* locweight, size_t numloc, std::vector<SKTRAN_PhaseMatrixScalar>::iterator cdf ) const
{
	bool				ok = true;
	
	std::vector< SKTRAN_PhaseMatrixScalar >::const_iterator  lookupSpaceEnd   = cdf + m_numAngles;

	//if(NULL == cdf){
	//	nxLog::Record( NXLOG_ERROR, "SKTRAN_TableOpticalProperties_MCBase::GetCosScattAngleCdf, user must provide an array in which to store the cos(scattAngle) cdf");	// This isn't really the correct message to give
	//	ok = false;
	//}

	cdf[0] = 0;
	for(size_t angleidx = 1; angleidx < m_numAngles; angleidx++)
	{
		cdf[angleidx] = 0;
		for( size_t altidx = 0; altidx < numalt; altidx++ )
		{
			for( size_t locidx = 0; locidx < numloc; locidx++ )
			{
				cdf[angleidx] += m_scatterCdfMatrices[altindex[altidx]*m_numAngles + locindex[locidx]*m_numAngles*m_numalt + angleidx] * altweight[altidx] * locweight[locidx];
			}
		}
	}

	return ok;	
}

bool SKTRAN_TableOpticalProperties_MCBase::GetCosScattAngleCdf(size_t* wavindex, double* wavweight, size_t numwav, size_t* altindex, double* altweight, size_t numalt, size_t* locindex, double* locweight, size_t numloc, std::vector<SKTRAN_PhaseMatrixScalar>::iterator cdf) const
{
	bool				ok = true;

	std::vector< SKTRAN_PhaseMatrixScalar >::const_iterator  lookupSpaceEnd = cdf + m_numAngles;

	//if(NULL == cdf){
	//	nxLog::Record( NXLOG_ERROR, "SKTRAN_TableOpticalProperties_MCBase::GetCosScattAngleCdf, user must provide an array in which to store the cos(scattAngle) cdf");	// This isn't really the correct message to give
	//	ok = false;
	//}

	cdf[0] = 0;
	for (size_t angleidx = 1; angleidx < m_numAngles; angleidx++)
	{
		cdf[angleidx] = 0;
		for (size_t wavidx = 0; wavidx < numwav; wavidx++)
		{
			for (size_t altidx = 0; altidx < numalt; altidx++)
			{
				for (size_t locidx = 0; locidx < numloc; locidx++)
				{
					cdf[angleidx] += m_scatterCdfMatrices[altindex[altidx] * m_numAngles + locindex[locidx] * m_numAngles*m_numalt + wavindex[wavidx] * m_numAngles*m_numalt*m_numloc + angleidx] * wavweight[wavidx] * altweight[altidx] * locweight[locidx];
				}
			}
		}
	}

	return ok;
}


bool SKTRAN_TableOpticalProperties_MCBase::GetCosScatteringAngleWeights( const double& randNum, SKTRAN_GridIndex& loindex, double& loweight, SKTRAN_GridIndex& hiindex, double& hiweight ) const {

	bool ok = true;
	double target;
	const size_t threadid = omp_get_thread_num();

	ok = ok && GetCosScattAngleCdf(loindex, loweight, hiindex, hiweight, m_lookupSpaceIndices[threadid] );
	target = randNum * (m_lookupSpaceIndices[threadid][m_numAngles-1]);

	// Should just have number of angles stored
	ok = ok && FindScatterAngleBoundingIndices( m_lookupSpaceIndices[threadid], target, loindex, hiindex, loweight, hiweight);

	return ok;
}

bool SKTRAN_TableOpticalProperties_MCBase::GetCosScatteringAngleWeights( double randNum, size_t* altindex, double* altweight, size_t numalt, size_t* locindex, double* locweight, size_t numloc, size_t* scatindex, double* scatweight ) const
{
	bool ok = true;
	const size_t threadid = omp_get_thread_num();

	ok = ok && GetCosScattAngleCdf( altindex, altweight, numalt, locindex, locweight, numloc, m_lookupSpaceIndices[threadid] );

	ok = ok && FindScatterAngleBoundingIndices( m_lookupSpaceIndices[threadid], randNum, scatindex[0], scatindex[1], scatweight[0], scatweight[1] );

	return ok;
}

bool SKTRAN_TableOpticalProperties_MCBase::GetCosScatteringAngleWeights(double randNum, size_t* wavindex, double* wavweight, size_t numwav, size_t* altindex, double* altweight, size_t numalt, size_t* locindex, double* locweight, size_t numloc, size_t* scatindex, double* scatweight) const
{
	bool ok = true;
	const size_t threadid = omp_get_thread_num();

	ok = ok && GetCosScattAngleCdf(wavindex, wavweight, numwav, altindex, altweight, numalt, locindex, locweight, numloc, m_lookupSpaceIndices[threadid]);

	ok = ok && FindScatterAngleBoundingIndices(m_lookupSpaceIndices[threadid], randNum, scatindex[0], scatindex[1], scatweight[0], scatweight[1]);

	return ok;
}



/* 1D code */

SKTRAN_TableOpticalProperties_1D_Height_MC::SKTRAN_TableOpticalProperties_1D_Height_MC( )
{
	//m_cdfLookupSpace     = NULL;
	//m_lookupSpaceIndices = NULL;
}

SKTRAN_TableOpticalProperties_1D_Height_MC::~SKTRAN_TableOpticalProperties_1D_Height_MC( )
{
	ReleaseResources();
}

void SKTRAN_TableOpticalProperties_1D_Height_MC::ReleaseResources( )
{
	//delete[] m_cdfLookupSpace;     m_cdfLookupSpace     = NULL;
	//delete[] m_lookupSpaceIndices; m_lookupSpaceIndices = NULL;
}

bool SKTRAN_TableOpticalProperties_1D_Height_MC::ConfigureGeometry( const SKTRAN_Specifications_Base* specs){
	bool ok = true;

	const SKTRAN_Specifications_MC* mcspecs = dynamic_cast<const SKTRAN_Specifications_MC*> (specs);
	ok = ok && NULL!=mcspecs;

	ok = ok && SKTRAN_TableOpticalProperties_1D_Height_V3::ConfigureGeometry( *mcspecs->GetScatterAngleGrid(), *mcspecs->GetOpticalPropRadii());
	ok = ok && allocateCdf( *mcspecs->GetScatterAngleGrid(), *mcspecs->GetOpticalPropRadii() );

	return ok;
}


bool SKTRAN_TableOpticalProperties_1D_Height_MC::ConfigureOptical( double wavelen, SKTRAN_AtmosphericOpticalState_V21& opticalState )
{

	bool ok = true;
	
	ok = ok && SKTRAN_TableOpticalProperties_1D_Height_V3::ConfigureOptical( wavelen, opticalState );
	ok = ok && makeScatterCdfs( m_scatprops, m_scatteranglegrid->NumAngles() );

	return ok;
}


bool SKTRAN_TableOpticalProperties_1D_Height_MC::GetCosScatteringAngle( const HELIODETIC_POINT& point, const double& randNum, double &cosScatAngle, SKTRAN_ScatMat_MIMSNC* pmatrix) const
{

	bool ok = true;

	SKTRAN_GridIndex loindex, hiindex;
	double loweight, hiweight;

	ok = ok && m_altitudegrid->FindBoundingIndices( point.Altitude(), SKTRAN_GridDefBase_V2::OUTOFBOUND_ERROR, &loindex, &loweight, &hiindex, &hiweight );
	ok = ok && GetCosScatteringAngleWeights( randNum, loindex, loweight, hiindex, hiweight );
	cosScatAngle = loweight*m_scatteranglegrid->At(loindex) + hiweight*m_scatteranglegrid->At(hiindex);

	if(nullptr!=pmatrix){
        ok = ok && GetScatteringMatrixCM2( point, cosScatAngle, *pmatrix );
        pmatrix->Normalize();
    }
    
	return ok;
}

bool SKTRAN_TableOpticalProperties_1D_Height_MC::AllocateCdfLookupSpace( size_t numThreads )
{
	bool ok = true;
	m_cdfLookupSpace.resize( numThreads * m_scatteranglegrid->NumAngles() );
	m_lookupSpaceIndices.resize( numThreads );

	//ok = ok && NULL!=m_cdfLookupSpace;
	//ok = ok && NULL!=m_lookupSpaceIndices;
	if(ok)
	{
		for(size_t tidx=0; tidx<numThreads; tidx++)
		{
			m_lookupSpaceIndices[tidx] = m_cdfLookupSpace.begin() + (tidx*m_numAngles);			// First element in this thread's lookup space
		}
	} else{
		ReleaseResources();
		nxLog::Record(NXLOG_ERROR, "SKTRAN_TableOpticalProperties_1D_Height_MC::AllocateCdfLookupSpace, Could not allocate space to perform CDF lookup.");
	}

	return ok;
}

bool SKTRAN_TableOpticalProperties_1D_Height_MC::MakeThreadsafeFor( size_t numThreads )
{
	bool ok = true;

	ok = ok && AllocateCdfLookupSpace( numThreads );
	return ok;
}

/* 3D code */
SKTRAN_TableOpticalProperties_3D_UnitSphere_MC::SKTRAN_TableOpticalProperties_3D_UnitSphere_MC( )
{
	//m_cdfLookupSpace     = NULL;
	//m_lookupSpaceIndices = NULL;
}

SKTRAN_TableOpticalProperties_3D_UnitSphere_MC::~SKTRAN_TableOpticalProperties_3D_UnitSphere_MC( )
{
	ReleaseResources();
}

void SKTRAN_TableOpticalProperties_3D_UnitSphere_MC::ReleaseResources( )
{
	//delete[] m_cdfLookupSpace;     m_cdfLookupSpace     = NULL;
	//delete[] m_lookupSpaceIndices; m_lookupSpaceIndices = NULL;
}

bool SKTRAN_TableOpticalProperties_3D_UnitSphere_MC::ConfigureGeometry( const SKTRAN_Specifications_Base* specs){
	bool ok = true;

	const SKTRAN_Specifications_MC* mcspecs = dynamic_cast<const SKTRAN_Specifications_MC*> (specs);
	ok = ok && NULL!=mcspecs;

	ok = ok && SKTRAN_TableOpticalProperties_3D_UnitSphere::SetAltitudes( *mcspecs->GetOpticalPropRadii() );
	ok = ok && SKTRAN_TableOpticalProperties_3D_UnitSphere::SetScatterGrid( *mcspecs->GetScatterAngleGrid() );
	ok = ok && SKTRAN_TableOpticalProperties_3D_UnitSphere::SetUnitSphere( *mcspecs->GetOpticalUnitSphere() );
	ok = ok && SKTRAN_TableOpticalProperties_3D_UnitSphere::SetWavelengthGrid( *mcspecs->GetWavelengthGridOpticalProperties() );
	ok = ok && SKTRAN_TableOpticalProperties_3D_UnitSphere::ConfigureGeometry( nullptr );

	ok = ok && allocateCdf( *mcspecs->GetScatterAngleGrid(), *mcspecs->GetOpticalPropRadii(), *mcspecs->GetOpticalUnitSphere(), *mcspecs->GetWavelengthGridOpticalProperties() );
	ok = ok && m_inelasticOptProp->ConfigureGeometry( this );

	return ok;
}


bool SKTRAN_TableOpticalProperties_3D_UnitSphere_MC::ConfigureOptical( double wavelen, SKTRAN_AtmosphericOpticalState_V21& opticalState ){

	bool ok = true;
	
	ok = ok && SKTRAN_TableOpticalProperties_3D_UnitSphere::ConfigureOptical( wavelen, opticalState );
	ok = ok && makeScatterCdfs( m_scatprops, m_scatteranglegrid->NumAngles() );
	ok = ok && m_inelasticOptProp->ConfigureOptical(wavelen, opticalState);
	//nxLog::Record( NXLOG_WARNING, "SKTRAN_TableOpticalProperties_3D_UnitSphere_MC::ConfigureOptical, If using precache this only fills cdf for wavidx=0." );

	return ok;
}


bool SKTRAN_TableOpticalProperties_3D_UnitSphere_MC::GetCosScatteringAngle( const HELIODETIC_POINT& point, const double& randNum, double &cosScatAngle, SKTRAN_ScatMat_MIMSNC* pmatrix) const
{

	bool ok;

	size_t numalt, numloc;
	size_t altindex[2];
	size_t locindex[3];
	size_t scatindex[2];
	double altweights[2];
	double locweights[3];
	double scatweights[2];
	ok = CalcAltIndices( point, altweights, altindex, numalt );
	ok = CalcSphereIndices( point, locweights, locindex, numloc );

	ok = ok && GetCosScatteringAngleWeights( randNum, altindex, altweights, numalt, locindex, locweights, numloc, scatindex, scatweights );
	cosScatAngle = scatweights[0]*m_scatteranglegrid->At(scatindex[0]) + scatweights[1]*m_scatteranglegrid->At(scatindex[1]);

	if(nullptr!=pmatrix){
        ok = ok && GetScatteringMatrixCM2( point, cosScatAngle, *pmatrix );
        pmatrix->Normalize();
    }

	return ok;
}

bool SKTRAN_TableOpticalProperties_3D_UnitSphere_MC::GetCosScatteringAngle(const double& wavelength, const HELIODETIC_POINT& point, const double& randNum, double &cosScatAngle, SKTRAN_ScatMat_MIMSNC* pmatrix) const
{

	bool ok;

	size_t numwav, numalt, numloc;
	size_t wavindex[2];
	size_t altindex[2];
	size_t locindex[3];
	size_t scatindex[2];
	double wavweights[2];
	double altweights[2];
	double locweights[3];
	double scatweights[2];
	ok = CalcWavelengthIndices(wavelength, wavweights, wavindex, numwav);
	ok = CalcAltIndices(point, altweights, altindex, numalt);
	ok = CalcSphereIndices(point, locweights, locindex, numloc);
	

	ok = ok && GetCosScatteringAngleWeights(randNum, wavindex, wavweights, numwav, altindex, altweights, numalt, locindex, locweights, numloc, scatindex, scatweights);
	cosScatAngle = scatweights[0] * m_scatteranglegrid->At(scatindex[0]) + scatweights[1] * m_scatteranglegrid->At(scatindex[1]);

	if (nullptr != pmatrix) {
		ok = ok && GetScatteringMatrixCM2(point, cosScatAngle, *pmatrix);
		pmatrix->Normalize();
	}

	return ok;
}

bool SKTRAN_TableOpticalProperties_3D_UnitSphere_MC::AllocateCdfLookupSpace( size_t numThreads )
{
	bool ok = true;
	m_cdfLookupSpace.resize( numThreads * m_scatteranglegrid->NumAngles() );
	m_lookupSpaceIndices.resize( numThreads );

	//ok = ok && NULL!=m_cdfLookupSpace;
	//ok = ok && NULL!=m_lookupSpaceIndices;
	if(ok)
	{
		for(size_t tidx=0; tidx<numThreads; tidx++)
		{
			m_lookupSpaceIndices[tidx] = m_cdfLookupSpace.begin() + (tidx*m_numAngles);			// First element in this thread's lookup space
		}
	} else{
		ReleaseResources();
		nxLog::Record(NXLOG_ERROR, "SKTRAN_TableOpticalProperties_1D_Height_MC::AllocateCdfLookupSpace, Could not allocate space to perform CDF lookup.");
	}

	return ok;
}

bool SKTRAN_TableOpticalProperties_3D_UnitSphere_MC::MakeThreadsafeFor( size_t numThreads )
{
	bool ok = true;

	ok = ok && AllocateCdfLookupSpace( numThreads );
	return ok;
}



/* 3D constant code */
SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant_MC::SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant_MC()
{
	//m_cdfLookupSpace     = NULL;
	//m_lookupSpaceIndices = NULL;
}

SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant_MC::~SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant_MC()
{
	ReleaseResources();
}

void SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant_MC::ReleaseResources()
{
	//delete[] m_cdfLookupSpace;     m_cdfLookupSpace     = NULL;
	//delete[] m_lookupSpaceIndices; m_lookupSpaceIndices = NULL;
}

bool SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant_MC::ConfigureGeometry(const SKTRAN_Specifications_Base* specs) {
	bool ok = true;

	const SKTRAN_Specifications_MC* mcspecs = dynamic_cast<const SKTRAN_Specifications_MC*> (specs);
	ok = ok && NULL != mcspecs;

	ok = ok && SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant::SetAltitudes(*mcspecs->GetOpticalPropRadii());
	ok = ok && SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant::SetScatterGrid(*mcspecs->GetScatterAngleGrid());
	ok = ok && SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant::SetUnitSphere(*mcspecs->GetOpticalUnitSphere());
	ok = ok && SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant::SetWavelengthGrid(*mcspecs->GetWavelengthGridOpticalProperties());
	ok = ok && SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant::ConfigureGeometry( nullptr );

	ok = ok && allocateCdf(*mcspecs->GetScatterAngleGrid(), *mcspecs->GetOpticalPropRadii(), *mcspecs->GetOpticalUnitSphere(), *mcspecs->GetWavelengthGridOpticalProperties());
	ok = ok && m_inelasticOptProp->ConfigureGeometry(this);

	return ok;
}

bool SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant_MC::ConfigureOptical(double wavelen, SKTRAN_AtmosphericOpticalState_V21& opticalState) {

	bool ok = true;

	ok = ok && SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant::ConfigureOptical(wavelen, opticalState);
	ok = ok && makeScatterCdfs(m_scatprops, m_scatteranglegrid->NumAngles());
	ok = ok && m_inelasticOptProp->ConfigureOptical(wavelen, opticalState);
	//nxLog::Record( NXLOG_WARNING, "SKTRAN_TableOpticalProperties_3D_UnitSphere_MC::ConfigureOptical, If using precache this only fills cdf for wavidx=0." );

	return ok;
}


bool SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant_MC::GetCosScatteringAngle(const HELIODETIC_POINT& point, const double& randNum, double &cosScatAngle, SKTRAN_ScatMat_MIMSNC* pmatrix) const
{

	bool ok;

	size_t numalt, numloc;
	size_t altindex[2];
	size_t locindex[3];
	size_t scatindex[2];
	double altweights[2];
	double locweights[3];
	double scatweights[2];
	ok = CalcAltIndices(point, altweights, altindex, numalt);
	ok = CalcSphereIndices(point, locweights, locindex, numloc);

	ok = ok && GetCosScatteringAngleWeights(randNum, altindex, altweights, numalt, locindex, locweights, numloc, scatindex, scatweights);
	cosScatAngle = scatweights[0] * m_scatteranglegrid->At(scatindex[0]) + scatweights[1] * m_scatteranglegrid->At(scatindex[1]);

	if (nullptr != pmatrix) {
		ok = ok && GetScatteringMatrixCM2(point, cosScatAngle, *pmatrix);
		pmatrix->Normalize();
	}

	return ok;
}

bool SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant_MC::GetCosScatteringAngle(const double& wavelength, const HELIODETIC_POINT& point, const double& randNum, double &cosScatAngle, SKTRAN_ScatMat_MIMSNC* pmatrix) const
{

	bool ok;

	size_t numwav, numalt, numloc;
	size_t wavindex[2];
	size_t altindex[2];
	size_t locindex[3];
	size_t scatindex[2];
	double wavweights[2];
	double altweights[2];
	double locweights[3];
	double scatweights[2];
	ok = CalcWavelengthIndices(wavelength, wavweights, wavindex, numwav);
	ok = CalcAltIndices(point, altweights, altindex, numalt);
	ok = CalcSphereIndices(point, locweights, locindex, numloc);


	ok = ok && GetCosScatteringAngleWeights(randNum, wavindex, wavweights, numwav, altindex, altweights, numalt, locindex, locweights, numloc, scatindex, scatweights);
	cosScatAngle = scatweights[0] * m_scatteranglegrid->At(scatindex[0]) + scatweights[1] * m_scatteranglegrid->At(scatindex[1]);

	if (nullptr != pmatrix) {
		ok = ok && GetScatteringMatrixCM2(point, cosScatAngle, *pmatrix);
		pmatrix->Normalize();
	}

	return ok;
}

bool SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant_MC::AllocateCdfLookupSpace(size_t numThreads)
{
	bool ok = true;
	m_cdfLookupSpace.resize(numThreads * m_scatteranglegrid->NumAngles());
	m_lookupSpaceIndices.resize(numThreads);

	//ok = ok && NULL!=m_cdfLookupSpace;
	//ok = ok && NULL!=m_lookupSpaceIndices;
	if (ok)
	{
		for (size_t tidx = 0; tidx < numThreads; tidx++)
		{
			m_lookupSpaceIndices[tidx] = m_cdfLookupSpace.begin() + (tidx*m_numAngles);			// First element in this thread's lookup space
		}
	}
	else {
		ReleaseResources();
		nxLog::Record(NXLOG_ERROR, "SKTRAN_TableOpticalProperties_1D_Height_MC::AllocateCdfLookupSpace, Could not allocate space to perform CDF lookup.");
	}

	return ok;
}

bool SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant_MC::MakeThreadsafeFor(size_t numThreads)
{
	bool ok = true;

	ok = ok && AllocateCdfLookupSpace(numThreads);
	return ok;
}



//
///***********
// * polarized
// ************/
//
//SKTRAN_TableOpticalProperties_1D_Height_Polarized_MC::SKTRAN_TableOpticalProperties_1D_Height_Polarized_MC( )
//{
//	m_cdfLookupSpace     = NULL;
//	m_lookupSpaceIndices = NULL;
//}
//
//SKTRAN_TableOpticalProperties_1D_Height_Polarized_MC::~SKTRAN_TableOpticalProperties_1D_Height_Polarized_MC( )
//{
//	ReleaseResources();
//}
//
//void SKTRAN_TableOpticalProperties_1D_Height_Polarized_MC::ReleaseResources( )
//{
//	delete[] m_cdfLookupSpace;     m_cdfLookupSpace     = NULL;
//	delete[] m_lookupSpaceIndices; m_lookupSpaceIndices = NULL;
//}
//
//bool SKTRAN_TableOpticalProperties_1D_Height_Polarized_MC::ConfigureGeometry( const SKTRAN_Specifications_Base* specs){
//	bool ok = true;
//
//	const SKTRAN_Specifications_MC* mcspecs = dynamic_cast<const SKTRAN_Specifications_MC*> (specs);
//	ok = ok && NULL!=mcspecs;
//
//	ok = ok && SKTRAN_TableOpticalProperties_1D_Height_Polarized_V3::ConfigureGeometry( *mcspecs->GetScatterAngleGrid(), *mcspecs->GetOpticalPropRadii());
//	ok = ok && allocateCdf( *mcspecs->GetScatterAngleGrid(), *mcspecs->GetOpticalPropRadii() );
//
//	return ok;
//}
//
//
//bool SKTRAN_TableOpticalProperties_1D_Height_Polarized_MC::ConfigureOptical( double wavelen, SKTRAN_AtmosphericOpticalState_V21& opticalState ){
//
//	bool ok = true;
//	
//	ok = ok && SKTRAN_TableOpticalProperties_1D_Height_Polarized_V3::ConfigureOptical( wavelen, opticalState );
//	ok = ok && makeScatterCdfs( m_singleScatt, m_scatteranglegrid->NumAngles() );
//
//	return ok;
//}
//
//
//bool SKTRAN_TableOpticalProperties_1D_Height_Polarized_MC::GetCosScatteringAngle( const HELIODETIC_POINT& point, const double& randNum, double &cosScatAngle, skRTPhaseMatrix* pmatrix) const
//{
//
//	bool ok;
//
//	SKTRAN_GridIndex loIndex, hiIndex;
//	SKTRAN_GridIndex loAltIndex, hiAltIndex;
//	double loAltWeight, hiAltWeight, loWeight, hiWeight;
//	std::vector<skRTPhaseMatrix>::const_iterator lowerAltZeroAngle;
//	std::vector<skRTPhaseMatrix>::const_iterator upperAltZeroAngle;
//
//	ok =       m_altitudegrid->FindBoundingIndices( point.Altitude(), SKTRAN_GridDefBase_V2::OUTOFBOUND_ERROR, &loIndex, &loWeight, &hiIndex, &hiWeight );
//	loAltIndex = loIndex;
//	hiAltIndex = hiIndex;
//	loAltWeight = loWeight;
//	hiAltWeight = hiWeight; 
//	lowerAltZeroAngle = m_muellerMatrices.begin() + loAltIndex*m_numAngles;
//	upperAltZeroAngle = m_muellerMatrices.begin() + hiAltIndex*m_numAngles;
//
//	ok = ok && GetCosScatteringAngleWeights( randNum, loIndex, loWeight, hiIndex, hiWeight );
//	cosScatAngle = loWeight*m_scatteranglegrid->At(loIndex) + hiWeight*m_scatteranglegrid->At(hiIndex);
//
//	// Photon really scatters from -oldDir to -newDir, but the angle is the same between oldDir and newDir
//	*pmatrix = (lowerAltZeroAngle[loIndex]*(loAltWeight*loWeight) + lowerAltZeroAngle[hiIndex]*(loAltWeight*hiWeight)) + (upperAltZeroAngle[loIndex]*(hiAltWeight*loWeight) + upperAltZeroAngle[hiIndex]*(hiAltWeight*hiWeight));
//	double norm = 1e-40<pmatrix->At(1,1) ? 1.0/pmatrix->At(1,1) : 0.0;
//	*pmatrix *= norm;
//
//	//if(NULL!=pmatrix)
//	//{
//	//	loIndex = loAltIndex;
//	//	hiIndex = hiAltIndex;
//	//	m_scatteranglegrid->FindBoundingIndices( -cosScatAngle, SKTRAN_GridDefBase_V2::OUTOFBOUND_ERROR, &loIndex, &loWeight, &hiIndex, &hiWeight ); // Need scatter from source to observer -- forward MC
//	//	*pmatrix = (lowerAltZeroAngle[loIndex]*(loAltWeight*loWeight) + lowerAltZeroAngle[hiIndex]*(loAltWeight*hiWeight)) + (upperAltZeroAngle[loIndex]*(hiAltWeight*loWeight) + upperAltZeroAngle[hiIndex]*(hiAltWeight*hiWeight));
//	//}
//
//	return ok;
//}
//
//bool SKTRAN_TableOpticalProperties_1D_Height_Polarized_MC::AllocateCdfLookupSpace( size_t numThreads )
//{
//	bool ok = true;
//	m_cdfLookupSpace     = new SKTRAN_PhaseMatrixScalar [numThreads * m_scatteranglegrid->NumAngles()];
//	m_lookupSpaceIndices = new SKTRAN_PhaseMatrixScalar*[numThreads];
//
//	ok = ok && NULL!=m_cdfLookupSpace;
//	ok = ok && NULL!=m_lookupSpaceIndices;
//	if(ok)
//	{
//		for(size_t tidx=0; tidx<numThreads; tidx++)
//		{
//			m_lookupSpaceIndices[tidx] = m_cdfLookupSpace + (tidx*m_numAngles);			// First element in this thread's lookup space
//		}
//	} else{
//		ReleaseResources();
//		nxLog::Record(NXLOG_ERROR, "SKTRAN_TableOpticalProperties_1D_Height_MC::AllocateCdfLookupSpace, Could not allocate space to perform CDF lookup.");
//	}
//
//	return ok;
//}
//
//bool SKTRAN_TableOpticalProperties_1D_Height_Polarized_MC::MakeThreadsafeFor( size_t numThreads )
//{
//	bool ok = true;
//
//	ok = ok && AllocateCdfLookupSpace( numThreads );
//	return ok;
//}
