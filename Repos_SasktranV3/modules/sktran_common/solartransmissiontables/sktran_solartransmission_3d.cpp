#include "../sktran_common.h"

void SKTRAN_SolarTransmission_3D::ReleaseResources()
{
	m_transmission.clear();
}

inline size_t SKTRAN_SolarTransmission_3D::sub2ind  ( size_t szaidx, size_t slonidx, size_t altidx ) const
{
	return ( ( (            szaidx
		   * m_numslon)  + slonidx )
		   * m_numalts ) +  altidx;
}

bool SKTRAN_SolarTransmission_3D::SlonWeights( double slon, double* slonweights, size_t* slonindex, size_t& numindex ) const
{
	bool ok = true;
	ok = ok && m_slongrid.FindBoundingIndices( slon, SKTRAN_GridDefBase_V2::OUTOFBOUND_TRUNCATE, &slonindex[0], &slonweights[0], &slonindex[1], &slonweights[1] );
	if( 1e-10 < slonweights[1]){
		numindex = 2;
	} else{
		numindex = 1;
	}
	return ok;
}

bool SKTRAN_SolarTransmission_3D::FillTableAtIndex( size_t szaidx, size_t slonidx, size_t altidx) const
{
	bool ok = true;

	double cossza = CosSZAGrid().CosineSZA( szaidx );
	double sinsza = sqrt(1 - cossza*cossza);

	double slon = m_slongrid.SLON( slonidx );
	double alt  = HeightGrid().At( altidx );
	double transmission;

	HELIODETIC_UNITVECTOR locunit;
	HELIODETIC_VECTOR     loc;

	locunit.SetCoords( sinsza * cos( slon ), sinsza * sin( slon ), cossza );
	loc.SetCoords( locunit, CoordinatesPtr()->AltitudeToRadius( alt ) );;
	ok = ok && CreateRayAndCalcTransmission( loc,transmission );
	m_transmission[ sub2ind( szaidx, slonidx, altidx ) ] = transmission;

	return ok;
}

std::unique_ptr<SKTRAN_RayOptical_Base> SKTRAN_SolarTransmission_3D::CreateRayAndCalcFullTransmission( const HELIODETIC_VECTOR& loc ) const
{
	bool ok;

	HELIODETIC_UNITVECTOR						sun;
	std::unique_ptr<SKTRAN_RayOptical_Base>		ray;

	sun.SetCoords(0,0,-1);

	ok =        RayFactory()->CreateRayObject( &ray);
	ok = ok &&  ray->MoveObserver( loc, sun );
	ok = ok &&  ray->TraceRay_NewMethod();
	if(ok)
	{
		Integrator()->CalculateRayScalarTransmission_withMinContainer( ray.get(), NULL, false, false );
	} 
	else
	{
		nxLog::Record(NXLOG_ERROR, "SKTRAN_SolarTransmission_3D::CreateRayAndCalcFullTransmission, Couldn't calculate transmission.");
	}

	return ray;
}


SKTRAN_SolarTransmission_3D::SKTRAN_SolarTransmission_3D(bool prefill, bool accountforrefraction)
{
	m_prefilltable = prefill;
	m_refractionenabled = accountforrefraction;
}
SKTRAN_SolarTransmission_3D::~SKTRAN_SolarTransmission_3D()
{
	ReleaseResources();
}

bool SKTRAN_SolarTransmission_3D::TransmissionAtPoint(	const HELIODETIC_POINT&					point,
													    double&									transmission) const
{
	bool ok = true;

	double cossza = point.CosSZA();
	double alt    = point.Altitude();
	double slon   = atan2( point.LongitudeY(), point.LongitudeX() );

	return ok = ok && Transmission_Interpolate( cossza, alt, slon, transmission );
}


bool SKTRAN_SolarTransmission_3D::DeflectionAtPoint(	const HELIODETIC_POINT&					point,
													    double&									deflection) const
{
	bool ok = true;

	double cossza = point.CosSZA();
	double alt    = point.Altitude();
	double slon   = atan2( point.LongitudeY(), point.LongitudeX() );

	return ok = ok && Deflection_Interpolate( cossza, alt, slon, deflection );
}

double SKTRAN_SolarTransmission_3D::CosAngleToSource( const HELIODETIC_UNITVECTOR& look, const HELIODETIC_POINT* location) const
{
	double deflection;
	bool ok = true;

	if (m_refractionenabled)
	{
		ok = ok && DeflectionAtPoint(*location, deflection);
	}
	else
	{
		deflection = 0.0;
	}

	// Nominal solar direction is (0, 0, 1), with deflection we get an outgoing direction of (sin(deflection), 0, cos(deflection))
	HELIODETIC_UNITVECTOR sun;
	sun.SetCoords(nxmath::sind(deflection), 0, nxmath::cosd(deflection));

	return sun & look;
}

bool SKTRAN_SolarTransmission_3D::TransmissionAtVector	(	const HELIODETIC_VECTOR&				point,
															double&									transmission ) const
{
	bool ok = true;

	double magnitude = point.Magnitude();
	double cossza = point.Z() * (magnitude>0.0?1.0/magnitude:0.0);
	double alt    = RayFactory()->CoordsPtr()->RadiusToAltitude( magnitude );
	double slon   = atan2( point.Y(), point.X() );

	return ok = ok && Transmission_Interpolate( cossza, alt, slon, transmission );
}

bool SKTRAN_SolarTransmission_3D::Transmission_Interpolate ( double cossza, double alt, double slon, double& transmission ) const
{
	bool ok = true;
	double				altweights[2];
	double				szaweights[2];
	double				slonweights[2];
	size_t				altindex[2];
	size_t				szaindex[2];
	size_t				slonindex[2];

	size_t				numalt, numsza, numslon;

	ok = ok && AltWeightsForProfile( alt, altweights, altindex, numalt );
	ok = ok && CosSzaWeights( cossza, szaweights, szaindex, numsza );
	ok = ok && SlonWeights( slon, slonweights, slonindex, numslon );

	transmission = 0;
	for( size_t szaidx = 0; szaidx < numsza; szaidx++ )
	{
		for( size_t slonidx = 0; slonidx < numslon; slonidx++ )
		{
			for( size_t altidx = 0; altidx < numalt; altidx++ )
			{
				double trans = m_transmission[ sub2ind( szaindex[szaidx], slonindex[slonidx], altindex[altidx] ) ];
				if( trans < 0.0 )
				{
					//std::mutex::scoped_lock lock; // This seems like it would be useful to prevent the same index from being filled by each 
					// thread, but it turns out that the locking does more harm than good. 
					ok = ok && FillTableAtIndex( szaindex[szaidx], slonindex[slonidx], altindex[altidx] );
					trans = m_transmission[ sub2ind( szaindex[szaidx], slonindex[slonidx], altindex[altidx] ) ];
				}
				transmission += std::log(trans) * altweights[altidx] * szaweights[szaidx] * slonweights[slonidx];
			}
		}
	}
	transmission = std::exp(transmission);

	return ok;
}

bool SKTRAN_SolarTransmission_3D::Deflection_Interpolate ( double cossza, double alt, double slon, double& deflection ) const
{
	bool ok = true;
	double				altweights[2];
	double				szaweights[2];
	double				slonweights[2];
	size_t				altindex[2];
	size_t				szaindex[2];
	size_t				slonindex[2];

	size_t				numalt, numsza, numslon;

	ok = ok && AltWeightsForProfile( alt, altweights, altindex, numalt );
	ok = ok && CosSzaWeights( cossza, szaweights, szaindex, numsza );
	ok = ok && SlonWeights( slon, slonweights, slonindex, numslon );

	deflection = 0;
	for( size_t szaidx = 0; szaidx < numsza; szaidx++ )
	{
		for( size_t slonidx = 0; slonidx < numslon; slonidx++ )
		{
			for( size_t altidx = 0; altidx < numalt; altidx++ )
			{
				double defl = m_deflectionangle[ sub2ind( szaindex[szaidx], slonindex[slonidx], altindex[altidx] ) ];
				deflection += deflection * altweights[altidx] * szaweights[szaidx] * slonweights[slonidx];
			}
		}
	}

	return ok;
}

bool SKTRAN_SolarTransmission_3D::FillTable_ClassSpecific ()
{
	m_transmission.assign( m_transmission.size(), -1.0 );

	if( m_prefilltable )
	{
		return PrefillTable();
	}
	else
	{
		return true;
	}
}

bool SKTRAN_SolarTransmission_3D::PrefillTable()
{
	bool ok = true;

	if(m_prefillgrid.NumAngles() == 0)
	{
		m_prefillgrid = CosSZAGrid();
	}

	HELIODETIC_UNITVECTOR locunit;
	HELIODETIC_VECTOR     loc;
	HELIODETIC_POINT      quadlocation;

	double sinsza, cossza, slon;

	std::vector<double> altitudevector;
	std::vector<double> cosszavector;
	std::vector<double> opticaldepthvector;
	std::vector<double> deflectionvector;

	// Treat every solar longitude separately
	for (int slonidx = 0; slonidx < m_slongrid.NumAngles(); slonidx++)
	{
		altitudevector.clear();
		cosszavector.clear();
		opticaldepthvector.clear();
		deflectionvector.clear();

		slon = m_slongrid.SLON( slonidx );

		std::vector<std::unique_ptr<SKTRAN_RayOptical_Base>> rays;
		rays.resize(m_prefillgrid.NumAngles());

		#pragma omp parallel for schedule(dynamic,1) private(locunit, loc, cossza, sinsza)
		for (int szaidx = 0; szaidx < m_prefillgrid.NumAngles(); szaidx++)
		{
			cossza = m_prefillgrid.CosineSZA( szaidx );
			sinsza = sqrt(1 - cossza*cossza);

			if(cossza < 0)
			{
				// Angles greater than 90 degrees we don't have to do
				continue;
			}

			locunit.SetCoords( sinsza * cos( slon ), sinsza * sin( slon ), cossza );
			// Slightly outside the atmosphere
			loc.SetCoords( locunit, CoordinatesPtr()->AltitudeToRadius( CoordinatesPtr()->TOAAltitude() + 1.0));

			rays[szaidx] = CreateRayAndCalcFullTransmission(loc);

		}
		for (int rayidx = 0; rayidx < rays.size(); rayidx++)
		{
			if(!rays[rayidx])
			{
				continue;
			}

			for(int quadidx = 0; quadidx < rays[rayidx]->GetNumQuadraturePoints(); quadidx++)
			{
				rays[rayidx]->Storage()->LocationOfPoint(quadidx, &quadlocation);

				auto raydirection = rays[rayidx]->Storage()->AverageLookVectorTowardsObserver(quadidx);
				nxVector dirvector(raydirection.X(), raydirection.Y(), raydirection.Z());

				cosszavector.push_back(quadlocation.CosSZA());
				altitudevector.push_back(max(quadlocation.Altitude(), CoordinatesPtr()->GroundAltitude()));
				opticaldepthvector.push_back(rays[rayidx]->OpticalDepthArray()[quadidx]);
				deflectionvector.push_back(dirvector.AngleTo(nxVector(0, 0, 1)));
			}
		}
		ok = ok && InterpolateAndFillSLON(slonidx, altitudevector, cosszavector, opticaldepthvector, deflectionvector);
	}

	// DumpTable();
	return ok;
}

bool SKTRAN_SolarTransmission_3D::InterpolateAndFillSLON( int slonidx, std::vector<double>& altitude, std::vector<double>& cossza, std::vector<double>& opticaldepth, std::vector<double>& deflection)
{
	#pragma omp parallel for schedule(dynamic,1)
	for (int shellaltidx = 0; shellaltidx < HeightGrid().NumShells(); shellaltidx++ )
	{
		std::vector<double> sortedcossza;
		std::vector<double> sortedod;
		std::vector<double> sorteddeflect;
		std::vector<std::tuple<double, double, double>> opticaldepthvector;

		double shellheight = HeightGrid().ShellHeight()[shellaltidx];

		for (int altidx = 0; altidx < altitude.size(); altidx++)
		{
			if (std::abs(altitude[altidx] - shellheight) < 0.1)
			{
				opticaldepthvector.push_back(std::tuple<double, double, double>(cossza[altidx], opticaldepth[altidx], deflection[altidx]));
			}
		}


		std::sort(std::begin(opticaldepthvector), std::end(opticaldepthvector), [](std::tuple<double, double, double> a, std::tuple<double, double, double> b){return std::get<0>(a) < std::get<0>(b);});
		
		sortedcossza.reserve(opticaldepthvector.size());
		sortedod.reserve(opticaldepthvector.size());


		for (int idx = 0; idx < opticaldepthvector.size(); idx++)
		{
			sortedcossza.push_back(std::get<0>(opticaldepthvector[idx]));
			sortedod.push_back(std::get<1>(opticaldepthvector[idx]));
			sorteddeflect.push_back(std::get<2>(opticaldepthvector[idx]));
		}

		for (int cosszaidx = 0; cosszaidx < CosSZAGrid().size(); cosszaidx++)
		{
			int transidx = sub2ind(cosszaidx, slonidx, shellaltidx);
			m_transmission[transidx] = nxLinearInterpolate::EvaluateYatX(CosSZAGrid().At(cosszaidx), sortedcossza, sortedod, nxLinearInterpolate::ENUM_MISSINGVALUE, 1000.0);
			m_deflectionangle[transidx] = nxLinearInterpolate::EvaluateYatX(CosSZAGrid().At(cosszaidx), sortedcossza, sorteddeflect, nxLinearInterpolate::ENUM_MISSINGVALUE, 0.0);
			m_transmission[transidx] = exp(-1.0 * m_transmission[transidx]);
		}
	}

	return true;
}


bool SKTRAN_SolarTransmission_3D::InitializeGeometry(	std::shared_ptr<const SKTRAN_CoordinateTransform_V2>&		coords,
														const SKTRAN_GridDefRayTracingShells_V21& altgrid,
														const SKTRAN_GridDefCosSZA_V21& cosszagrid,
														const SKTRAN_GridDefSLON_V21& slongrid )
{
	bool ok = true;

	ok = SKTRAN_SolarTransmission_2D::SetGeometry( coords, altgrid, cosszagrid );
	m_slongrid.DeepCopy( slongrid );
	m_numcosangles = CosSZAGrid().NumAngles();
	m_numslon      = m_slongrid.NumAngles();
	m_numalts      = HeightGrid().NumShells();
	m_transmission.resize( m_numcosangles * m_numslon * m_numalts );
	m_deflectionangle.resize( m_numcosangles * m_numslon * m_numalts );

	return ok;
}

void SKTRAN_SolarTransmission_3D::DumpTable()
{
	std::ofstream file;
	file.open("/home/dannyz/tmp/solartable.txt");

	int c = 0;
	for (int aidx = 0; aidx < CosSZAGrid().NumAngles(); aidx++)
	{
		for(int hidx = 0; hidx < HeightGrid().NumCells(); hidx++)
		{
			int transidx = sub2ind(aidx, 0, hidx);
			file << HeightGrid().At(hidx) << " " << CosSZAGrid().CosineSZA(aidx) << " " << m_transmission[transidx] << " " << m_deflectionangle[transidx] << std::endl;
		}
	}

}