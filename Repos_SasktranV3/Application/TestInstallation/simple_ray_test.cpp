#include <sasktran.h>
#include "include/simple_ray_test.h"
#include <memory>


static void MakeRayTracingGrid( SKTRAN_GridDefRayTracingShells_V21* shellgrid )
{
	double shellheights[1000];
	
	for(int i = 0; i < 1000; i++)
	{
		shellheights[i] = 100+i*100;
	}
	shellgrid->ConfigureHeights( shellheights, 1000 );
}

static void MakeRayTracer( SKTRAN_RayTracer_Shells& raytracer )
{
	std::shared_ptr<SKTRAN_GridDefRayTracingShells_V21> shellgrid ( new SKTRAN_GridDefRayTracingShells_V21);


	MakeRayTracingGrid( shellgrid.get() );
	raytracer.Initialize( shellgrid );
}

void TestRay()
{
	const double Rearth = 6372.0*1000;
	nxVector look(-1,0,0);
	nxVector obs(Rearth+100*1000,0,Rearth+40*1000);
	nxVector sun(0,0,1);

	HELIODETIC_VECTOR obs_v;
	obs_v.SetCoords(obs.X(), obs.Y(), obs.Z());
	HELIODETIC_UNITVECTOR look_v;
	look_v.SetCoords(look.X(),look.Y(),look.Z());

	SKTRAN_CoordinateTransform_V2 coords;
	coords.ConfigureCoordinates(obs,look,53000.0,sun);

//	SKTRAN_RayTracer_Shells raytracer;
//	MakeRayTracer( raytracer );

//	SKTRAN_RayGeometry_Straight	raygeo;
//	SKTRAN_RayOptical_Straight	rayoptical;

//	raygeo.Initialize( coords, obs_v, look_v);
//	rayoptical.SetGeometryRay( raygeo );

//	raytracer.CreateRay( rayoptical );

	
}
