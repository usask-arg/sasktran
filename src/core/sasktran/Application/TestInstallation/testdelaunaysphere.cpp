#include <sasktran.h>
#include "include/testDelaunaySphere.h"

static void testDelaunaySphere_loadVerts(nxVector* verts)
{
verts[0].SetCoords(7.211550e-001,-6.700301e-001,1.760544e-001);
verts[1].SetCoords(2.783061e-001,-7.427018e-002,9.576167e-001);
verts[2].SetCoords(-6.800764e-001,2.768693e-001,6.788516e-001);
verts[3].SetCoords(-1.452466e-001,7.232536e-001,6.751353e-001);
verts[4].SetCoords(-6.918556e-001,6.031491e-001,3.969218e-001);
verts[5].SetCoords(-8.615858e-001,1.583484e-001,4.822818e-001);
verts[6].SetCoords(-4.248235e-001,8.884607e-001,1.736736e-001);
verts[7].SetCoords(4.519162e-001,8.429297e-001,2.919610e-001);
verts[8].SetCoords(7.739878e-001,-1.911579e-001,6.036568e-001);
verts[9].SetCoords(1.469153e-001,4.265497e-001,8.924524e-001);
verts[10].SetCoords(2.576160e-001,2.094600e-001,9.432712e-001);
verts[11].SetCoords(-3.418158e-001,2.944951e-001,8.924318e-001);
verts[12].SetCoords(3.974158e-001,-6.494992e-001,6.482372e-001);
verts[13].SetCoords(1.329560e-001,-4.437941e-001,8.862108e-001);
verts[14].SetCoords(4.425646e-001,-5.217816e-001,7.293014e-001);
verts[15].SetCoords(-3.883911e-001,-9.124181e-001,1.290179e-001);
verts[16].SetCoords(2.235720e-001,4.587386e-001,8.599852e-001);
verts[17].SetCoords(-2.132260e-001,-7.284580e-001,6.510635e-001);
verts[18].SetCoords(-1.466553e-001,9.317043e-001,3.322941e-001);
verts[19].SetCoords(-4.349674e-001,-3.412865e-001,8.332628e-001);
verts[20].SetCoords(-9.055946e-001,3.275963e-001,2.694051e-001);
verts[21].SetCoords(7.489971e-001,-1.956700e-001,6.330217e-001);
verts[22].SetCoords(4.355426e-001,6.280441e-001,6.448746e-001);
verts[23].SetCoords(6.884511e-001,-4.366748e-003,7.252697e-001);
verts[24].SetCoords(-6.134717e-001,4.045163e-002,7.886800e-001);
verts[25].SetCoords(7.689809e-001,-1.660477e-001,6.173301e-001);
verts[26].SetCoords(3.356863e-001,-9.246494e-001,1.798284e-001);
verts[27].SetCoords(-6.645668e-002,-3.814261e-001,9.220074e-001);
verts[28].SetCoords(-7.118514e-002,5.364934e-001,8.408969e-001);
verts[29].SetCoords(-1.935484e-001,4.314456e-001,8.811321e-001);
verts[30].SetCoords(8.800141e-002,9.576796e-001,2.740540e-001);
verts[31].SetCoords(2.712852e-001,5.977259e-002,9.606412e-001);
verts[32].SetCoords(2.180362e-001,9.131705e-001,3.443543e-001);
verts[33].SetCoords(-5.674453e-001,-1.761411e-001,8.043508e-001);
verts[34].SetCoords(5.750737e-001,-4.790714e-001,6.631597e-001);
verts[35].SetCoords(7.374524e-001,6.357482e-001,2.280092e-001);
verts[36].SetCoords(2.869505e-001,-9.546999e-001,7.878758e-002);
verts[37].SetCoords(-1.146453e-002,-3.189397e-001,9.477057e-001);
verts[38].SetCoords(3.845163e-001,9.183578e-001,9.362772e-002);
verts[39].SetCoords(1.465275e-003,-6.627679e-001,7.488234e-001);
verts[40].SetCoords(-6.590304e-001,7.773743e-002,7.480882e-001);
verts[41].SetCoords(-8.447049e-001,-5.039623e-001,1.802656e-001);
verts[42].SetCoords(6.466154e-001,-2.201983e-002,7.624983e-001);
verts[43].SetCoords(4.260224e-001,-7.333112e-001,5.298675e-001);
verts[44].SetCoords(8.816671e-001,2.962859e-001,3.672571e-001);
verts[45].SetCoords(7.056712e-001,5.401769e-001,4.585162e-001);
verts[46].SetCoords(5.100350e-001,8.550948e-001,9.315115e-002);
verts[47].SetCoords(8.496354e-001,-5.382789e-003,5.273430e-001);
verts[48].SetCoords(9.497799e-002,5.928494e-001,7.996929e-001);
verts[49].SetCoords(5.557151e-001,-5.585787e-001,6.157683e-001);
}

int testDelaunaySphere()
{	

    SKTRAN_UnitSphere_Delaunay_nonTabledLookup sphere;
	nxVector openAxis;
    const size_t numVecs = 50;
    nxVector unitVecs[numVecs];

    testDelaunaySphere_loadVerts(unitVecs);
    openAxis.SetCoords(0.0, 0.0, -1.0);
	sphere.CreateTriangulation(unitVecs, numVecs, &openAxis);

	nxVector test(-0.262749243333333,-0.121910366666667,0.891133433333333);
	size_t testindex[3];
	double testweight[3];

	size_t speedHelper = 30;
	sphere.Triangulate(test, testindex, testweight, 3);
	sphere.Triangulate(test, testindex, testweight, 3, speedHelper);

    return 0;
}

