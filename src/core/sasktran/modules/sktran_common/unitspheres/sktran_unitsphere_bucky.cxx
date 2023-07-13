#include "../sktran_common.h"

/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereBucky_V2::SKTRAN_UnitSphereBucky_V2		2007-12-17*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_UnitSphereBucky_V2::SKTRAN_UnitSphereBucky_V2()
{
	bool	ok;
	size_t	unitCtr;

	ok = AllocateVertices( 92 );
	m_phi  = (sqrt(5.0)+1)/ 2.0;
	if (ok)
	{

		UnitVectorAtVar( 0).SetCoords(	0,	 1,	 3*m_phi);
		UnitVectorAtVar( 1).SetCoords(	0,	 1,	-3*m_phi);
		UnitVectorAtVar( 2).SetCoords(	0,	-1,	 3*m_phi);
		UnitVectorAtVar( 3).SetCoords(	0,	-1,	-3*m_phi);

		UnitVectorAtVar( 4).SetCoords(	 1,	 3*m_phi,	0);
		UnitVectorAtVar( 5).SetCoords(	 1,	-3*m_phi,	0);
		UnitVectorAtVar( 6).SetCoords(	-1,	 3*m_phi,	0);
		UnitVectorAtVar( 7).SetCoords(	-1,	-3*m_phi,	0);

		UnitVectorAtVar( 8).SetCoords(	 3*m_phi,	0,	 1);
		UnitVectorAtVar( 9).SetCoords(	 3*m_phi,	0,	-1);
		UnitVectorAtVar(10).SetCoords(	-3*m_phi,	0,	 1);
		UnitVectorAtVar(11).SetCoords(	-3*m_phi,	0,	-1);

		UnitVectorAtVar(12).SetCoords(	 2,	 (1+2*m_phi),	 m_phi);
		UnitVectorAtVar(13).SetCoords(	 2,	 (1+2*m_phi),	-m_phi);
		UnitVectorAtVar(14).SetCoords(	 2,	-(1+2*m_phi),	 m_phi);
		UnitVectorAtVar(15).SetCoords(	 2,	-(1+2*m_phi),	-m_phi);
		UnitVectorAtVar(16).SetCoords(	-2,	 (1+2*m_phi),	 m_phi);
		UnitVectorAtVar(17).SetCoords(	-2,	 (1+2*m_phi),	-m_phi);
		UnitVectorAtVar(18).SetCoords(	-2,	-(1+2*m_phi),	 m_phi);
		UnitVectorAtVar(19).SetCoords(	-2,	-(1+2*m_phi),	-m_phi);

		UnitVectorAtVar(20).SetCoords(	 (1+2*m_phi),	 m_phi,	 2);
		UnitVectorAtVar(21).SetCoords(	 (1+2*m_phi),	 m_phi,	-2);
		UnitVectorAtVar(22).SetCoords(	 (1+2*m_phi),	-m_phi,	 2);
		UnitVectorAtVar(23).SetCoords(	 (1+2*m_phi),	-m_phi,	-2);
		UnitVectorAtVar(24).SetCoords(	-(1+2*m_phi),	 m_phi,	 2);
		UnitVectorAtVar(25).SetCoords(	-(1+2*m_phi),	 m_phi,	-2);
		UnitVectorAtVar(26).SetCoords(	-(1+2*m_phi),	-m_phi,	 2);
		UnitVectorAtVar(27).SetCoords(	-(1+2*m_phi),	-m_phi,	-2);

		UnitVectorAtVar(28).SetCoords(	 m_phi,	 2,	(1+2*m_phi));
		UnitVectorAtVar(29).SetCoords(	 m_phi,	 2,	-(1+2*m_phi));
		UnitVectorAtVar(30).SetCoords(	 m_phi,	-2,	 (1+2*m_phi));
		UnitVectorAtVar(31).SetCoords(	 m_phi,	-2,	-(1+2*m_phi));
		UnitVectorAtVar(32).SetCoords(	-m_phi,	 2,	 (1+2*m_phi));
		UnitVectorAtVar(33).SetCoords(	-m_phi,	 2,	-(1+2*m_phi));
		UnitVectorAtVar(34).SetCoords(	-m_phi,	-2,	 (1+2*m_phi));
		UnitVectorAtVar(35).SetCoords(	-m_phi,	-2,	-(1+2*m_phi));

		UnitVectorAtVar(36).SetCoords(	 1,	 (2+m_phi),	 2*m_phi);
		UnitVectorAtVar(37).SetCoords(	 1,	 (2+m_phi),	-2*m_phi);
		UnitVectorAtVar(38).SetCoords(	 1,	-(2+m_phi),	 2*m_phi);
		UnitVectorAtVar(39).SetCoords(	 1,	-(2+m_phi),	-2*m_phi);
		UnitVectorAtVar(40).SetCoords(	-1,	 (2+m_phi),	 2*m_phi);
		UnitVectorAtVar(41).SetCoords(	-1,	 (2+m_phi),	-2*m_phi);
		UnitVectorAtVar(42).SetCoords(	-1,	-(2+m_phi),	 2*m_phi);
		UnitVectorAtVar(43).SetCoords(	-1,	-(2+m_phi),	-2*m_phi);

		UnitVectorAtVar(44).SetCoords(	 (2+m_phi),	 2*m_phi,	 1);
		UnitVectorAtVar(45).SetCoords(	 (2+m_phi),	 2*m_phi,	-1);
		UnitVectorAtVar(46).SetCoords(	 (2+m_phi),	-2*m_phi,	 1);
		UnitVectorAtVar(47).SetCoords(	 (2+m_phi),	-2*m_phi,	-1);
		UnitVectorAtVar(48).SetCoords(	-(2+m_phi),	 2*m_phi,	 1);
		UnitVectorAtVar(49).SetCoords(	-(2+m_phi),	 2*m_phi,	-1);
		UnitVectorAtVar(50).SetCoords(	-(2+m_phi),	-2*m_phi,	 1);
		UnitVectorAtVar(51).SetCoords(	-(2+m_phi),	-2*m_phi,	-1);

		UnitVectorAtVar(52).SetCoords(	2*m_phi,	 1,	 (2+m_phi));
		UnitVectorAtVar(53).SetCoords(	2*m_phi,	 1,	-(2+m_phi));
		UnitVectorAtVar(54).SetCoords(	2*m_phi,	-1,	 (2+m_phi));
		UnitVectorAtVar(55).SetCoords(	2*m_phi,	-1,	-(2+m_phi));
		UnitVectorAtVar(56).SetCoords(	-2*m_phi,	 1,	 (2+m_phi));
		UnitVectorAtVar(57).SetCoords(	-2*m_phi,	 1,	-(2+m_phi));
		UnitVectorAtVar(58).SetCoords(	-2*m_phi,	-1,	 (2+m_phi));
		UnitVectorAtVar(59).SetCoords(	-2*m_phi,	-1,	-(2+m_phi));
		// Unitize the first 60 points.
		for (unitCtr = 0; unitCtr < 60; ++unitCtr ) UnitVectorAtVar(unitCtr) *= 1.0/sqrt(9*m_phi+10);

		// 12 points on the pentagon centers.
		UnitVectorAtVar(60).SetCoords(	0,	 1,	 m_phi);
		UnitVectorAtVar(61).SetCoords(	0,	 1,	-m_phi);
		UnitVectorAtVar(62).SetCoords(	0,	-1,	 m_phi);
		UnitVectorAtVar(63).SetCoords(	0,	-1,	-m_phi);

		UnitVectorAtVar(64).SetCoords(	 1,	 m_phi,	0);
		UnitVectorAtVar(65).SetCoords(	 1,	-m_phi,	0);
		UnitVectorAtVar(66).SetCoords(	-1,	 m_phi,	0);
		UnitVectorAtVar(67).SetCoords(	-1,	-m_phi,	0);

		UnitVectorAtVar(68).SetCoords(	 m_phi,	0,	 1);
		UnitVectorAtVar(69).SetCoords(	 m_phi,	0,	-1);
		UnitVectorAtVar(70).SetCoords(	-m_phi,	0,	 1);
		UnitVectorAtVar(71).SetCoords(	-m_phi,	0,	-1);

		// Unitize these 12 points.
		for( unitCtr = 60; unitCtr < 72; ++unitCtr ) UnitVectorAtVar(unitCtr) *= 1.0/sqrt(m_phi+2);

		// 20 points on the hexagon centers.
		UnitVectorAtVar(72).SetCoords(	 1,	 1,	 1);
		UnitVectorAtVar(73).SetCoords(	 1,	 1,	-1);
		UnitVectorAtVar(74).SetCoords(	 1,	-1,	 1);
		UnitVectorAtVar(75).SetCoords(	 1,	-1,	-1);
		UnitVectorAtVar(76).SetCoords(	-1,	 1,	 1);
		UnitVectorAtVar(77).SetCoords(	-1,	 1,	-1);
		UnitVectorAtVar(78).SetCoords(	-1,	-1,	 1);
		UnitVectorAtVar(79).SetCoords(	-1,	-1,	-1);

		UnitVectorAtVar(80).SetCoords(	0,	 m_phi,	 1/m_phi);
		UnitVectorAtVar(81).SetCoords(	0,	 m_phi,	-1/m_phi);
		UnitVectorAtVar(82).SetCoords(	0,	-m_phi,	 1/m_phi);
		UnitVectorAtVar(83).SetCoords(	0,	-m_phi,	-1/m_phi);

		UnitVectorAtVar(84).SetCoords(	 m_phi,	 1/m_phi,	0);
		UnitVectorAtVar(85).SetCoords(	 m_phi,	-1/m_phi,	0);
		UnitVectorAtVar(86).SetCoords(	-m_phi,	 1/m_phi,	0);
		UnitVectorAtVar(87).SetCoords(	-m_phi,	-1/m_phi,	0);

		UnitVectorAtVar(88).SetCoords(	 1/m_phi,	0,	 m_phi);
		UnitVectorAtVar(89).SetCoords(	 1/m_phi,	0,	-m_phi);
		UnitVectorAtVar(90).SetCoords(	-1/m_phi,	0,	 m_phi);
		UnitVectorAtVar(91).SetCoords(	-1/m_phi,	0,	-m_phi);
		// Unitize these 20 points.
		for( unitCtr = 72; unitCtr < 92; ++unitCtr ) UnitVectorAtVar(unitCtr) *= 1.0/sqrt(3.0);
	}
	InitializeLookupTable();
	NXTRACE_ONCEONLY(firsttime,("***** TODO ***** SKTRAN_UnitSphereBucky_V2::SKTRAN_UnitSphereBucky_V2, Must implement a cubature weights algorithm before using this code\n"));
	NXASSERT((false));
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_UnitSphereBucky_V2, Error creating Bucky Sphere vertices");
	}
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereBucky_V2::~SKTRAN_UnitSphereBucky_V2		2008-1-4*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_UnitSphereBucky_V2::~SKTRAN_UnitSphereBucky_V2()
{
}
