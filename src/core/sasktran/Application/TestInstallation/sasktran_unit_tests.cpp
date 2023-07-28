
#ifdef SKTRAN_CATCH2_VERSION3
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_all.hpp>

#else
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#endif

//#define _CRTDBG_MAP_ALLOC // uncomment for memory dump
#include <stdlib.h>
#if defined (NX_WINDOWS)
	#include <crtdbg.h>
#endif
#include <sasktran.h>
#include "include/short_test_suite.h"
#include "include/test_mcamf.h"

extern int testDelaunaySphere();
extern void TestRay();
extern void TestSingleScatter();
extern void TestLatLonUnitSphere();
extern bool test_hitran();
extern bool test_so2();
extern bool test_o2_o2();
extern bool test_linearcombo();
extern bool test_msisatomicO();

double g_errortolerance = 1.0E-4;


TEST_CASE( "Standard HR Test ", "[HR]" )
{
	SKTRAN_Short_Test_HR	hr(g_errortolerance);
    REQUIRE( hr.RunStandardTest() );
}

TEST_CASE( "Standard MC Test ", "[MC]" )
{
	SKTRAN_Short_Test_MC	mc(g_errortolerance);
    REQUIRE( mc.RunStandardTest() );
}

TEST_CASE("Inelastic MC Test", "[MC][INELASTIC]")
{
	SKTRAN_Short_Test_Inelastic_MC mc(g_errortolerance, true, true, true);
	REQUIRE(mc.RunStandardTest());
}


TEST_CASE( "SO2", "[SO2]")
{
	REQUIRE( test_so2() );
}

TEST_CASE( "O2-O2 CIA", "[O2-O2]")
{
	REQUIRE(test_o2_o2());
}

TEST_CASE( "LINEAR COMBO", "[LINEARCOMBO]")
{
	REQUIRE( test_linearcombo() );
}

TEST_CASE( "MSIS Atomic Oxygen", "[MSIS O]")
{
	REQUIRE( test_msisatomicO() );
}

TEST_CASE(" MC AMF Test", "[MC][AMF]")
{
	SKTRAN_MCAMF_Test mcamf(false, g_errortolerance);
	REQUIRE(mcamf.RunAllTests());
}





#include "include/simple_ray_test.h"
#include "include/compare_sktran_to_hr.h"
#include "include/test_hr.h"
#include "include/test_mc.h"
#include "include/test_do.h"
#include "include/test_mcamf.h"
#include <boost/timer/timer.hpp>
#include "include/testscattermatrices.h"
#include "include/test_curvedrays.h"


int	g_testdescription = 0;
bool g_doso           = true;
bool g_domc           = true;
bool g_dohr           = true;

/*-----------------------------------------------------------------------------
 *					GroundBaseTest		2014-4-8*/
/** **/
/*---------------------------------------------------------------------------*/

void GroundBasedTest()
{
	SKTRAN_Short_Test	tests;

	nxLog::Record(NXLOG_INFO, "Ground Based Test. This will test a ground based geometry with all 3 engines");
	tests.RunGroundBasedTests(true, false, false);
	nxLog::Record(NXLOG_INFO, "Ground Based Test  Completed.");
}


/*---------------------------------------------------------------------------
 *						CheckOptions
 *-------------------------------------------------------------------------*/

//
//bool CheckOptions( int argc, char** argv)
//{
//
//	nxGetOpt		  opt;
//	nxGetOptEntry*	  entry;
//	bool		      ok;
//	const char*		  s;
//	char			  c;
//	size_t			  i;
//
//	static const char*	  shortopts = "t:e:hHMS";
//	static struct option  long_options[] =
//	{
//		{"test",             1, 0, 't'},
//		{"help",             0, 0, 'h'},
//		{"HR",               0, 0, 'H'},
//		{"MC",               0, 0, 'M'},
//		{"SO",               0, 0, 'S'},
//		{"error",            1, 0, 'e'},
//		{0,                  0, 0, 0  }
//    };
//
//	opt.SetPrintErrors(true);
//	ok = opt.Parse( argc, argv, shortopts, long_options);
//	if (ok)
//	{
//		for (i=0; i < opt.NumOptions(); i++)
//		{
//			if (opt.GetNextOption(&entry))
//			{
//
//				s = entry->GetKeyStr();
//				c = s[0];
//				switch (c)
//				{
//				case 't' :  g_testdescription = atoi(entry->GetValue());
//							break;
//				case 'e' :  g_errortolerance  = atof(entry->GetValue());
//							break;
//
//				case 'H' :  g_dohr = false; break;
//				case 'M' :  g_domc = false; break;
//				case 'S' :  g_doso = false; break;
//
//				};
//			}
//		}
//		ok = (opt.NumArguments() == 0);
//		if (!ok)
//		{
//			nxLog::Record(NXLOG_WARNING,"CheckOptions, you cannot have any command line arguments if you specify any options");
//		}
//	}
////	if (!ok) Usage(1);
//	return ok;
//}
//

/*-----------------------------------------------------------------------------
 *					TestAdapative		 2014- 12- 10*/
/** **/
/*---------------------------------------------------------------------------*/

void TestHR()
{
	 SKTRAN_HR_Test hrt(false, 1e-5);
	 hrt.RunAllTests();
}


/*-----------------------------------------------------------------------------
 *					TestMC		 2014- 12- 11*/
/** **/
/*---------------------------------------------------------------------------*/

void TestMC()
{
	SKTRAN_MC_Test mct( false, 1e-7 );
	mct.RunAllTests();
}


/*-----------------------------------------------------------------------------
 *					TestPolarization		 2016- 12- 5*/
/** **/
/*---------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------
 *					TestMC		 2014- 12- 11*/
 /** **/
 /*---------------------------------------------------------------------------*/

void TestMCAMF()
{
	SKTRAN_MCAMF_Test mct(false, 1e-7);
	mct.RunAllTests();
}

/*-------------------------------------------------------------------------------+
 *					themain		2010-12-16*/
/** **/
/*---------------------------------------------------------------------------*/

//int themain ( int argc, char* argv[])
//{
////	std::timer::auto_cpu_timer	t;
//	nxLogConsole					thelog;
//	bool							ok;
//
//	ok   = CheckOptions(argc, argv);
//	if (ok)
//	{
//		switch (g_testdescription)
//		{
////		case  0: StandardTest();                             break;
//		case  1: GroundBasedTest();                          break;
//		case  2: testDelaunaySphere();                       break;
//		case  3: TestRay();                                  break;
//		case  4: TestSingleScatter();                        break;
//		case  5: TestLatLonUnitSphere();                     break;
//		case  6: TestHR();                                   break;
//		case  7: TestMC();                                   break;
//		#if defined (NX_WINDOWS)
//		case  8: Run_SKTRAN_Emission_Test();                 break;
//        case  9: TestPolarization();                         break;
//        case 10: TestEmission();                             break;
//		case 11: TestCurvedRays();							 break;
//		#endif
//		case 12: TestMCAMF();                                break;
//		case 13: test_hitran();								 break;
//		case 14: test_so2();								 break;
//		case 15: test_o2_o2();								 break;
//		case 16: test_linearcombo();						 break;
//		case 17: test_msisatomicO();						 break;
//
//		  //		case 9: TestDO(argc, argv);							break;
////		case 999: TestHR(); StandardTest();                 break;
//		// case 2: TestSingleScatter();						break;
//		// case 3: TestMonteCarloEngineClass();				break;
//		// case 4: TestLatLonUnitSphere();					break;
//		// case 5: TestSimpleInversionFramework(0.001);		break;
//		};
//	}
//	return 0;
//}


//
///*-----------------------------------------------------------------------------
// *					main		2014-4-8*/
///** **/
///*---------------------------------------------------------------------------*/
//
//int old_main(int argc, char* argv[])
//{
//    //_CrtSetBreakAlloc( 13589 );
//
//    {
//	nxLogConsole	thelog;
//
//	thelog.SetVerbose();
//	boost::timer::auto_cpu_timer t;
//
//	themain(argc, argv);
//
//    }
//
//#if defined (NX_WINDOWS)
//	_CrtDumpMemoryLeaks();
//#endif
//	return 0;
//}
//
//
