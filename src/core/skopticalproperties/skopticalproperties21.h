#pragma once
#if !defined(SKOPTICALPROPERTIES_H_INCLUDED)
#define SKOPTICALPROPERTIES_H_INCLUDED

#define SKOPTICALPROPERTIES_VERSION 220

#include <math.h>
#include <nxbase_core.h>
#include <nxbase_math.h>
#include <nxbase_threads.h>
#include <limits>
#include <sasktranif.h>
#include <skclimatology21.h>

#if !defined(SKCLIMATOLOGY_VERSION) || (SKCLIMATOLOGY_VERSION < 210)
#error The skclimatology version does not meet the requirements of skopticalproperties. Please upgrade your version of skclimatology from SVN.
#endif

#include "sources/emissions/skemission.h"
#include "sources/skrefractiveindex.h"
#include "sources/brdf/sktran_brdf.h"
#include "sources/brdf/sktran_albedo.h"
#include "sources/skparticledist.h"
#include "sources/skphasematrix.h"
#include "sources/skopticalpropertylistentry.h"
#include "sources/skmiesphericalparticle.h"
#include "sources/sknonsphericalparticle.h"
#include "sources/skrtscatter_phasematrix.h"
#include "sources/skabsorptiontable.h"
#include "sources/skabsorptiontablepressure.h"
#include "sources/skrtextinction_tabulatedphasematrix.h"
#include "sources/baum_bulk_icecloudmodel/skrtscatter_tabulatedcirrusproperties.h"
#include "sources/baum_bulk_icecloudmodel/baumicecrystals_database2014.h"
#include "sources/emissions/sktranatmosphericemission.h"
#include "sources/sktranatmosphericstate.h"
#include "sources/skconvolvedopticalproperties.h"
#include "sources/skrtextinctionmarthybridprofile.h"
#include "sources/skspectralline/skspectralline.h"
#include "sources/solarspectrum/sksolarspectrum.h"
#include "sources/emissions/skemissiontabulatedheightwavelength.h"
#include "sources/emissions/skemissionthermal.h"
#include "sources/emissions/hitranline_upperstates.h"
#include "sources/o2-o2/o4_cross_sections.h"
#include "sources/skuserscatterconstantheight.h"


#endif

