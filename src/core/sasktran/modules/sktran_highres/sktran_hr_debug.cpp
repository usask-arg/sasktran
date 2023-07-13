#include "include/sktran_hr_internals.h"


void DumpRay( SKTRAN_RayOptical_Base* baseray, const SKTRAN_TableOpticalProperties_Base& opttable, const SKTRAN_Source_Term& source )
{
	bool ok = true;
	SKTRAN_RayOptical_Straight*	ray = reinterpret_cast<SKTRAN_RayOptical_Straight*>(baseray);
	HELIODETIC_UNITVECTOR		look;
	double						ds;
	double						s = 0;
    HELIODETIC_POINT            obspt;
	HELIODETIC_POINT			centerpoint;
	HELIODETIC_POINT			startpoint;
	double						sourceval;
	double						sourcemid;
	double k;
	double opticaldepth;
	double fullsource;

	std::ofstream				outfile;
	std::stringstream			sstream;
	sstream << SKTRAN_HR_DEBUG_DIR << "ray_source.txt";
	outfile.open(sstream.str().c_str());
	ok = ok && ray->Coordinates()->HelioVectorToHelioPoint(ray->GetObserver(), &obspt );
    if(ok){
        SKTRAN_SourceTermQueryObject_Simple qStart (obspt, ray->LookVector());
        SKTRAN_SourceTermQueryObject_Simple qCenter(obspt, ray->LookVector());
	    for( size_t cellidx = 0; cellidx < ray->GetNumCells(); cellidx++ )
	    {
		    look  = ray->Storage()->AverageLookVectorAwayFromObserver( cellidx);
		    ds    = ray->Storage()->CellLength( cellidx);
		    s     = ray->Storage()->DistanceOfPointFromOrigin( cellidx);

		    ok    = ok && ray->Storage()->LocationOfPoint( cellidx, &startpoint );
		    ok    = ok && ray->Storage()->CellMidPoint   ( cellidx, &centerpoint );
            qStart .UpdateQuery( startpoint,  look );
            qCenter.UpdateQuery( centerpoint, look );
		    ok    = ok && source.SourceTermAtPoint( qStart,  &sourceval );
		    ok    = ok && source.SourceTermAtPoint( qCenter, &sourcemid  );

		    //s += ds/2;

		    k = opttable.TotalExtinctionPerCM( centerpoint ) * 100;
		    opticaldepth =  ray->OpticalDepthArray().at( cellidx);
		    fullsource   = (sourcemid * ( 1 - exp(-1.0*k*ds)) / k) * exp(-1.0*opticaldepth);
		    outfile.precision(50);
		    outfile.setf(std::ios::fixed);
		    outfile.setf(std::ios::showpoint);
		    outfile << s << " " << ds << " " << sourceval << " " << opticaldepth << " " << fullsource << std::endl;
	    }

	    outfile.close();
	} else{
        nxLog::Record( NXLOG_WARNING, "DumpRay, Couldn't dump ray.");
	}
}
