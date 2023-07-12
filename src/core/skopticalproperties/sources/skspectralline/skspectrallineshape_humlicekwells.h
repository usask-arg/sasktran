/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_HumlicekWells		2013-3-7*/
/** \ingroup spectralline
 *	[DEPRECATED}. This class is no longer supported and may not work properly.
 *	Use class skSpectralLineShape_VoigtKuntz.
 *
 *	This class is a wrapper for the skSpectralLineShape interface and does
 *	not contain any info about the actual Voigt Humlik algorithm. The work
 *	class is stored inside the storage buffer, see CreateStorageBuffer, and
 *	can be kept out of the interface include files.
 *
 *	The class implements a quick and accurate Voigt profile based on
 *	the HUMLIK.FOR code derived from the paper: 
 *	"Rapid Approximation to the Voigt/Faddeeva function and its derivatives",
 *	R. J Wells, JQSRT 62, 1999, 29-48
 *
 **/
/*---------------------------------------------------------------------------*/

class skSpectralLineShape_HumlicekWells : public skSpectralLineShape
{
	public:
							skSpectralLineShape_HumlicekWells		( );
		virtual			   ~skSpectralLineShape_HumlicekWells		( );
		virtual bool		LineShapeFunction						( double nu, double* uservalue, const skSpectralLine* spectralline, double maxlinestrength, skSpectralLineShapeStorageBuffer* storagebuffer ) const;
		virtual bool		ConfigureLineParameters					( const skSpectralLine* spectralline,  double	temperature, double	pressure, const GEODETIC_INSTANT& geopt, skClimatology* atmopshericstate, skSpectralLineShapeStorageBuffer*	storagebuffer );
		virtual bool		CreateStorageBuffer						( skSpectralLineShapeStorageBuffer** storagebuffer );
};

