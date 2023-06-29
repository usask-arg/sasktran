
/*-----------------------------------------------------------------------------
 *					SKTRANSO_RayStorageBaseGeometry		2007-11-16*/
/** @ingroup rays
 *	This is a class used to trace straight line rays through the atmosphere.
 *	It generates an array of elements along a ray trajectory corresponding
 *	to the atmospheric cells that the ray passes through. The properties of
 *	each ray element are defined by class SKTRAN_ElementBaseGeometry_V2. 
 *	This class does not allocate the storage for the cell elements but provides a
 *	purely abstract function, AllocatePathElements, which is implemented in
 *	derived classes. The derived class allocates storage and provides an
 *	implentation of RayElement to return a pointer to the specific cell.
 *
 *	\p Usage
 *	The derived class will normally call InitializeGeometry folloewed
 *	by a call to TraceRaysThroughShells.  The derived class will also provide
 *	implementations of AllocatePathElements and RayElement()
 *
 *	\p Testing
 *	This code was tested on 2007-11-16 against Chris Roths code for approximately
 *	150000 geometries both inside and above the atmosphere that exercised all 5
 *	geometries.  In all cases this code produced results exactly identical (all
 *	bits of double precision number identical) to Chris's code.
 *
 **/
/*---------------------------------------------------------------------------*/

class SKTRANSO_RayStorageBaseGeometry : public SKTRAN_RayStorage_Straight
{		


	public:
													SKTRANSO_RayStorageBaseGeometry			(std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords);
		virtual									   ~SKTRANSO_RayStorageBaseGeometry			();

};


