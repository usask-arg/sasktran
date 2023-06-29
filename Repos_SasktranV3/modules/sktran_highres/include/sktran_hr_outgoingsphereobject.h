//#include "sktran_hr_internals.h"

//class SKTRAN_HR_OutgoingSphereObject_Base;
//class SKTRAN_HR_OutgoingSphereObject_NoRotation;
//class SKTRAN_HR_OutgoingSphereObject_AziRotation;

class SKTRAN_HR_OutgoingSphereObject_Base : public nxUnknown
{
	protected:
		const SKTRAN_UnitSphere_V2*							m_outgoingsphere;

	private:
		void ReleaseResources();

	public:
		         SKTRAN_HR_OutgoingSphereObject_Base ( );
		virtual ~SKTRAN_HR_OutgoingSphereObject_Base ( );
		void     SetOutgoingSphere                   ( const SKTRAN_UnitSphere_V2* sphere );
		bool     TriangulateOnOutgoing               ( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t maxvertices) const;
		size_t   NumOutgoingRays                     ( ) const;
		void     OutgoingRayLocalCoords                         ( size_t idx, nxVector& outray ) const;
		double   OutgoingCubatureWeight              ( size_t outidx ) const;
};

//class SKTRAN_HR_OutgoingSphereObject_NoRotation;
//class SKTRAN_HR_OutgoingSphereObject_AziRotation;
