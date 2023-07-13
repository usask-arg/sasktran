class SKTRAN_RayTracer_Curved_Shells;

class SKTRAN_RayOptical_Curved : public SKTRAN_RayOptical_StraightQuadrature_Base
{
private:
	std::unique_ptr<SKTRAN_RayStorage_CurvedPiecewise>				m_trajectorystorageobject;	
	SKTRAN_RayStorage_CurvedPiecewise*								m_trajectorystorage;		
	std::shared_ptr<SKTRAN_RayTracer_Curved_Shells>					m_raytracer;


public:
	SKTRAN_RayOptical_Curved(std::unique_ptr<SKTRAN_RayStorage_CurvedPiecewise> trajectorystorage, std::shared_ptr< SKTRAN_RayTracer_Curved_Shells> raytracer);
	virtual											   ~SKTRAN_RayOptical_Curved();
	SKTRAN_RayStorage_CurvedPiecewise*					StorageVar() { return m_trajectorystorage; }			//!< Straight rays need a pointer to staright ray storage but the base class holds the unique_ptr

	bool												LocationAlongRayAsVector(const double& distancealongray, HELIODETIC_VECTOR* pt) const;
	bool												GetQuadratureInterpParams_startPoint(size_t quadraturepoint, size_t cellidx, double* r0, double* t0, double* rt, HELIODETIC_POINT* startpoint) const;
	bool												GetQuadratureInterpParams(size_t cellidx, double* r0, double* r1, double* t0, double* t1, double* rt,
		HELIODETIC_POINT* startpoint, HELIODETIC_POINT* endpoint) const;

	virtual bool										TraceRay_NewMethod() override;
	virtual void										NotifyDerived_RayInvalid() override;


};