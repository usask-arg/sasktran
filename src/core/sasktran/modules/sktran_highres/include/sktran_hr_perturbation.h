//#include "sktran_hr_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Perturbation_Base		2014-10-30*/
/** Defines a singular perturbation in the atmosphere in which a weighting
 *  function is to be calculated from.  
 *
 *  For example, a perturbation could be a linear increase of absorption extinction
 *  from 11km to 12km with the peak at 11.5km, or it could be a constant increase
 *  of absorption between 11km and 12km.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_HR_Perturbation_Base
{
	private:
		double m_pertval;
	private:
	public:
		SKTRAN_HR_Perturbation_Base() { SetPertVal(1.0); };
		virtual ~SKTRAN_HR_Perturbation_Base() { };

		void SetPertVal(double val) { m_pertval = val; }
		double GetPertVal() const { return m_pertval; }

		virtual HELIODETIC_POINT PerturbationLocation( const SKTRAN_CoordinateTransform_V2& coords ) const = 0;
		virtual double PerturbationAltitudeWidth() const = 0;
		virtual double PerturbationAltitudeLower() const = 0;
		virtual double PerturbationAltitudeUpper() const = 0;
		virtual size_t NumBoundingGeometry() const = 0;
		virtual std::unique_ptr< SKTRAN_GeometryObject > BoundingGeometryObject( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords, size_t idx ) const = 0;
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Perturbation_Absorption		2014-10-30*/
/** Defines a constant increase of an atmospheric paramater between
 *  two altitudes.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_HR_Perturbation_Absorption : public SKTRAN_HR_Perturbation_Base
{
	private:
		double			m_upperbound;
		double			m_lowerbound;
	public:
		SKTRAN_HR_Perturbation_Absorption() { };
		virtual ~SKTRAN_HR_Perturbation_Absorption() { };

		bool		 Initialize( double upperheight,
								 double lowerheight,
								 double pertval );

		virtual bool ExtinctionPerturbation( const HELIODETIC_POINT& location,
												   bool& isperturbation,
												   double& value) const;

		virtual HELIODETIC_POINT PerturbationLocation( const SKTRAN_CoordinateTransform_V2& coords ) const;
		virtual double PerturbationAltitudeWidth() const { return m_upperbound - m_lowerbound; }
		virtual double PerturbationAltitudeLower() const { return m_lowerbound; }
		virtual double PerturbationAltitudeUpper() const { return m_upperbound; }
		virtual size_t NumBoundingGeometry() const override { return 3; };
		virtual std::unique_ptr< SKTRAN_GeometryObject > BoundingGeometryObject( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords, size_t idx ) const override;
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Perturbation_Linear		2014-10-30*/
/**  Defines a perturbation centered at an altitude which decreases linearly in radius
 *   from the center value to a value of 0.  This best models SASKTRAN calculations
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_HR_Perturbation_Absorption_Linear : public SKTRAN_HR_Perturbation_Base
{
	private:
		double		m_center;
		double		m_distancetozeroabove;
		double		m_distancetozerobelow;
	public:
		bool		Initialize( double center,
								double distancetozerobelow,
								double distancetozeroabove,
								double pertval );
		virtual bool ExtinctionPerturbation( const HELIODETIC_POINT& location,
												   bool& isperturbation,
												   double& value) const;
		virtual HELIODETIC_POINT PerturbationLocation( const SKTRAN_CoordinateTransform_V2& coords ) const;

		virtual size_t NumBoundingGeometry() const override { return 3; };
		virtual double PerturbationAltitudeWidth() const { return m_distancetozeroabove + m_distancetozerobelow; }
		virtual double PerturbationAltitudeLower() const { return m_center - m_distancetozerobelow; }
		virtual double PerturbationAltitudeUpper() const { return m_center + m_distancetozeroabove; }
		virtual std::unique_ptr< SKTRAN_GeometryObject > BoundingGeometryObject( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords, size_t idx ) const override;
		
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Perturbation_Absorption_Box		2014-10-30*/
/** Defines a perturbation bounded by two spherical shells and two planes,
 *  interpolation is bilinear from the center of the box in altitude and angle
 *  in the normal direction of the planes
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_HR_Perturbation_Absorption_Box : public SKTRAN_HR_Perturbation_Base
{
	private:
		double		m_centerradius;
		double		m_radiusdistancetozero;
		double		m_anglebetweennormals;
		nxVector	m_normal1;
		nxVector	m_normal2;
		
	public:
		bool		 Initialize( double center,
								 double distancetozero,
								 const nxVector& normal1,
								 const nxVector& normal2,
								 double pertval );

		bool ExtinctionPerturbation( const HELIODETIC_POINT& location,
												   bool& isperturbation,
												   double& value) const
		{
			bool ok = true;

			double alt = location.Altitude();
			double disttocenter = abs(alt - m_centerradius);

			if (disttocenter > m_radiusdistancetozero)
			{
				isperturbation = false;
				value = 0.0;
				return ok;
			}

			const HELIODETIC_UNITVECTOR& heliounit = location.Vector().UnitVector();
			nxVector locunit(heliounit.X(), heliounit.Y(), heliounit.Z());

			double angleToNormal1 = abs(abs(nxmath::acosd(locunit & m_normal1)) - 90);

			if (angleToNormal1 <= m_anglebetweennormals && abs(abs(nxmath::acosd(locunit & m_normal2)) - 90) <= m_anglebetweennormals)
			{
				isperturbation = true;
				value = GetPertVal() * (1 - (disttocenter / m_radiusdistancetozero));
				// linear interpolation in angle from the middle
				value *= (m_anglebetweennormals / 2 - abs(angleToNormal1 - m_anglebetweennormals / 2)) / (m_anglebetweennormals / 2);
			}
			else
			{
				isperturbation = false;
				value = 0.0;
			}

			return ok;
		}
		virtual HELIODETIC_POINT PerturbationLocation( const SKTRAN_CoordinateTransform_V2& coords ) const;
		virtual double PerturbationAltitudeWidth() const { return 2*m_radiusdistancetozero; }
		double PerturbationCenterAltitude() const { return m_centerradius; }
		virtual double PerturbationAltitudeLower() const { return m_centerradius - m_radiusdistancetozero; }
		virtual double PerturbationAltitudeUpper() const { return m_centerradius + m_radiusdistancetozero; }

		virtual size_t NumBoundingGeometry() const override { return 4; };
		virtual std::unique_ptr< SKTRAN_GeometryObject > BoundingGeometryObject( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords, size_t idx ) const override;
};