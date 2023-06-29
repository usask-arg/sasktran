//#include "sktran_common_internals.h"

class SKTRAN_MCPhoton_Base;

/*-----------------------------------------------------------------------------
 *					SKTRAN_SourceTermQueryObject_Base		 2016- 7- 12*/
/** @ingroup stokes
 *  A class that manages the different axes associated with a ray of light. It includes
 *	the following axes, unit vectors and vectors:
 *	-# The look unit vector of the ray away from the observer/origin/diffuse point. See method GetLook
 *	-# The vector location of this point along the ray trajectory. See method GetPoint
 *	-# The three basis vectors of this ray, See method GetBasis()
 *		-# basis[0] = unit vector towards observer (anti-parallel to method GetLook())
 *		-# basis[1] = unit vector perpendicular to ray propagation
 *		-# basis[2] = unit vector forming rhs orthogonal systemm.
 *
 *	See: Page 25 Mischenko, Scattering, Absorption and Emission of Light by Small particles. 
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_SourceTermQueryObject_Base
{

    protected:
		SKTRAN_SourceTermQueryObject_Base (){;}

    public:
		SKTRAN_SourceTermQueryObject_Base ( const HELIODETIC_POINT& observer, const HELIODETIC_UNITVECTOR& look ) 
		{
//			NXASSERT( (observer.Altitude() < 150000) );
		}
		virtual ~SKTRAN_SourceTermQueryObject_Base() {;}

        virtual void  UpdateQuery ( const HELIODETIC_POINT& pt, const HELIODETIC_UNITVECTOR& look ) = 0;
		virtual const HELIODETIC_POINT&      GetPoint   () const = 0; 
		virtual const HELIODETIC_UNITVECTOR& GetLookAway()   const = 0;
        virtual const HELIODETIC_BASIS       GetBasis   ()  const = 0; // Let RVO save copy

};



/*-----------------------------------------------------------------------------
 *					SKTRAN_SourceTermQueryObject_Simple		 2016- 7- 12*/
/**	 @ingroup stokes
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_SourceTermQueryObject_Simple : public SKTRAN_SourceTermQueryObject_Base
{

    protected:
		SKTRAN_SourceTermQueryObject_Simple () : SKTRAN_SourceTermQueryObject_Base() {;}

    private:
        HELIODETIC_POINT      m_pt;
        HELIODETIC_UNITVECTOR m_look;

    public:
        SKTRAN_SourceTermQueryObject_Simple( const HELIODETIC_POINT& observer, const HELIODETIC_UNITVECTOR& look )
			: SKTRAN_SourceTermQueryObject_Base ( observer, look ), m_pt(observer), m_look(look) {;}
		virtual ~SKTRAN_SourceTermQueryObject_Simple(){;}

		virtual void  UpdateQuery ( const HELIODETIC_POINT& pt, const HELIODETIC_UNITVECTOR& look ) override { m_pt=pt; m_look=look;}
		virtual const HELIODETIC_POINT&      GetPoint   () const override   { return m_pt;}
		virtual const HELIODETIC_UNITVECTOR& GetLookAway() const override { return m_look;}
		virtual const HELIODETIC_BASIS       GetBasis   () const override  { return  HELIODETIC_BASIS(m_pt,m_look);}
};



/*-----------------------------------------------------------------------------
 *					SKTRAN_SourceTermQueryObject_StraightPolarized		 2016- 7- 12*/
/** @ingroup stokes
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_SourceTermQueryObject_StraightPolarized : public SKTRAN_SourceTermQueryObject_Base
{

    protected:
		SKTRAN_SourceTermQueryObject_StraightPolarized () : SKTRAN_SourceTermQueryObject_Base() {;}

    protected:
        HELIODETIC_BASIS       m_basis;
        HELIODETIC_POINT       m_pt;
        HELIODETIC_UNITVECTOR  m_look;

    public:
        SKTRAN_SourceTermQueryObject_StraightPolarized( const HELIODETIC_POINT& observer, const HELIODETIC_UNITVECTOR& look )
			: SKTRAN_SourceTermQueryObject_Base ( observer, look ), m_basis(observer,look), m_pt(observer), m_look(look) 
		{
		;
		}
		virtual ~SKTRAN_SourceTermQueryObject_StraightPolarized(){;}

		virtual void  UpdateQuery ( const HELIODETIC_POINT& pt, const HELIODETIC_UNITVECTOR& look ) override { m_pt=pt; m_look=look;}
		virtual const HELIODETIC_POINT&      GetPoint () const override { return m_pt;}
		virtual const HELIODETIC_UNITVECTOR& GetLookAway()   const override { return m_look;}
		virtual const HELIODETIC_BASIS       GetBasis()  const override { return  m_basis;}
};



/*-----------------------------------------------------------------------------
 *					SKTRAN_SourceTermQueryObject_SolarScatterPolarized		 2016- 7- 12*/
/** @ingroup stokes
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_SourceTermQueryObject_SolarScatterPolarized : public SKTRAN_SourceTermQueryObject_Base
{

    protected:
		SKTRAN_SourceTermQueryObject_SolarScatterPolarized () : SKTRAN_SourceTermQueryObject_Base() {;}
		void ProduceBasisFromCurrentGeometry( ){
			m_basis.x.SetCoords( -m_look.X(), -m_look.Y(), -m_look.Z() );
			HELIODETIC_VECTOR temp;
			temp.SetCoords( -m_look.Z()*m_look.X(), -m_look.Y()*m_look.Z(), 1.0 - m_look.Z()*m_look.Z() );
			if(1e-6<temp.Magnitude()){
				m_basis.y = temp.UnitVector();
				temp.SetCoords(m_basis.x.Y()*m_basis.y.Z() - m_basis.x.Z()*m_basis.y.Y(), m_basis.x.Z()*m_basis.y.X() - m_basis.x.X()*m_basis.y.Z(), m_basis.x.X()*m_basis.y.Y() - m_basis.x.Y()*m_basis.y.X() );
				m_basis.z = temp.UnitVector();
			} else{
				m_basis.ProduceBasis( m_pt, m_look );
			}	
		}
		
		size_t BranchUsedForLook( const HELIODETIC_UNITVECTOR& look ) const{
			size_t bidx = 0;
			HELIODETIC_VECTOR temp;
			temp.SetCoords( -look.Z()*look.X(), -look.Y()*look.Z(), 1.0 - pow(look.Z(),2.0) );
			if(1e-6<temp.Magnitude()){
					bidx = 0;
				} else{
					bidx = 1;
				}	
			return bidx;
		}
		void ProduceBasisFromCurrentGeometry_UsingBranch ( size_t bidx ){
			m_basis.x.SetCoords( -m_look.X(), -m_look.Y(), -m_look.Z() );
			HELIODETIC_VECTOR temp;
			temp.SetCoords( -m_look.Z()*m_look.X(), -m_look.Y()*m_look.Z(), 1.0 - m_look.Z()*m_look.Z() );
			switch(bidx){
			case 0:
				if(1e-6<temp.Magnitude()){
					m_basis.y = temp.UnitVector();
					temp.SetCoords(m_basis.x.Y()*m_basis.y.Z() - m_basis.x.Z()*m_basis.y.Y(), m_basis.x.Z()*m_basis.y.X() - m_basis.x.X()*m_basis.y.Z(), m_basis.x.X()*m_basis.y.Y() - m_basis.x.Y()*m_basis.y.X() );
					m_basis.z = temp.UnitVector();
				} else{
					//nxLog::Record(NXLOG_WARNING, "SKTRAN_SourceTermQueryObject_SolarScatterPolarized::ProduceBasisFromCurrentGeometry_UsingBranch, Unsafe branch, mag= %1.16e.", temp.Magnitude());
					m_basis.ProduceBasis( m_pt, m_look );
				}
				break;
			case 1:
				m_basis.ProduceBasis( m_pt, m_look );
				break;
			default:
				nxLog::Record(NXLOG_ERROR, "SKTRAN_SourceTermQueryObject_SolarScatterPolarized::ProduceBasisFromCurrentGeometry_UsingBranch, Unrecognized branch id %u.", bidx);
			}	
		}

    protected:
        HELIODETIC_BASIS       m_basis;
        HELIODETIC_POINT       m_pt;
        HELIODETIC_UNITVECTOR  m_look;

    public:
        SKTRAN_SourceTermQueryObject_SolarScatterPolarized( const HELIODETIC_POINT& observer, const HELIODETIC_UNITVECTOR& look )
			: SKTRAN_SourceTermQueryObject_Base ( observer, look ), m_pt(observer), m_look(look) 
			{ 
				ProduceBasisFromCurrentGeometry( );
		    }
		virtual ~SKTRAN_SourceTermQueryObject_SolarScatterPolarized(){;}
		
		virtual void  UpdateQuery ( const HELIODETIC_POINT& pt, const HELIODETIC_UNITVECTOR& look ) override { m_pt=pt; m_look=look; ProduceBasisFromCurrentGeometry();}
		virtual void  UpdateQuery_SameBranch ( const HELIODETIC_POINT& pt, const HELIODETIC_UNITVECTOR& look, const HELIODETIC_UNITVECTOR& lookToMatchBranch ) {
			m_pt=pt;
			m_look=look;
			size_t bidx = BranchUsedForLook( lookToMatchBranch );
			ProduceBasisFromCurrentGeometry_UsingBranch( bidx );
		}
		virtual const HELIODETIC_POINT&      GetPoint () const override { return m_pt;}
		virtual const HELIODETIC_UNITVECTOR& GetLookAway()   const override { return m_look;}
		virtual const HELIODETIC_BASIS       GetBasis()  const override { return  m_basis;}
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_SourceTermQueryObject_ModifiablePolarized		 2016- 7- 12*/
/** @ingroup stokes
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_SourceTermQueryObject_ModifiablePolarized : public SKTRAN_SourceTermQueryObject_StraightPolarized
{

    protected:
		SKTRAN_SourceTermQueryObject_ModifiablePolarized () : SKTRAN_SourceTermQueryObject_StraightPolarized() {;}
        
    public:
        SKTRAN_SourceTermQueryObject_ModifiablePolarized( const HELIODETIC_POINT& observer, const HELIODETIC_UNITVECTOR& look )
			: SKTRAN_SourceTermQueryObject_StraightPolarized( observer, look ) {;}
		virtual ~SKTRAN_SourceTermQueryObject_ModifiablePolarized(){;}

		void UpdateBasis( const HELIODETIC_BASIS& basis ) { m_basis = basis;}
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_GeographicBasisGenerator		 2016- 7- 12*/
/** @ingroup stokes
**/
/*---------------------------------------------------------------------------*/

template< typename queryObjectType >
class SKTRAN_GeographicBasisGenerator
{
	public:
		static bool Generate( const nxVector& obs, const nxVector& look, GEOGRAPHIC_BASIS& gb ) {
			HELIODETIC_VECTOR     v;
			HELIODETIC_POINT      p;
			HELIODETIC_UNITVECTOR u;
			v.SetCoords( obs.X(), obs.Y(), obs.Z() );
			p.FromVector( v, nullptr );
			u.SetCoords( look.X(), look.Y(), look.Z() );
			HELIODETIC_BASIS temp = queryObjectType( p, u ).GetBasis();
			gb =  GEOGRAPHIC_BASIS(  nxVector(temp.X().X(), temp.X().Y(), temp.X().Z()), 
								  	 nxVector(temp.Y().X(), temp.Y().Y(), temp.Y().Z()), 
									 nxVector(temp.Z().X(), temp.Z().Y(), temp.Z().Z()) );
			return true;
		}
};

template<>
class SKTRAN_GeographicBasisGenerator< SKTRAN_SourceTermQueryObject_SolarScatterPolarized >
{
	public:
		static bool Generate( const nxVector& obs, const nxVector& look, GEOGRAPHIC_BASIS& gb ) {
			nxLog::Record( NXLOG_WARNING, "SKTRAN_GeographicBasisGenerator< SKTRAN_SourceTermQueryObject_SolarScatterPolarized >::Generate, Would need to have an implementation that receives coords (or sun direction) to make this basis.");
			return false;
		}
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_Source_Term		2014-2-7*/
/** @ingroup stokes
**/
/*---------------------------------------------------------------------------*/


class SKTRAN_Source_Term : public nxUnknown
{
	private:

	public:
		
		/* The argument #ray provides a coordinate system and, for the polarized functions, a basis for the polarization components. Ray should maybe provide a "CoordInfo" class that
		   gets passed around instead -- it's not clear here whether more information is needed from the ray. */


		/*-----------------------------------------------------------------------------
		 *					SourceTermAtPoint		 2016- 11- 25*/
		/** Calculate the "Source Term" radiance at the point location defined by qobj
		 *	in the look direction away from the observer/diffuse point specified by qobj. 
		 *	Return radiance in photons/cm2/sec/steradian/nm possibly normalized to solar irradiance
		**/
		/*---------------------------------------------------------------------------*/

		virtual bool			SourceTermAtPoint			   ( const SKTRAN_SourceTermQueryObject_Base& qobj, double*           source   ) const = 0;
		virtual bool			SourceTermAtPoint			   ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC* source   ) const = 0;

		virtual bool			GroundSourceAtPoint			   ( const SKTRAN_SourceTermQueryObject_Base& qobj, double*           source   ) const = 0;
		virtual bool			GroundSourceAtPoint			   ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC* source   ) const = 0;
		
		virtual bool            MonteCarlo_SingleScatteredRadianceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, double&           radiance ) const { return false;}
		virtual bool            MonteCarlo_SingleScatteredRadianceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC& radiance ) const { return false;}

		virtual bool            MonteCarlo_GroundScatteredRadianceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, double&           radiance ) const { return false;}
		virtual bool            MonteCarlo_GroundScatteredRadianceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC& radiance ) const { return false;}


		// these four should probably be pure virtual functions like their wavelength-independent counterparts if we implement wavelength dependence everywhere
		virtual bool			SourceTermAtPoint( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, double*           source) const { return false; }
		virtual bool			SourceTermAtPoint( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC* source) const { return false; }

		virtual bool			GroundSourceAtPoint( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, double*           source) const { return false; }
		virtual bool			GroundSourceAtPoint( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC* source) const { return false; }

 		virtual bool			MonteCarlo_SingleScatteredRadianceAtPoint( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, double&		   radiance) const { return false; }
		virtual bool			MonteCarlo_SingleScatteredRadianceAtPoint( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC& radiance) const { return false; }

		virtual bool            MonteCarlo_GroundScatteredRadianceAtPoint( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, double&           radiance) const { return false; }
		virtual bool            MonteCarlo_GroundScatteredRadianceAtPoint( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC& radiance) const { return false; }

		// these could potentially replace all MonteCarlo overloads
		virtual bool			MonteCarlo_SingleScatteredRadianceAtPoint( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_MCPhoton_Base* photon ) const { return false; }
		virtual bool			MonteCarlo_GroundScatteredRadianceAtPoint( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_MCPhoton_Base* photon ) const { return false; }
};

