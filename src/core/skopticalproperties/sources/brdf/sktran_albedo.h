#pragma once


/*-----------------------------------------------------------------------------
 *					class skRTAlbedo 					2008-3-20*/
/** \ingroup skalbedo
 *	A class for retrieving the albedo at a given wavelength and geographic
 *	variation. This is still in development. The call to GetAlbedo must be
 *	thread safe for Sasktran version 2.1.
**/
/*---------------------------------------------------------------------------*/

class skRTAlbedo_DEPRECATED : public nxUnknown
{
	public:
		virtual bool			GetAlbedo( double wavelennm, const GEODETIC_INSTANT& pt, double* albedo)  const= 0; /* MUST BE Thread safe for Sasktran version 2.1 */
		virtual bool			CreateClone( skRTAlbedo_DEPRECATED** clone ) const  = 0;
		virtual bool			IsConstantAlbedo( ) const = 0;
};


/*-----------------------------------------------------------------------------
 *					skRTAlbedo_Constant		2008-3-20*/
/** \ingroup skalbedo
	A class that implements a "grey/white" albedo that is wavelength independent albedo
 */
/*---------------------------------------------------------------------------*/

class skRTAlbedo_Constant_DEPRECATE : public skRTAlbedo_DEPRECATED
{
	private:
//		static size_t			m_numinstances;
		double					m_albedo;

	public:
								skRTAlbedo_Constant_DEPRECATE()															{ /*m_numinstances++;*/ m_albedo = 0;}
								skRTAlbedo_Constant_DEPRECATE( double albedo )											{ /*m_numinstances++;*/ m_albedo = 0; SetAlbedo(albedo);}
		virtual				   ~skRTAlbedo_Constant_DEPRECATE()															{ /*m_numinstances--;*/}
//		static size_t			NumInstances ()																	{ return m_numinstances;}
		bool					SetAlbedo( double albedo )														{ m_albedo = albedo; return true;}
		virtual bool			GetAlbedo( double /*wavelennm*/, const GEODETIC_INSTANT& /*pt*/, double* albedo) const	{ *albedo = m_albedo; return true;}
		virtual bool			CreateClone( skRTAlbedo_DEPRECATED** clone ) const;
		virtual bool			IsConstantAlbedo( ) const {return true;}


};

class skRTAlbedo_Plane_DEPRECATE : public skRTAlbedo_DEPRECATED
{
	private:
		skClimatology*			m_clim;
	public:
			skRTAlbedo_Plane_DEPRECATE ( skClimatology* clim )
			{
				m_clim = clim;
				clim->AddRef();
			};
			virtual ~skRTAlbedo_Plane_DEPRECATE()
			{
				m_clim->Release();
			};
			bool GetAlbedo( double /*wavelennm*/, const GEODETIC_INSTANT& pt, double* albedo ) const
			{
				return m_clim->GetParameter( SKCLIMATOLOGY_ALBEDO, pt, albedo, true );
			};
			bool CreateClone( skRTAlbedo_DEPRECATED** /*clone*/ ) const { return false; }
			bool IsConstantAlbedo() const { return true; }

};



/*-----------------------------------------------------------------------------
 *					skRTBRDF_AlbedoPlane		 2016- 12- 21*/
/** **/
/*---------------------------------------------------------------------------*/

class skBRDF_AlbedoPlane : public skBRDF
{
	private:
		skClimatology*			m_clim;
	public:
			skBRDF_AlbedoPlane ( skClimatology* clim )
			{
				m_clim = clim;
				clim->AddRef();
			};
			virtual ~skBRDF_AlbedoPlane()
			{
				m_clim->Release();
			};

			bool BRDF( double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double DPHI, double* brdf) const
			{
				double	albedo;
				bool	ok;

				ok = m_clim->GetParameter( SKCLIMATOLOGY_ALBEDO, pt, &albedo, false );
				*brdf = albedo/nxmath::Pi;
				return ok;
			};

			virtual bool	IsLambertian() const override {
				return true;
			}
};

class SKTRAN_BRDF_UserDefinedLatLon : public skBRDF
{
	private:
		nx2dArray<skBRDF*> m_brdfs;

		//SKTRAN_GridDefBase_V2 m_latitudes;
		//SKTRAN_GridDefBase_V2 m_longitudes;

		std::vector<double>	  m_latitudes;
		std::vector<double>	  m_longitudes;

		void LatitudeInterpolate(double latitude, std::array<size_t, 2>& indicies, std::array<double, 2>& weights, int& numweights) const;
		void LongitudeInterpolate(double longitude, std::array<size_t, 2>& indicies, std::array<double, 2>& weights, int& numweights) const;

	public:
		SKTRAN_BRDF_UserDefinedLatLon(const std::vector<double>& latitudes, const std::vector<double>& longitudes, const nx2dArray<skBRDF*>& brdfs);
		SKTRAN_BRDF_UserDefinedLatLon();
		virtual ~SKTRAN_BRDF_UserDefinedLatLon();

		bool Assign(const std::vector<double>& latitudes, const std::vector<double>& longitudes, const nx2dArray<skBRDF*>& brdfs);


		bool BRDF(double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double DPHI, double* brdf) const;

		virtual bool IsLambertian() const override {
			return false;
		}
};

class SKTRAN_BRDF_SpectralVarying : public skBRDF
{
private:
	nx1dArray<skBRDF*> m_brdfs;

	std::vector<double> m_wavelengths;

	void Interpolate(double wavelength, std::array<size_t, 2>& indicies, std::array<double, 2>& weights, int& numweights) const;

public:
	SKTRAN_BRDF_SpectralVarying() {};
	SKTRAN_BRDF_SpectralVarying(const std::vector<double>& wavelengths, const nx1dArray<skBRDF*>& brdfs);
	virtual ~SKTRAN_BRDF_SpectralVarying();

	bool Assign(const std::vector<double>& wavelengths, const nx1dArray<skBRDF*>& brdfs);

	bool BRDF(double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double DPHI, double* brdf) const;

	virtual bool IsLambertian() const override {
		bool lambertian = false;
		for (auto& brdf : m_brdfs) {
			lambertian = lambertian && brdf->IsLambertian();
		}
		return lambertian;
	}
};


/*-----------------------------------------------------------------------------
 *					skRTAlbedo_Variable		2010-3-18*/
/**  \ingroup skalbedo
 *	A class that provides linear interpolation of the albedo with
 *	wavelength
 **/
/*---------------------------------------------------------------------------*/

class skRTAlbedo_Variable_DEPRECATE : public skRTAlbedo_DEPRECATED
{

	private:
		nx1dArray<double>		m_wavelengths;
		nx1dArray<double>		m_albedo;

	public:
								skRTAlbedo_Variable_DEPRECATE(){}
		virtual				   ~skRTAlbedo_Variable_DEPRECATE(){}
		bool					SetAlbedo( const double* albedo, const double* wavelennm, size_t npts );
		virtual bool			GetAlbedo( double wavelennm, const GEODETIC_INSTANT& pt, double* albedo) const;
		virtual bool			CreateClone( skRTAlbedo_DEPRECATED** clone ) const;
		virtual bool			IsConstantAlbedo( ) const {return false;}
};

/*-----------------------------------------------------------------------------
 *					skRTAlbedo_Variable		2010-3-18*/
/**  \ingroup skalbedo
 *	A class that provides linear interpolation of the albedo with
 *	wavelength
 **/
/*---------------------------------------------------------------------------*/

class skBRDF_VariableAlbedo : public skBRDF
{

	private:
		nx1dArray<double>		m_wavelengths;
		nx1dArray<double>		m_albedo;

	public:
								skBRDF_VariableAlbedo(){}
		virtual				   ~skBRDF_VariableAlbedo(){}
		bool					SetAlbedo( const double* albedo, const double* wavelennm, size_t npts );
//		bool					GetAlbedo( double wavelennm, const GEODETIC_INSTANT& pt, double* albedo) const;
//		virtual bool			CreateClone( skRTAlbedo_DEPRECATED** clone ) const;
//		virtual bool			IsConstantAlbedo( ) const {return false;}
		virtual bool			BRDF( double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* brdf)  const override; /* MUST BE Thread safe for Sasktran version 2.1 */

		virtual bool	IsLambertian() const override {
			return true;
		}

};



