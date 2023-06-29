#pragma once


/*-----------------------------------------------------------------------------
*					SKTRAN_TableOpticalProperties_PointCache		 2016-05-17*/
/** @ingroup optprop
*  Allows optical properties to be cached for a single point, so that lookups
*  into all but the scatter angle grid can be avoided. 
**/
/*---------------------------------------------------------------------------*/
class SKTRAN_TableOpticalProperties_PointCache
{
    private:
    std::vector<SKTRAN_ScatMat_MIMSNC>        m_scatterMatricesPerMeter;
    const SKTRAN_GridDefScatterAngle_V21*   m_scatteranglegrid;			//!< The grid that defines the angular grid for specifying scattered rays.
    double m_singleScatterAlbedo;
    HELIODETIC_POINT m_fillPoint;

    private:
    void Release(){ 
        if(nullptr!=m_scatteranglegrid) m_scatteranglegrid->Release(); 
        m_scatteranglegrid=nullptr; 
    }

    public:
    SKTRAN_TableOpticalProperties_PointCache ( ) { m_scatteranglegrid=nullptr; }
    ~SKTRAN_TableOpticalProperties_PointCache ( ) { Release(); }

    void InjectTable( std::vector<SKTRAN_ScatMat_MIMSNC>& scatmatsPerMeter, const SKTRAN_GridDefScatterAngle_V21* anglegrid, double singleScatterAlbedo, const HELIODETIC_POINT& fillPoint ) // Destroys #scatmats 
    {
        m_scatterMatricesPerMeter.swap(scatmatsPerMeter);
        anglegrid->AddRef();
        if(nullptr!=m_scatteranglegrid) m_scatteranglegrid->Release();
        m_scatteranglegrid = anglegrid;

        m_singleScatterAlbedo = singleScatterAlbedo;
        m_fillPoint                 = fillPoint;
    }

    bool Interpolate_kscattPerM (double cosScattAngle, SKTRAN_ScatMat_MIMSNC& scatmat) const {
        bool ok = true;
        SKTRAN_GridIndex loAngIndex, hiAngIndex;
        double loAngWeight, hiAngWeight;
        ok = ok && m_scatteranglegrid->FindBoundingIndices(cosScattAngle, SKTRAN_GridDefBase_V2::OUTOFBOUND_ERROR, &loAngIndex, &loAngWeight, &hiAngIndex, &hiAngWeight );
        scatmat.SetTo(0.0);
        scatmat.AddToThis( m_scatterMatricesPerMeter[loAngIndex], loAngWeight);
        scatmat.AddToThis( m_scatterMatricesPerMeter[hiAngIndex], hiAngWeight);
        return ok;
    }

    double GetSingleScatterAlbedo        ( ) const { return m_singleScatterAlbedo; }
    const HELIODETIC_POINT& GetFillPoint ( ) const { return m_fillPoint;           }

};