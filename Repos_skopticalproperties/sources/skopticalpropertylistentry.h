#pragma once


/*---------------------------------------------------------------------------
 *            skOpticalProperties_ListEntry            2020-01-02 */
/** **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_ListEntry
{
	public:
		skOpticalProperties*	m_opticalproperty;
		double					m_startwavenumber;
		double					m_endwavenumber;

	public:
								skOpticalProperties_ListEntry( double startwavenumber, double endwavenumber,  skOpticalProperties* optprop);
								skOpticalProperties_ListEntry( const skOpticalProperties_ListEntry& other );
							   ~skOpticalProperties_ListEntry();
		bool					operator <					(double wavenum) const { return m_endwavenumber < wavenum;}
};


/*---------------------------------------------------------------------------
 *            Class skOpticalProperties_ListEntries            2019-12-05 */
/** A class used to create a single optical property object from a list
 *  of individual optical property objects. This is useful for implementing
 *	cross-sections of species that have bbeen measured in diffrent spectral
 *	regions with diffrent instruments. For example the Hitran O4 cross-section.
 *
 *	We currently support only absorbing species with non-overlapping 
 *	spectral regions
 **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_ListEntries : public skOpticalProperties
{
	private:
		        std::vector< skOpticalProperties_ListEntry >					m_entries;
		typedef std::vector< skOpticalProperties_ListEntry >::iterator			iterator;

	private:
											skOpticalProperties_ListEntries					(const skOpticalProperties_ListEntries& other);	// Dont allow copy constructor
											skOpticalProperties_ListEntries&	operator =	(const skOpticalProperties_ListEntries& other);	// Dont allow assignment operator
		skOpticalProperties*				FindEntry										( double wavenum );
		bool								CheckNonOverlapping								( double wavenum );

	public:
											skOpticalProperties_ListEntries					();
		virtual							   ~skOpticalProperties_ListEntries					() override;
		bool								AddEntry										( double startwavenum, double endwavenum, skOpticalProperties* optprop);

	public:
		virtual bool						SetAtmosphericState								( skClimatology* neutralatmosphere) override;
		virtual bool						SetLocation										( const GEODETIC_INSTANT& pt, bool* crosssectionschanged ) override;
		virtual bool						InternalClimatology_UpdateCache					( const GEODETIC_INSTANT& pt) override;
		virtual bool						CalculateCrossSections							( double wavenumber, double* absxs, double* extxs, double* scattxs ) override;
		virtual bool						CalculatePhaseMatrix							( double wavenumber, double cosscatterangle, skRTPhaseMatrix* phasematrix) override;
		virtual bool						IsScatterer										() const override	{ return false;}
		virtual bool						IsAbsorber										() const override	{ return true;}
		virtual double						DeltaFunctionForwardScatterFraction				() const			{ return 0.0;}

};


