#include <skopticalproperties21.h>

//#include "pch.h"
//#include "hitranline_upperstates.h"


struct lessthan_Eupper : public binary_function <const skSpectralLine_HitranLine*, const skSpectralLine_HitranLine*, bool>
{
    bool operator()( const skSpectralLine_HitranLine* Left, const skSpectralLine_HitranLine* Right) const 
	{ 
		return Left->EUpper() < Right->EUpper(); 
	}

};


/*---------------------------------------------------------------------------
 *            HitranLine_UpperState::BoltzmannExponent            2020-07-24 */
/** **/
/*---------------------------------------------------------------------------*/

double	HitranLine_UpperState::BoltzmannExponent( double T) const
{
	double a = nxSI::HCOK*100.0;			// Multiply by 100 to convert wavenumber to m-1 from cm-1
	double	frac;
	double	val;

	frac = a*m_E/T;
	val  = m_degeneracy*exp( -frac );
	return val;
}


/*---------------------------------------------------------------------------
 *        HitranLine_UpperStates::~HitranLine_UpperStates         2020-07-28 */
/** **/
/*---------------------------------------------------------------------------*/

HitranLine_UpperStates::~HitranLine_UpperStates()
{
	for ( auto iter = m_upperstates.begin(); iter != m_upperstates.end(); ++iter)
	{
		(*iter)->Release();
	}
};


/*---------------------------------------------------------------------------
 *          HitranLine_UpperStates::CreateNewUpperState           2020-07-24 */
/** **/
/*---------------------------------------------------------------------------*/

bool HitranLine_UpperStates::CreateNewUpperState( HitranLine_UpperState** newupperstate, double E, double degeneracy )
{
	HitranLine_UpperState*	entry;

	entry = new  HitranLine_UpperState(E, degeneracy);
	entry->AddRef();
	m_upperstates.push_back( entry );
	*newupperstate = entry;
	return true;
}


/*---------------------------------------------------------------------------
 *            HitranLine_UpperStates::UpdatePartitions            2020-07-24 */
/** **/
/*---------------------------------------------------------------------------*/

bool HitranLine_UpperStates::UpdatePartitions( double T)
{
	double	Z = 0.0;
	double	fraction;

	if (T != m_current_temperature)
	{
		for (auto & entry : m_upperstates)
		{
			Z  += entry->BoltzmannExponent(T);
		}
		for (auto entry = m_upperstates.begin(); entry != m_upperstates.end(); entry++)
		{
			fraction = (*entry)->BoltzmannExponent(T)/Z;
			(*entry)->SetFractionalOccupation(fraction);
		}
		m_current_temperature = T;
	}
	return true;
}


 size_t skSpectralLine_HitranEmission::g_numinstances = 0;

/*---------------------------------------------------------------------------
 *  skSpectralLine_HitranEmission::skSpectralLine_HitranEmission  2020-07-28 */
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLine_HitranEmission::skSpectralLine_HitranEmission( const skSpectralLine_HitranLine* hitranline, const HitranLine_UpperState* upperstate	)
	:skSpectralLine_HitranLine(*hitranline)
{
	g_numinstances++;
	m_upperstate         = upperstate;
	if (m_upperstate    != nullptr ) m_upperstate->AddRef();
	SetParentMolecule( hitranline->ParentMolecule());
}

/*---------------------------------------------------------------------------
 * skSpectralLine_HitranEmission::~skSpectralLine_HitranEmission  2020-07-28 */
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLine_HitranEmission::~skSpectralLine_HitranEmission()
{
	g_numinstances--;
	if (m_upperstate         != nullptr) m_upperstate->Release();
}

/*---------------------------------------------------------------------------
 *     skSpectralLine_HitranEmission::CalculateLineIntensity      2020-07-28 */
/** Calculkates the line intensity of the emission. Note that the temperature
 *	dependence of the upper state partition function is not calculated in this 
 *  function. It is called in the parent container in function, skSpectralEmissionCollection_HitranIsotope::SetAtmosphericState.
 *  This allows the parent container to do all of partition function for the
 *	upper states and just set the FractionalOccupation of each of the individual upper states.
 **/
/*---------------------------------------------------------------------------*/

bool skSpectralLine_HitranEmission::CalculateLineIntensity(  double T)
{
	double emission = EinsteinA()*m_upperstate->FractionalOccupation();
	SetCurrentLineIntensity( emission, T);
	return true;
}

/*---------------------------------------------------------------------------
 * Hitran_UpperState_SpectralLine_Entry::Hitran_UpperState_SpectralLine_Entry2020-07-27 */
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLineEntry_HitranEmission::skSpectralLineEntry_HitranEmission	( )
{
}

/*---------------------------------------------------------------------------
 * Hitran_UpperState_SpectralLine_Entry::Hitran_UpperState_SpectralLine_Entry2020-07-27 */
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLineEntry_HitranEmission::skSpectralLineEntry_HitranEmission	( skSpectralLine_HitranEmission* line, skSpectralLineShape* lineshapeobject 	)
{
	SetSpectralLine( line );
	SetLineShapeObject( lineshapeobject );
}

/*---------------------------------------------------------------------------
 * Hitran_UpperState_SpectralLine_Entry::Hitran_UpperState_SpectralLine_Entry2020-07-27 */
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLineEntry_HitranEmission::skSpectralLineEntry_HitranEmission	( const skSpectralLineEntry_HitranEmission& other )
	: skSpectralLineEntry( other )
{
}


/*---------------------------------------------------------------------------
 * Hitran_UpperState_SpectralLine_Entry::~Hitran_UpperState_SpectralLine_Entry2020-07-27 */
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLineEntry_HitranEmission::~skSpectralLineEntry_HitranEmission	( )
{
}

/*---------------------------------------------------------------------------
 *      HitranLine_SelectUpperStatesFromLines::median_value       2020-07-24 */
/** **/

/*---------------------------------------------------------------------------*/

double HitranLine_SelectUpperStatesFromLines::median_value( std::vector<double>& sorted_mu_array  )
{
	size_t n = sorted_mu_array.size()/2;

	return sorted_mu_array.at(n); 

}


/*---------------------------------------------------------------------------
 *    HitranLine_UpperStates::make_list_of_unique_upperstates     2020-07-24 */
/** This is a method that finds the list of unique upper state levels from the
 *	HITRAN spectroscopic data. The energy of the upper  state is found for 
 *	each line and sorted into ascending order. Upper states with the same energy 
 *	and the same degeneragcy are considered to be the same uppe state. The energy of the
 *	upper state is the "same" if two upper state energies are within a pre-set
 *	threshold given by m_delta_E, eg 0.0005 cm-1.
 *
 *	The algorithm, seemed to work on O2 A-Band upper states. However it is not 
 *	completely robust and it may be possible it gets quietly confused if the energy value
 *	of one upper state from one line is sorted into the middle of energy values
 *	of another upper state from other spectral lines.
 *
 *	The code is given a list of Hitran Lines read from the Hitran database. It calculates
 *	and returns all of the unique upper states represented by the list of lines. It also
 *	returns the original hitran lines but converted to an Hitran Emissions.
 *
 **/
/*---------------------------------------------------------------------------*/

bool HitranLine_SelectUpperStatesFromLines::make_list_of_unique_upperstates	(	std::list< const skSpectralLine_HitranLine*>&	all_spectral_lines,
																				HitranLine_UpperStates*							upperstates,
																				std::list<skSpectralLine_HitranEmission*>*		spectrallines)
{
	bool							ok  = true;
	bool							ok1;
	bool							samestate;
	double							lastdegeneracy = -1.0;
	double							lastE         = -9999.0;
	double							degeneracy;
	double							E;
	HitranLine_UpperState*			upperstate = nullptr;
	std::vector<double>				state_energy;

	state_energy.reserve(100);																							// Reserve some space to store all the lines associoated with current uppe level state
	all_spectral_lines.sort( 	lessthan_Eupper() );																	//  sort the Raw HITRAN lines into ascending wavenumber
	spectrallines->clear();																								// Clear out the list 

	for (auto& entry : all_spectral_lines)
	{
		E          = entry->EUpper();																					// Get the Energy of the upper state of this Hitran spectral line
		degeneracy = entry->HitranEntry().m_upperStatWt;																// Get the degeneracy of this Hitran line
		samestate = (upperstate != nullptr) && ( degeneracy == lastdegeneracy ) && ( ( E - lastE ) < m_delta_E);		// Is this still the same upper state as the last line we processed
		if (!samestate)																									// Nope. We have a new upper state
		{																												// So
			//printf("Finishing off last state\n");
			if (upperstate  != nullptr) upperstate->SetEnergyLevel( median_value( state_energy ) );						// If the last upper state is defined then set its energy to the median value of the upepr state energies 
			state_energy.resize(0);																						// Reset our list of upper state energies
			ok1 = upperstates->CreateNewUpperState( &upperstate, E, degeneracy);										// And create a new upper state object
			lastE = E;
			lastdegeneracy = degeneracy;
			ok = ok  && ok1;																					
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"HitranLine_SelectUpperStatesFromLines::make_list_of_unique_upperstates, Error creating a new upper state. Thats not good");
			}
		}
		//printf("%15.8f %5.1f\n",(double)E, (double)degeneracy);
		skSpectralLine_HitranEmission* spectralemission = new skSpectralLine_HitranEmission( entry, upperstate );		// CReate the skSpectralLine_HitranEmission object. Dont add a reference. 
		spectrallines->push_back( spectralemission );																	// And push it that is all we need for the moment, so just push it ontop the back of our list of lines
		state_energy.push_back( E );
	}
	if (upperstate  != nullptr) upperstate->SetEnergyLevel( median_value( state_energy ) );
	return ok;
}
