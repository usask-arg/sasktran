#include "include/sktran_montecarlo_internals.h"

SKTRAN_RNG::SKTRAN_RNG() : m_randGen(*(new SKTRAN_RNG_generator_type), SKTRAN_RNG_globalDistribution)
{
}

SKTRAN_RNG::~SKTRAN_RNG()
{
    SKTRAN_RNG_generator_type* createdInInitializer = &m_randGen.engine();
	delete createdInInitializer;
}

double SKTRAN_RNG::operator()() const
{
	return m_randGen();
}

void SKTRAN_RNG:: SetSeed( const std::uint32_t seed )
{
	m_randGen.engine().seed(seed);
}
