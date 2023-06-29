#pragma once

// Include Boost's random number generator
#include <ctime>            // std::time
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
// Sun CC doesn't handle std::iterator_adaptor yet
#if !defined(__SUNPRO_CC) || (__SUNPRO_CC > 0x530)
#include <boost/generator_iterator.hpp>
#endif

#ifdef BOOST_NO_STDC_NAMESPACE
namespace std {
  using ::time;
}
#endif

typedef boost::mt11213b SKTRAN_RNG_generator_type;
typedef boost::uniform_real<> SKTRAN_RNG_distribution_type;
typedef boost::variate_generator<SKTRAN_RNG_generator_type&, SKTRAN_RNG_distribution_type > SKTRAN_rng;	// A random number generator for u(0,1)

static SKTRAN_RNG_distribution_type SKTRAN_RNG_globalDistribution(0,1);

class SKTRAN_RNG 
{
	private:
		//typedef std::mt11213b base_generator_type;
		//typedef std::variate_generator<base_generator_type&, std::uniform_real<> > SKTRAN_rng;	// A random number generator for u(0,1)

		mutable SKTRAN_rng m_randGen;

	public:
		 SKTRAN_RNG();
		~SKTRAN_RNG();
        //SKTRAN_RNG( const SKTRAN_RNG& ) = delete; // May want to enable copy constructor to change how multithreading is handled

		double operator()() const;                     // Get a random number
		void   SetSeed( const std::uint32_t seed );  // Set generator seed value
}; 


// Below is just commented-out test code from Boost

//#include <iostream>
//#include <fstream>
//#include <ctime>            // std::time
//
//#include <boost/random/linear_congruential.hpp>
//#include <boost/random/uniform_int.hpp>
//#include <boost/random/uniform_real.hpp>
//#include <boost/random/variate_generator.hpp>
//#include "boost/random.hpp"
//#include "boost/generator_iterator.hpp"
//
//// Sun CC doesn't handle std::iterator_adaptor yet
//#if !defined(__SUNPRO_CC) || (__SUNPRO_CC > 0x530)
//#include <boost/generator_iterator.hpp>
//#endif
//
//#ifdef BOOST_NO_STDC_NAMESPACE
//namespace std {
//  using ::time;
//}
//#endif
//
//// This is a typedef for a random number generator.
//// Try std::mt19937 or std::ecuyer1988 instead of std::minstd_rand
//typedef std::mt19937 base_generator_type;
//
//// This is a reproducible simulation experiment.  See main().
//void experiment(base_generator_type & generator)
//{
//  // Define a uniform random number distribution of integer values between
//  // 1 and 6 inclusive.
//  typedef std::uniform_int<> distribution_type;
//  typedef std::variate_generator<base_generator_type&, distribution_type> gen_type;
//  gen_type die_gen(generator, distribution_type(1, 6));
//
//#if !defined(__SUNPRO_CC) || (__SUNPRO_CC > 0x530)
//  // If you want to use an STL iterator interface, use iterator_adaptors.hpp.
//  // Unfortunately, this doesn't work on SunCC yet.
//  std::generator_iterator<gen_type> die(&die_gen);
//  for(int i = 0; i < 10; i++)
//    std::cout << *die++ << " ";
//  std::cout << '\n';
//#endif
//}
//
//double randomNumber()
//{
//  // Define a random number generator and initialize it with a reproducible
//  // seed.
//  // (The seed is unsigned, otherwise the wrong overload may be selected
//  // when using mt19937 as the base_generator_type.)
//	base_generator_type generator((const std::uint32_t)std::time(0));
//
//  // std::cout << "10 samples of a uniform distribution in [0..1):\n";
//
//  // Define a uniform random number distribution which produces "double"
//  // values between 0 and 1 (0 inclusive, 1 exclusive).
//  std::uniform_real<> uni_dist(0,1);
//  std::variate_generator<base_generator_type&, std::uniform_real<> > uni(generator, uni_dist);
//
//  std::cout.setf(std::ios::fixed);
//  // You can now retrieve random numbers from that distribution by means
//  // of a STL Generator interface, i.e. calling the generator as a zero-
//  // argument function.
//  // for(int i = 0; i < 10; i++)
//  //  std::cout << uni() << '\n';
//
//  return uni();
//  /*
//   * Change seed to something else.
//   *
//   * Caveat: std::time(0) is not a very good truly-random seed.  When
//   * called in rapid succession, it could return the same values, and
//   * thus the same random number sequences could ensue.  If not the same
//   * values are returned, the values differ only slightly in the
//   * lowest bits.  A linear congruential generator with a small factor
//   * wrapped in a uniform_smallint (see experiment) will produce the same
//   * values for the first few iterations.   This is because uniform_smallint
//   * takes only the highest bits of the generator, and the generator itself
//   * needs a few iterations to spread the initial entropy from the lowest bits
//   * to the whole state.
//   */
//  /*
//
//  generator.seed(static_cast<unsigned int>(std::time(0)));
//
//  std::cout << "\nexperiment: roll a die 10 times:\n";
//
//  // You can save a generator's state by copy construction.
//  base_generator_type saved_generator = generator;
//
//  // When calling other functions which take a generator or distribution
//  // as a parameter, make sure to always call by reference (or pointer).
//  // Calling by value invokes the copy constructor, which means that the
//  // sequence of random numbers at the caller is disconnected from the
//  // sequence at the callee.
//  experiment(generator);
//
//  std::cout << "redo the experiment to verify it:\n";
//  experiment(saved_generator);
//
//  // After that, both generators are equivalent
//  assert(generator == saved_generator);
//
//  // as a degenerate case, you can set min = max for uniform_int
//  std::uniform_int<> degen_dist(4,4);
//  std::variate_generator<base_generator_type&, std::uniform_int<> > deg(generator, degen_dist);
//  std::cout << deg() << " " << deg() << " " << deg() << std::endl;
//  
//#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
//  {
//    // You can save the generator state for future use.  You can read the
//    // state back in at any later time using operator>>.
//    std::ofstream file("rng.saved", std::ofstream::trunc);
//    file << generator;
//  }
//#endif
//  // Some compilers don't pay attention to std:3.6.1/5 and issue a
//  // warning here if "return 0;" is omitted.
//  return 0.0;
//
//  */
//}
