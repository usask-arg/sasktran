#pragma once
#include "sktran_disco/sktran_do.h"

namespace sasktran_disco
{
	// Get quadrature angles and weights.
	// Note : SASKTRAN-DO uses the double quadrature scheme common to most
	// discrete-ordinate algorithms. Essentially integration over [-1, 1] is 
	// handled as separate integrals over [-1, 0] and [0, 1].
	void getStreamsAndWeights(uint num_streams, 
							  VectorDim1<double>& angles, 
							  VectorDim1<double>& weights);
	void getStreamsAndWeights(const std::vector<double>& given_angle,
							  const std::vector<double>& given_weights,
							  VectorDim1<double>& angles,
							  VectorDim1<double>& weights);

	// Get absolute value of high-precision quadrature angles. Supported sizes
	// are: 64, 128, 256, 512, 1024. Size of returned is nterms/2.
	const double* getQuadratureAbscissae(uint nterms);

	// Get the weights corresponding to positive quadrature angles. Size of 
	// returned is nterms/2.
	const double* getQuadratureWeights(uint nterms);

	namespace tables {
		extern const std::array<double, 109> quadrature_angles;
		extern const std::array<double, 109> quadrature_weights;
		extern const std::array<double, 32 > gq64_nodes;
		extern const std::array<double, 32 > gq64_weights;
		extern const std::array<double, 64 > gq128_nodes;
		extern const std::array<double, 64 > gq128_weights;
		extern const std::array<double, 128 > gq256_nodes;
		extern const std::array<double, 128 > gq256_weights;
		extern const std::array<double, 256 > gq512_nodes;
		extern const std::array<double, 256 > gq512_weights;
		extern const std::array<double, 512 > gq1024_nodes;
		extern const std::array<double, 512 > gq1024_weights;
		extern const std::map<uint, const double*> gqnodes;
		extern const std::map<uint, const double*> gqweights;
	}

	double gq_integral(double start, double end, std::function<double(double)> func);
}
