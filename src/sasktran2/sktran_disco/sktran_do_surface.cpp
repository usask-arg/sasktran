#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_surface.h"

void  sasktran_disco::Albedo::configure(AEOrder m,
                                        const std::vector<LineOfSight>& los,
                                        const std::vector<double>& streams,
                                        double csz,
                                        BRDF_Base* brdf,
                                        uint nterms)
{
	// Resize vector members and set invalid locations nan
	m_brdf = brdf;

	// get quadrature specs
	m_nterms = nterms;
	m_gq_angles = getQuadratureAbscissae(m_nterms);
	m_gq_weights = getQuadratureWeights(m_nterms);

	// cache streams from streams
	const uint NSTR = static_cast<const uint>(streams.size());

	m_brdf_los_stream.resize(los.size(), std::vector<double>(NSTR));
	m_brdf_stream_stream.resize(NSTR/2, std::vector<double>(NSTR));
	m_brdf_los_solar.resize(los.size());
	m_brdf_stream_solar.resize(NSTR/2);

	// resize

	for(auto& a_los : m_brdf_los_stream) {
		std::fill_n(a_los.begin(), NSTR / 2, std::nan("1")); // upwelling albedo source (directly) is not possible so set these cases to nan
	}
	for(auto& an_up_stream : m_brdf_stream_stream) {
		std::fill_n(an_up_stream.begin(), NSTR / 2, std::nan("1")); // same as above
	}

	// compute 
	for(const auto& l : los) {
		double out = l.coszenith;

		// compute from solar source
		m_brdf_los_solar[l.unsorted_index] = computeBDR(m, out, -csz);

		auto& curr_out_los = m_brdf_los_stream[l.unsorted_index];
		for(uint i = NSTR / 2; i < NSTR; ++i) {
			double in = streams[i];
			curr_out_los[i] = computeBDR(m, out, in);
		}
	}

	for(uint str_idx = 0; str_idx < NSTR / 2; ++str_idx) {
		double out = streams[str_idx];
		m_brdf_stream_solar[str_idx] = computeBDR(m, out, -csz);

		auto& curr_out_stream = m_brdf_stream_stream[str_idx];
		for(uint i = NSTR / 2; i < NSTR; ++i) {
			double in = streams[i];
			curr_out_stream[i] = computeBDR(m, out, in);
		}
	}

}


double sasktran_disco::Albedo::computeBDR(AEOrder m, double outgoing, double incoming) const {
	if(isLambertian()) {
		return (m == 0) ? m_brdf->brdf(outgoing, incoming, PI) : 0;
	}

	double integral = 0;
	for(uint i = 0; i < m_nterms / 2; ++i) { 
		// inplace double gauss quadrature rule
		double a1 = 0.5 * m_gq_angles[i] + 0.5;
		double a2 = -0.5 * m_gq_angles[i] + 0.5;
		double a3 = 0.5 * m_gq_angles[i] - 0.5;
		double a4 = -0.5 * m_gq_angles[i] - 0.5;
		double w = 0.5 * m_gq_weights[i];

		integral += w * m_brdf->brdf(outgoing, incoming, PI * a1) * cos(m * PI * a1);
		integral += w * m_brdf->brdf(outgoing, incoming, PI * a2) * cos(m * PI * a2);
		integral += w * m_brdf->brdf(outgoing, incoming, PI * a3) * cos(m * PI * a3);
		integral += w * m_brdf->brdf(outgoing, incoming, PI * a4) * cos(m * PI * a4);
	}
	return 0.5 * (2.0 - kronDelta(m, 0)) * integral;
}