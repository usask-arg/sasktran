#pragma once
#include "modules/sktran_do_deprecated/include/sktran_do.h"

namespace sktran_do_detail
{
	// BRDF interface SKDO 
	class BRDF_Base
	{
	public:
		// BRDF arguments:
		// coszen_out : cosine of outgoing ray zenith angle
		// coszen_in : cosine of incoming ray zenith angle
		// az_diff : difference in azimuth angles [units: radians] 
		// return: BRDF value, no pi normalization, integral over BRDF should -> 1
		virtual double brdf(double coszen_out, double coszen_in, double az_diff) const = 0;
		virtual bool isLambertian() const = 0;
	};

	// Support for the BRDF_Base interface for skBRDF objects.
	class Wrapped_skBRDF : public BRDF_Base
	{
	public:
		Wrapped_skBRDF(double wavlen, const skBRDF& skbrdf, GEODETIC_INSTANT ref_inst):
			m_brdf(skbrdf),
			m_wavlen(wavlen),
			m_reference_instant(ref_inst),
			m_is_lambertian(skbrdf.IsLambertian())
		{
			// empty
		}

		virtual double brdf(double coszen_out, double coszen_in, double az_diff) const override {
			double reflected;
			m_brdf.BRDF(m_wavlen, m_reference_instant, abs(coszen_in), abs(coszen_out), cos(PI - az_diff), &reflected);

			// Imporant note: SKTRAN normalized BRDF w/ the div by Pi. SKDO doesn't => we multiply by Pi.
			return PI * reflected;
		}

		virtual bool isLambertian() const override {
			return m_is_lambertian;
		}

	private:
		const skBRDF& m_brdf;
		const double m_wavlen;
		GEODETIC_INSTANT m_reference_instant;
		const bool m_is_lambertian;
	};

	// Allows tests to be setup. Overrides standard configuration so writing 
	// tests is easy.
	class TestBRDF: public BRDF_Base
	{
	public:
		TestBRDF() {}
		TestBRDF(std::function<double(double, double, double)> brdf) {
			setBRDF(brdf, false);
		}
		TestBRDF(double lambertian) {
			setBRDF([=](double, double, double) { return lambertian; }, true);
		}
		void setBRDF(std::function<double(double, double, double)> brdf, bool lambertian) {
			m_brdf = brdf;
			m_is_lambertian = lambertian;
		}

		virtual double brdf(double coszen_out, double coszen_in, double az_diff) const override {
			return m_brdf(coszen_out, coszen_in, az_diff);
		}

		virtual bool isLambertian() const override {
			return m_is_lambertian;
		}

	private:
		std::function<double(double, double, double)> m_brdf;
		bool m_is_lambertian;

	};

	template<class NonLambertianBRDF>
	class PretendNotLambertian : public NonLambertianBRDF
	{
	public:
		using NonLambertianBRDF::NonLambertianBRDF;
		virtual bool isLambertian() const override {
			return false;
		}
	};

	// Surface object used internally by SASKTRAN-Disco. Handles the expansion
	// of the BRDF function, and caches all reflection angles which will be 
	// requested internally by SKDO.
	class Albedo
	{
	public:
		Albedo() {}

		void configure(AEOrder m, const std::vector<LineOfSight>& los, const std::vector<double>& streams,
					   double csz, 
					   BRDF_Base* brdf,
					   uint nterms);

		bool isLambertian() const {
			return m_brdf->isLambertian();
		}

		inline const std::vector<double>& losBDRFromStreams(uint los_idx) const {
			return m_brdf_los_stream[los_idx];
		}

		inline const std::vector<double>& streamBDRFromStreams(uint out_stream_idx) const {
			return m_brdf_stream_stream[out_stream_idx];
		}

		inline double losBDRFromSun(uint los_idx) const {
			return m_brdf_los_solar[los_idx];
		}
		inline double streamBDRFromSun(StreamIndex str_idx) const {
			return m_brdf_stream_solar[str_idx];
		}

		inline double d_streamBRDFFromSun(AEOrder m, StreamIndex str_idx) const {
			if (isLambertian()) {
				return kronDelta(m, 0);
			}
			else {
				return m_brdf_stream_solar[str_idx];
			}
		}

		
	private:

		double computeBDR(AEOrder m, double outgoing, double incoming) const;

	private: // members
		VectorDim2<double> m_brdf_los_stream;
		VectorDim2<double> m_brdf_stream_stream;
		VectorDim1<double> m_brdf_los_solar;
		VectorDim1<double> m_brdf_stream_solar;

		bool m_configured;
		const double* m_gq_angles;
		const double* m_gq_weights;
		uint m_nterms;

		BRDF_Base* m_brdf;
	};

	class SurfaceEmission {
	public:
		SurfaceEmission(const std::vector<double>& wavelengths,
			const std::vector<double>& emissions) :
		m_wavelengths(wavelengths),
		m_emissions(emissions)
		{};

		double emission(double wavelength) const {
			// Interpolate in wavelength, set to 0 outside
			return nxLinearInterpolate::EvaluateYatX(wavelength, m_wavelengths, m_emissions, nxLinearInterpolate::EnumOutOfBoundAction::ENUM_MISSINGVALUE, 0.0);
		}

	private:
		std::vector<double> m_wavelengths;
		std::vector<double> m_emissions;
	};
}