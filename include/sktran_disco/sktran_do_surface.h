#pragma once
#include "sktran_disco/sktran_do.h"

namespace sasktran_disco
{
	// BRDF interface SKDO 
	class BRDF_Base
	{
	public:
        virtual ~BRDF_Base() {};

        // BRDF arguments:
		// coszen_out : cosine of outgoing ray zenith angle
		// coszen_in : cosine of incoming ray zenith angle
		// az_diff : difference in azimuth angles [units: radians] 
		// return: BRDF value, no pi normalization, integral over BRDF should -> 1
		virtual double brdf(double coszen_out, double coszen_in, double az_diff) const = 0;
		virtual bool isLambertian() const = 0;
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
            // TODO: Fix this
            return 0.0;
			//return nxLinearInterpolate::EvaluateYatX(wavelength, m_wavelengths, m_emissions, nxLinearInterpolate::EnumOutOfBoundAction::ENUM_MISSINGVALUE, 0.0);
		}

	private:
		std::vector<double> m_wavelengths;
		std::vector<double> m_emissions;
	};
}