#pragma once
#include "sktran_disco/sktran_do.h"
#include "sktran_do_specs.h"

namespace sasktran_disco
{
	namespace testing
	{
		struct TestLayerSpecBase {
			virtual ~TestLayerSpecBase() {};
		};

		struct TestLayerSpecHG : public TestLayerSpecBase
		{
			double optical_depth;
			double ssa;
			double hg_asym;

			TestLayerSpecHG(double optical_depth, double ssa, double hg_asym) {
				this->optical_depth = optical_depth;
				this->ssa = ssa;
				this->hg_asym = hg_asym;
			}
		};

		struct TestLayerSpecRayleigh : public TestLayerSpecBase
		{
			double optical_depth;
			double ssa;
			double depol;

			TestLayerSpecRayleigh(double optical_depth, double ssa, double depol) {
				this->optical_depth = optical_depth;
				this->ssa = ssa;
				this->depol = depol;
			}
		};

		struct TestLayerSpecSiewert : public TestLayerSpecBase
		{
			double optical_depth;
			double ssa;

			TestLayerSpecSiewert(double optical_depth, double ssa) {
				this->optical_depth = optical_depth;
				this->ssa = ssa;
			}
		};

		template <int NSTOKES, int CNSTR=-1>
		class TestLayer
		{
		public:
			TestLayer(uint NSTR, const TestLayerSpecBase& spec) {
				if (dynamic_cast<const TestLayerSpecHG*>(&spec) != nullptr) {
					assignHGLayer(NSTR, dynamic_cast<const TestLayerSpecHG&>(spec));
				}
				else if (dynamic_cast<const TestLayerSpecRayleigh*>(&spec) != nullptr) {
					assignRayleighLayer(NSTR, dynamic_cast<const TestLayerSpecRayleigh&>(spec));
				}
				else if (dynamic_cast<const TestLayerSpecSiewert*>(&spec) != nullptr) {
					assignSiewertLayer(NSTR, dynamic_cast<const TestLayerSpecSiewert&>(spec));
				}
			}

			double optical_depth;
			double ssa;
			std::vector<LegendreCoefficient<NSTOKES>> lephasef;

		private:
			void assignHGLayer(uint NSTR, const TestLayerSpecHG& spec);
			void assignRayleighLayer(uint NSTR, const TestLayerSpecRayleigh& spec);
			void assignSiewertLayer(uint NSTR, const TestLayerSpecSiewert& spec);


		};

		struct TestTOAIntensity
		{
			double direct;
			double diffuse;
		};
		struct TestSolarSpec
		{
			double csz;
			double saz;
			TestTOAIntensity intensities;
		};

		template <int NSTOKES, int CNSTR=-1>
		class TestAtmosphere
		{
		public:
			TestAtmosphere(std::vector<TestLayerSpecHG>&& layers) {
				m_layershg.reserve(layers.size());
				for(TestLayerSpecHG& layer : layers) m_layershg.emplace_back(layer);
			}
			TestAtmosphere(std::vector<TestLayerSpecRayleigh>&& layers) {
				m_layersray.reserve(layers.size());
				for (TestLayerSpecRayleigh& layer : layers) m_layersray.emplace_back(layer);
			}
			TestAtmosphere(std::vector<TestLayerSpecSiewert>&& layers) {
				m_layersray.reserve(layers.size());
				for (TestLayerSpecSiewert& layer : layers) m_layerssiewert.emplace_back(layer);
			}

			std::vector<TestLayer<NSTOKES>> build(uint NSTR) {
				std::vector<TestLayer<NSTOKES>> layers; 
				layers.reserve(m_layershg.size() + m_layersray.size());
				for(const TestLayerSpecHG& spec : m_layershg)
					layers.push_back(TestLayer<NSTOKES>(NSTR, spec));
				for (const TestLayerSpecRayleigh& spec : m_layersray)
					layers.push_back(TestLayer<NSTOKES>(NSTR, spec));
				for (const TestLayerSpecSiewert& spec : m_layerssiewert)
					layers.push_back(TestLayer<NSTOKES>(NSTR, spec));
				return layers;
			}
		private:
			std::vector<TestLayerSpecHG> m_layershg;
			std::vector<TestLayerSpecRayleigh> m_layersray;
			std::vector<TestLayerSpecSiewert> m_layerssiewert;

		};


		struct TestLOS
		{
			double coszen;
			double az;
		};
		
		template <int NSTOKES, int CNSTR=-1>
		struct TestCase
		{
			TestCase(uint NSTR, 
					 TestSolarSpec& solar_spec, 
					 TestAtmosphere<NSTOKES> atmo, 
					 double ALBEDO, 
					 const std::vector<TestLOS>& los, 
					 const std::vector<double>& correct)
			{
				linesofsight = los;
				nstr = NSTR;
				solar = solar_spec;
				lambertian = ALBEDO;
				is_lambertian = true;
				layers = atmo.build(NSTR);
				nlyr = (uint) layers.size();
				correct_radiances = &correct;
			}
			TestCase(uint NSTR,
					 TestSolarSpec& solar_spec,
					 TestAtmosphere<NSTOKES>& atmo,
					 std::function<double(double, double, double)> the_brdf,
					 const std::vector<TestLOS>& los,
					 const std::vector<double>& correct)
			{
				linesofsight = los;
				nstr = NSTR;
				solar = solar_spec;
				brdf = the_brdf;
				is_lambertian = false;
				layers = atmo.build(NSTR);
				nlyr = (uint) layers.size();
				correct_radiances = &correct;
			}

			uint nstr;
			uint nlyr;
			TestSolarSpec solar;


			double lambertian;
			std::function<double(double, double, double)> brdf;
			bool is_lambertian;

			std::unique_ptr<TestBRDF> getBRDF() const {
				if(is_lambertian) return std::unique_ptr<TestBRDF>(new TestBRDF(lambertian));
				else return std::unique_ptr<TestBRDF>(new TestBRDF(brdf));
			}
			
			std::vector<TestLayer<NSTOKES>> layers;
			std::vector<TestLOS> linesofsight;
			const std::vector<double>* correct_radiances;
		};
	}

	template <int NSTOKES, int CNSTR=-1>
	class SKTRAN_DO_TestSpec: public SKTRAN_DO_UserSpec
	{
	public:

		struct TestWF: WeightingFunctionSpec
		{
			TestWF() = default;
			TestWF(double a_ssa, double a_optd, double alb, uint layer) {
				ssa = a_ssa;
				optd = a_optd;
				albedo = alb;
				layer_index = layer;
			}
			TestWF(double a_ssa, double a_optd, double alb, uint layer, const std::vector<double>& legendre) {
				ssa = a_ssa;
				optd = a_optd;
				albedo = alb;
				layer_index = layer;
				legendrecoeff = legendre;
			}

            WeightingFunctionType type() override { return WeightingFunctionSpec::TestWF; }


            double ssa;
			double optd;
			double albedo;
			LayerIndex layer_index;
			std::vector<double> legendrecoeff;
		};

		void configure(const testing::TestCase<NSTOKES>& testcase)
		{
			setNumberOfStreams(testcase.nstr);
			cacheLPOfStreamAngles();
			setNumberOfLayers(testcase.nlyr);
			setTOAIntensities(testcase.solar.intensities.direct);
			setWFReturnForm(WeightingFunctionForm::dI_dLogX);
			m_testcase = &testcase;
		}
		const testing::TestCase<NSTOKES>* getTestCase() const {
			return m_testcase;
		}

	private:
		const testing::TestCase<NSTOKES>* m_testcase;
	};
}
