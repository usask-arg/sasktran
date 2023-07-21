#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_testing.h"

template <int NSTOKES, int CNSTR>
void sasktran_disco::testing::TestLayer<NSTOKES, CNSTR>::assignHGLayer(uint NSTR, const TestLayerSpecHG& spec) {
	optical_depth = spec.optical_depth;
	ssa = spec.ssa;
	lephasef.resize(NSTR);
	for (LPOrder l = 0; l < NSTR; ++l) {
		lephasef[l].a1 = (2 * l + 1) * pow(spec.hg_asym, l);
	}
}

template <>
void sasktran_disco::testing::TestLayer<1>::assignRayleighLayer(uint NSTR, const TestLayerSpecRayleigh& spec) {
	lephasef.resize(NSTR);

	double d = spec.depol;

	optical_depth = spec.optical_depth;
	ssa = spec.ssa;

	lephasef[0].a1 = 1;
	lephasef[2].a1 = (1 - d) / (2 + d);
}

template <>
void sasktran_disco::testing::TestLayer<3>::assignRayleighLayer(uint NSTR, const TestLayerSpecRayleigh& spec) {
	lephasef.resize(NSTR);

	double d = spec.depol;

	optical_depth = spec.optical_depth;
	ssa = spec.ssa;

	lephasef[0].a1 = 1;
	lephasef[2].a2 = 6 * (1 - d) / (2 + d);
	lephasef[2].a1 = (1 - d) / (2 + d);
	lephasef[2].b1 = sqrt(6.0) * (1 - d) / (2 + d);
}

template <>
void sasktran_disco::testing::TestLayer<4>::assignRayleighLayer(uint NSTR, const TestLayerSpecRayleigh& spec) {
    lephasef.resize(NSTR);

    double d = spec.depol;

    optical_depth = spec.optical_depth;
    ssa = spec.ssa;

    lephasef[0].a1 = 1;
    lephasef[2].a2 = 6 * (1 - d) / (2 + d);
    lephasef[2].a1 = (1 - d) / (2 + d);
    lephasef[2].b1 = sqrt(6.0) * (1 - d) / (2 + d);

    lephasef[1].a4 = 3 * (1-2*d) / (2+d);
}


template <>
void sasktran_disco::testing::TestLayer<1>::assignSiewertLayer(uint NSTR, const TestLayerSpecSiewert& spec) {
}

template <>
void sasktran_disco::testing::TestLayer<3>::assignSiewertLayer(uint NSTR, const TestLayerSpecSiewert& spec) {
    lephasef.resize(NSTR);

    std::vector<LegendreCoefficient<3>> temp;
    temp.resize(12);

    optical_depth = spec.optical_depth;
    ssa = spec.ssa;

    temp[0].a1 = 1.0;
    temp[1].a1 = 2.104031;
    temp[2].a1 = 2.095158;
    temp[3].a1 = 1.414939;
    temp[4].a1 = 0.703593;
    temp[5].a1 = 0.235001;
    temp[6].a1 = 0.064039;
    temp[7].a1 = 0.012837;
    temp[8].a1 = 0.002010;
    temp[9].a1 = 0.000246;
    temp[10].a1 = 0.000024;
    temp[11].a1 = 0.000002;

    temp[0].a2 = 0.0;
    temp[1].a2 = 0.0;
    temp[2].a2 = 3.726079;
    temp[3].a2 = 2.202868;
    temp[4].a2 = 1.190694;
    temp[5].a2 = 0.391203;
    temp[6].a2 = 0.105556;
    temp[7].a2 = 0.020484;
    temp[8].a2 = 0.003097;
    temp[9].a2 = 0.000366;
    temp[10].a2 = 0.000035;
    temp[11].a2 = 0.000003;

    temp[0].b1 = 0.0;
    temp[1].b1 = 0.0;
    temp[2].b1 = -0.116688;
    temp[3].b1 = -0.209370;
    temp[4].b1 = -0.227137;
    temp[5].b1 = -0.144524;
    temp[6].b1 = -0.052640;
    temp[7].b1 = -0.012400;
    temp[8].b1 = -0.002093;
    temp[9].b1 = -0.000267;
    temp[10].b1 = -0.000027;
    temp[11].b1 = -0.000002;

    temp[0].a3 = 0.0;
    temp[1].a3 = 0.0;
    temp[2].a3 = 3.615946;
    temp[3].a3 = 2.240516;
    temp[4].a3 = 1.139473;
    temp[5].a3 = 0.365605;
    temp[6].a3 = 0.082779;
    temp[7].a3 = 0.013649;
    temp[8].a3 = 0.001721;
    temp[9].a3 = 0.000172;
    temp[10].a3 = 0.000014;
    temp[11].a3 = 0.000001;

    for (int i = 0; i < std::min(int(NSTR), int(12)); ++i) {
        lephasef[i] = temp[i];
    }
}

template <>
void sasktran_disco::testing::TestLayer<4>::assignSiewertLayer(uint NSTR, const TestLayerSpecSiewert& spec) {
	lephasef.resize(NSTR);

	std::vector<LegendreCoefficient<4>> temp;
	temp.resize(12);

	optical_depth = spec.optical_depth;
	ssa = spec.ssa;

	temp[0].a1 = 1.0;
	temp[1].a1 = 2.104031;
	temp[2].a1 = 2.095158;
	temp[3].a1 = 1.414939;
	temp[4].a1 = 0.703593;
	temp[5].a1 = 0.235001;
	temp[6].a1 = 0.064039;
	temp[7].a1 = 0.012837;
	temp[8].a1 = 0.002010;
	temp[9].a1 = 0.000246;
	temp[10].a1 = 0.000024;
	temp[11].a1 = 0.000002;

	temp[0].a2 = 0.0;
	temp[1].a2 = 0.0;
	temp[2].a2 = 3.726079;
	temp[3].a2 = 2.202868;
	temp[4].a2 = 1.190694;
	temp[5].a2 = 0.391203;
	temp[6].a2 = 0.105556;
	temp[7].a2 = 0.020484;
	temp[8].a2 = 0.003097;
	temp[9].a2 = 0.000366;
	temp[10].a2 = 0.000035;
	temp[11].a2 = 0.000003;

	temp[0].a4 = 0.915207;
	temp[1].a4 = 2.095727;
	temp[2].a4 = 2.008624;
	temp[3].a4 = 1.436545;
	temp[4].a4 = 0.706244;
	temp[5].a4 = 0.238475;
	temp[6].a4 = 0.056448;
	temp[7].a4 = 0.009703;
	temp[8].a4 = 0.001267;
	temp[9].a4 = 0.000130;
	temp[10].a4 = 0.000011;
	temp[11].a4 = 0.000001;

	temp[0].b1 = 0.0;
	temp[1].b1 = 0.0;
	temp[2].b1 = -0.116688;
	temp[3].b1 = -0.209370;
	temp[4].b1 = -0.227137;
	temp[5].b1 = -0.144524;
	temp[6].b1 = -0.052640;
	temp[7].b1 = -0.012400;
	temp[8].b1 = -0.002093;
	temp[9].b1 = -0.000267;
	temp[10].b1 = -0.000027;
	temp[11].b1 = -0.000002;

	temp[0].a3 = 0.0;
	temp[1].a3 = 0.0;
	temp[2].a3 = 3.615946;
	temp[3].a3 = 2.240516;
	temp[4].a3 = 1.139473;
	temp[5].a3 = 0.365605;
	temp[6].a3 = 0.082779;
	temp[7].a3 = 0.013649;
	temp[8].a3 = 0.001721;
	temp[9].a3 = 0.000172;
	temp[10].a3 = 0.000014;
	temp[11].a3 = 0.000001;
	
	for (int i = 0; i < std::min(int(NSTR), int(12)); ++i) {
		lephasef[i] = temp[i];
	}

}

template class sasktran_disco::testing::TestLayer<1>;
template class sasktran_disco::testing::TestLayer<3>;
template class sasktran_disco::testing::TestLayer<4>;
