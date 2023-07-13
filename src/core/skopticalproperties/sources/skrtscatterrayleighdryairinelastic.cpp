#include <skopticalproperties21.h>
#include <skopticalproperties21.h>

static skOpticalProperties_RayleighDryAir_Inelastic::sk_raman_partition_data n2partitiondata[] = {
	{ 0,  0, 6,    0.0000},
	{ 1,  1, 3,    3.9791},
	{ 2,  2, 6,   11.9373},
	{ 3,  3, 3,   23.8741},
	{ 4,  4, 6,   39.7892},
	{ 5,  5, 3,   59.6821},
	{ 6,  6, 6,   83.5521},
	{ 7,  7, 3,  111.3983},
	{ 8,  8, 6,  143.2197},
	{ 9,  9, 3,  179.0154},
	{10, 10, 6,  218.7839},
	{11, 11, 3,  262.5240},
	{12, 12, 6,  310.2341},
	{13, 13, 3,  361.9126},
	{14, 14, 6,  417.5576},
	{15, 15, 3,  477.1673},
	{16, 16, 6,  540.7395},
	{17, 17, 3,  608.2722},
	{18, 18, 6,  679.7628},
	{19, 19, 3,  755.2090},
	{20, 20, 6,  834.6081},
	{21, 21, 3,  917.9574},
	{22, 22, 6, 1005.2540},
	{23, 23, 3, 1096.4948},
	{24, 24, 6, 1191.6766},
	{25, 25, 3, 1290.7963},
	{26, 26, 6, 1393.8503},
	{27, 27, 3, 1500.8350},
	{28, 28, 6, 1611.7467}
};

static skOpticalProperties_RayleighDryAir_Inelastic::sk_raman_xs_data n2linedata[] = {
	{-194.3015, 25, 3, 1290.7963, 0.3601},
	{-186.4226, 24, 6, 1191.6766, 0.3595},
	{-178.5374, 23, 3, 1096.4948, 0.3589},
	{-170.6459, 22, 6, 1005.2540, 0.3581},
	{-162.7484, 21, 3,  917.9574, 0.3573},
	{-154.8453, 20, 6,  834.6081, 0.3565},
	{-146.9368, 19, 3,  755.2090, 0.3555},
	{-139.0233, 18, 6,  679.7628, 0.3544},
	{-131.1049, 17, 3,  608.2722, 0.3532},
	{-123.1819, 16, 6,  540.7395, 0.3519},
	{-115.2547, 15, 3,  477.1673, 0.3504},
	{-107.3235, 14, 6,  417.5576, 0.3487},
	{ -99.3886, 13, 3,  361.9126, 0.3467},
	{ -91.4502, 12, 6,  310.2341, 0.3443},
	{ -83.5086, 11, 3,  262.5240, 0.3416},
	{ -75.5642, 10, 6,  218.7839, 0.3383},
	{ -67.6171,  9, 3,  179.0154, 0.3344},
	{ -59.6676,  8, 6,  143.2197, 0.3294},
	{ -51.7162,  7, 3,  111.3983, 0.3231},
	{ -43.7629,  6, 6,   83.5521, 0.3147},
	{ -35.8080,  5, 3,   59.6821, 0.3030},
	{ -27.8519,  4, 6,   39.7892, 0.2857},
	{ -19.8950,  3, 3,   23.8741, 0.2571},
	{ -11.9373,  2, 6,   11.9373, 0.2000},
	{  11.9373,  0, 6,    0.0000, 1.0000},
	{  19.8950,  1, 3,    3.9791, 0.6000},
	{  27.8519,  2, 6,   11.9373, 0.5143},
	{  35.8080,  3, 3,   23.8741, 0.4762},
	{  43.7629,  4, 6,   39.7892, 0.4545},
	{  51.7162,  5, 3,   59.6821, 0.4406},
	{  59.6676,  6, 6,   83.5521, 0.4308},
	{  67.6171,  7, 3,  111.3983, 0.4235},
	{  75.5642,  8, 6,  143.2197, 0.4180},
	{  83.5086,  9, 3,  179.0154, 0.4135},
	{  91.4502, 10, 6,  218.7839, 0.4099},
	{  99.3886, 11, 3,  262.5240, 0.4070},
	{ 107.3235, 12, 6,  310.2341, 0.4044},
	{ 115.2547, 13, 3,  361.9126, 0.4023},
	{ 123.1819, 14, 6,  417.5576, 0.4004},
	{ 131.1049, 15, 3,  477.1673, 0.3988},
	{ 139.0233, 16, 6,  540.7395, 0.3974},
	{ 146.9368, 17, 3,  608.2722, 0.3961},
	{ 154.8453, 18, 6,  679.7628, 0.3950},
	{ 162.7484, 19, 3,  755.2090, 0.3940},
	{ 170.6459, 20, 6,  834.6081, 0.3931},
	{ 178.5374, 21, 3,  917.9574, 0.3922},
	{ 186.4226, 22, 6, 1005.2540, 0.3915},
	{ 194.3015, 23, 3, 1096.4948, 0.3908}

};

static skOpticalProperties_RayleighDryAir_Inelastic::sk_raman_partition_data o2partitiondata[] = {
	{ 1,  0, 1,    0.0000},
	{ 1,  2, 1,    2.0843},
	{ 1,  1, 1,    3.9611},
	{ 3,  2, 1,   16.2529},
	{ 3,  4, 1,   16.3876},
	{ 3,  3, 1,   18.3372},
	{ 5,  5, 1,   44.2117},
	{ 5,  6, 1,   42.2240},
	{ 5,  4, 1,   42.2001},
	{ 7,  7, 1,   81.5805},
	{ 7,  6, 1,   79.6070},
	{ 7,  8, 1,   79.5646},
	{ 9,  8, 1,  128.4921},
	{ 9, 10, 1,  128.3978},
	{ 9,  9, 1,  130.4376},
	{11, 11, 1,  190.7749},
	{11, 12, 1,  188.7135},
	{11, 10, 1,  188.8532},
	{13, 14, 1,  260.5011},
	{13, 12, 1,  260.6826},
	{13, 13, 1,  262.5829},
	{15, 14, 1,  343.9697},
	{15, 15, 1,  345.8500},
	{15, 16, 1,  343.7484},
	{17, 16, 1,  438.7015},
	{17, 17, 1,  440.5620},
	{17, 18, 1,  438.4418},
	{19, 19, 1,  546.7050},
	{19, 18, 1,  544.8628},
	{19, 20, 1,  544.5658},
	{21, 20, 1,  662.4368},
	{21, 22, 1,  662.1030},
	{21, 21, 1,  664.2610},
	{23, 23, 1,  793.2100},
	{23, 24, 1,  791.0344},
	{23, 22, 1,  791.4045},
	{25, 24, 1,  931.7450},
	{25, 25, 1,  933.5330},
	{25, 26, 1,  931.3390},
	{27, 28, 1, 1082.9941},
	{27, 27, 1, 1085.2060},
	{27, 26, 1, 1083.4356},
	{29, 29, 1, 1248.2040},
	{29, 28, 1, 1246.4518},
	{29, 30, 1, 1245.9750},
	{31, 31, 1, 1422.5020},
	{31, 30, 1, 1420.7672},
	{31, 32, 1, 1420.2552},
	{33, 34, 1, 1605.8064},
	{33, 33, 1, 1608.0710},
	{33, 32, 1, 1606.3533},
	{35, 36, 1, 1802.5983},
	{35, 34, 1, 1803.1802},
	{35, 35, 1, 1804.8810},
};

static skOpticalProperties_RayleighDryAir_Inelastic::sk_raman_xs_data o2linedata[] = {
	{-185.5861, 32, 1, 1606.3533, 0.3630},
	{-185.5690, 33, 1, 1608.0710, 0.3630},
	{-185.5512, 34, 1, 1605.8064, 0.3637},
	{-174.3154, 30, 1, 1420.7672, 0.3622},
	{-174.2980, 31, 1, 1422.5020, 0.3622},
	{-174.2802, 32, 1, 1420.2552, 0.3630},
	{-163.0162, 28, 1, 1246.4518, 0.3612},
	{-162.9980, 29, 1, 1248.2040, 0.3613},
	{-162.9809, 30, 1, 1245.9750, 0.3622},
	{-151.6906, 26, 1, 1083.4356, 0.3602},
	{-151.6730, 27, 1, 1085.2060, 0.3602},
	{-151.6551, 28, 1, 1082.9941, 0.3612},
	{-140.3405, 24, 1,  931.7450, 0.3589},
	{-140.3230, 25, 1,  933.5330, 0.3589},
	{-140.3046, 26, 1,  931.3390, 0.3602},
	{-128.9677, 22, 1,  791.4045, 0.3574},
	{-128.9490, 23, 1,  793.2100, 0.3574},
	{-128.9314, 24, 1,  791.0344, 0.3589},
	{-117.5740, 20, 1,  662.4368, 0.3555},
	{-117.5560, 21, 1,  664.2610, 0.3556},
	{-117.5372, 22, 1,  662.1030, 0.3574},
	{-108.2632, 19, 1,  546.7050, 0.0020},
	{-106.1613, 18, 1,  544.8628, 0.3533},
	{-106.1430, 19, 1,  546.7050, 0.3534},
	{-106.1240, 20, 1,  544.5658, 0.3555},
	{-104.3008, 18, 1,  544.8628, 0.0023},
	{ -96.8136, 17, 1,  440.5620, 0.0025},
	{ -94.7318, 16, 1,  438.7015, 0.3504},
	{ -94.7120, 17, 1,  440.5620, 0.3506},
	{ -94.6934, 18, 1,  438.4418, 0.3533},
	{ -92.8515, 16, 1,  438.7015, 0.0029},
	{ -85.3489, 15, 1,  345.8500, 0.0032},
	{ -83.2871, 14, 1,  343.9697, 0.3467},
	{ -83.2671, 15, 1,  345.8500, 0.3471},
	{ -83.2473, 16, 1,  343.7484, 0.3505},
	{ -81.3868, 14, 1,  343.9697, 0.0037},
	{ -73.8694, 13, 1,  262.5829, 0.0042},
	{ -71.8294, 12, 1,  260.6826, 0.3417},
	{ -71.8080, 13, 1,  262.5829, 0.3422},
	{ -71.7876, 14, 1,  260.5011, 0.3468},
	{ -69.9077, 12, 1,  260.6826, 0.0051},
	{ -62.3771, 11, 1,  190.7749, 0.0058},
	{ -60.3611, 10, 1,  188.8532, 0.3345},
	{ -60.3373, 11, 1,  190.7749, 0.3354},
	{ -60.3157, 12, 1,  188.7135, 0.3418},
	{ -58.4156, 10, 1,  188.8532, 0.0073},
	{ -50.8730,  9, 1,  130.4376, 0.0086},
	{ -48.8851,  8, 1,  128.4921, 0.3233},
	{ -48.8571,  9, 1,  130.4376, 0.3251},
	{ -48.8332, 10, 1,  128.3978, 0.3347},
	{ -46.9116,  8, 1,  128.4921, 0.0113},
	{ -40.1158,  4, 1,   42.2001, 0.0010},
	{ -39.3565,  7, 1,   81.5805, 0.0139},
	{ -37.4069,  6, 1,   79.6070, 0.3037},
	{ -37.3688,  7, 1,   81.5805, 0.3077},
	{ -37.3406,  8, 1,   79.5646, 0.3236},
	{ -35.3953,  6, 1,   79.6070, 0.0198},
	{ -27.8241,  5, 1,   44.2117, 0.0261},
	{ -25.9472,  4, 1,   42.2001, 0.2599},
	{ -25.8745,  5, 1,   44.2117, 0.2727},
	{ -25.8364,  6, 1,   42.2240, 0.3045},
	{ -25.8125,  4, 1,   42.2001, 0.0015},
	{ -23.8629,  4, 1,   42.2001, 0.0434},
	{ -16.2529,  3, 1,   18.3372, 0.0660},
	{ -16.2529,  2, 1,   16.2529, 0.0923},
	{ -14.3761,  3, 1,   18.3372, 0.1714},
	{ -14.3033,  4, 1,   16.3876, 0.2627},
	{ -14.1686,  2, 1,   16.2529, 0.0177},
	{ -12.2918,  2, 1,   16.2529, 0.1615},
	{  -2.1392, 19, 1,  546.7050, 0.0019},
	{  -2.1202, 17, 1,  440.5620, 0.0024},
	{  -2.1016, 15, 1,  345.8500, 0.0030},
	{  -2.0843,  3, 1,   18.3372, 0.0769},
	{  -2.0843,  2, 1,    2.0843, 0.1077},
	{  -2.0818, 13, 1,  262.5829, 0.0039},
	{  -2.0614, 11, 1,  190.7749, 0.0054},
	{  -2.0398,  9, 1,  130.4376, 0.0078},
	{  -2.0159,  7, 1,   81.5805, 0.0122},
	{  -2.0116,  5, 1,   44.2117, 0.0284},
	{  -1.9877,  5, 1,   44.2117, 0.0221},
	{  -1.9735,  7, 1,   81.5805, 0.0147},
	{  -1.9496,  3, 1,   18.3372, 0.0513},
	{  -1.9455,  9, 1,  130.4376, 0.0090},
	{  -1.9217, 11, 1,  190.7749, 0.0060},
	{  -1.9003, 13, 1,  262.5829, 0.0043},
	{  -1.8803, 15, 1,  345.8500, 0.0033},
	{  -1.8768,  1, 1,    3.9611, 0.2308},
	{  -1.8605, 17, 1,  440.5620, 0.0026},
	{  -1.8422, 19, 1,  546.7050, 0.0020},
	{  -0.1347,  4, 1,   16.3876, 0.0002},
	{   0.1347,  2, 1,   16.2529, 0.0004},
	{   1.8422, 18, 1,  544.8628, 0.0022},
	{   1.8605, 16, 1,  438.7015, 0.0027},
	{   1.8768,  2, 1,    2.0843, 0.1385},
	{   1.8803, 14, 1,  343.9697, 0.0035},
	{   1.9003, 12, 1,  260.6826, 0.0047},
	{   1.9217, 10, 1,  188.8532, 0.0066},
	{   1.9455,  8, 1,  128.4921, 0.0100},
	{   1.9496,  4, 1,   16.3876, 0.0399},
	{   1.9735,  6, 1,   79.6070, 0.0170},
	{   1.9877,  6, 1,   42.2240, 0.0187},
	{   2.0116,  4, 1,   42.2001, 0.0347},
	{   2.0159,  8, 1,   79.5646, 0.0108},
	{   2.0398, 10, 1,  128.3978, 0.0070},
	{   2.0614, 12, 1,  188.7135, 0.0049},
	{   2.0818, 14, 1,  260.5011, 0.0036},
	{   2.0843,  2, 1,   16.2529, 0.1077},
	{   2.0843,  0, 1,    0.0000, 0.5383},
	{   2.1016, 16, 1,  343.7484, 0.0028},
	{   2.1202, 18, 1,  438.4418, 0.0022},
	{   2.1392, 20, 1,  544.5658, 0.0018},
	{  12.2918,  1, 1,    3.9611, 0.2692},
	{  14.1686,  2, 1,    2.0843, 0.0177},
	{  14.3033,  2, 1,    2.0843, 0.4729},
	{  14.3761,  1, 1,    3.9611, 0.4000},
	{  16.2529,  2, 1,    2.0843, 0.0923},
	{  16.2529,  0, 1,    0.0000, 0.4617},
	{  23.8629,  3, 1,   18.3372, 0.0558},
	{  25.8125,  4, 1,   16.3876, 0.0015},
	{  25.8364,  4, 1,   16.3876, 0.4398},
	{  25.8745,  3, 1,   18.3372, 0.4286},
	{  25.9472,  2, 1,   16.2529, 0.4678},
	{  27.8241,  4, 1,   16.3876, 0.0319},
	{  35.3953,  5, 1,   44.2117, 0.0234},
	{  37.3406,  6, 1,   42.2240, 0.4232},
	{  37.3688,  5, 1,   44.2117, 0.4196},
	{  37.4069,  4, 1,   42.2001, 0.4387},
	{  39.3565,  6, 1,   42.2240, 0.0160},
	{  40.1158,  2, 1,    2.0843, 0.0019},
	{  46.9116,  7, 1,   81.5805, 0.0128},
	{  48.8332,  8, 1,   79.5646, 0.4134},
	{  48.8571,  7, 1,   81.5805, 0.4118},
	{  48.8851,  6, 1,   79.6070, 0.4228},
	{  50.8730,  8, 1,   79.5646, 0.0096},
	{  58.4156,  9, 1,  130.4376, 0.0080},
	{  60.3157, 10, 1,  128.3978, 0.4069},
	{  60.3373,  9, 1,  130.4376, 0.4060},
	{  60.3611,  8, 1,  128.4921, 0.4132},
	{  62.3771, 10, 1,  128.3978, 0.0064},
	{  69.9077, 11, 1,  190.7749, 0.0055},
	{  71.7876, 12, 1,  188.7135, 0.4022},
	{  71.8080, 11, 1,  190.7749, 0.4017},
	{  71.8294, 10, 1,  188.8532, 0.4068},
	{  73.8694, 12, 1,  188.7135, 0.0045},
	{  81.3868, 13, 1,  262.5829, 0.0040},
	{  83.2473, 14, 1,  260.5011, 0.3988},
	{  83.2671, 13, 1,  262.5829, 0.3985},
	{  83.2871, 12, 1,  260.6826, 0.4022},
	{  85.3489, 14, 1,  260.5011, 0.0034},
	{  92.8515, 15, 1,  345.8500, 0.0031},
	{  94.6934, 16, 1,  343.7484, 0.3961},
	{  94.7120, 15, 1,  345.8500, 0.3959},
	{  94.7318, 14, 1,  343.9697, 0.3988},
	{  96.8136, 16, 1,  343.7484, 0.0026},
	{ 104.3008, 17, 1,  440.5620, 0.0024},
	{ 106.1240, 18, 1,  438.4418, 0.3940},
	{ 106.1430, 17, 1,  440.5620, 0.3938},
	{ 106.1613, 16, 1,  438.7015, 0.3961},
	{ 108.2632, 18, 1,  438.4418, 0.0021},
	{ 115.7318, 19, 1,  546.7050, 0.0019},
	{ 117.5372, 20, 1,  544.5658, 0.3922},
	{ 117.5560, 19, 1,  546.7050, 0.3921},
	{ 117.5740, 18, 1,  544.8628, 0.3940},
	{ 119.6952, 20, 1,  544.5658, 0.0017},
	{ 128.9314, 22, 1,  662.1030, 0.3908},
	{ 128.9490, 21, 1,  664.2610, 0.3907},
	{ 128.9677, 20, 1,  662.4368, 0.3922},
	{ 140.3046, 24, 1,  791.0344, 0.3896},
	{ 140.3230, 23, 1,  793.2100, 0.3895},
	{ 140.3405, 22, 1,  791.4045, 0.3908},
	{ 151.6551, 26, 1,  931.3390, 0.3884},
	{ 151.6730, 25, 1,  933.5330, 0.3885},
	{ 151.6906, 24, 1,  931.7450, 0.3896},
	{ 162.9809, 28, 1, 1082.9941, 0.3876},
	{ 162.9980, 27, 1, 1085.2060, 0.3876},
	{ 163.0162, 26, 1, 1083.4356, 0.3884},
	{ 174.2802, 30, 1, 1245.9750, 0.3868},
	{ 174.2980, 29, 1, 1248.2040, 0.3868},
	{ 174.3154, 28, 1, 1246.4518, 0.3876},
	{ 185.5512, 32, 1, 1420.2552, 0.3861},
	{ 185.5690, 31, 1, 1422.5020, 0.3861},
	{ 185.5861, 30, 1, 1420.7672, 0.3868},
	{ 196.7919, 34, 1, 1605.8064, 0.3855},
	{ 196.8100, 33, 1, 1608.0710, 0.3855},
	{ 196.8269, 32, 1, 1606.3533, 0.3861}
};


/*---------------------------------------------------------------------------
 *'					skOpticalProperties_RayleighDryAir_Inelastic::skOpticalProperties_RayleighDryAir_Inelastic		2020-03-09 */
 /**	Construct the Rayleigh, Dry Air scattering object.
  *-------------------------------------------------------------------------*/

skOpticalProperties_RayleighDryAir_Inelastic::skOpticalProperties_RayleighDryAir_Inelastic()
{
	m_O2mix = 0.209476;
	m_N2mix = 0.78084;
	m_CO2mix = 0.000314;
	m_Armix = 0.00934;
	m_XXmix = 1.0 - (m_O2mix + m_N2mix + m_CO2mix + m_Armix);
	m_backgroundatmosphere = nullptr;
	m_temperaturekelvin = 250.0; // assume 250 K as a default
	CalculatePartitionFunctions(); 
}

/*---------------------------------------------------------------------------
 *'					skOpticalProperties_RayleighDryAir_Inelastic::~skOpticalProperties_RayleighDryAir_Inelastic		2020-03-09
 *-------------------------------------------------------------------------*/

skOpticalProperties_RayleighDryAir_Inelastic::~skOpticalProperties_RayleighDryAir_Inelastic()
{
	if (m_backgroundatmosphere != nullptr) m_backgroundatmosphere->Release();
}


 /*-----------------------------------------------------------------------------
  *					skOpticalProperties_RayleighDryAir_Inelastic::SetAtmosphericState		2020-03-09*/
  /** 
  * We are currently neglecting temperature and pressure effects.
   **/
   /*---------------------------------------------------------------------------*/

bool skOpticalProperties_RayleighDryAir_Inelastic::SetAtmosphericState(skClimatology* neutralatmosphere)
{
	if (neutralatmosphere != NULL) neutralatmosphere->AddRef();
	if (m_backgroundatmosphere != NULL) m_backgroundatmosphere->Release();
	m_backgroundatmosphere = neutralatmosphere;
	return true;
}


bool skOpticalProperties_RayleighDryAir_Inelastic::SetLocation(const GEODETIC_INSTANT& pt, bool* crosssectionschanged)
{
	bool ok = true;
	ok = ok && m_backgroundatmosphere != NULL;
	ok = ok && m_backgroundatmosphere->GetParameter(SKCLIMATOLOGY_TEMPERATURE_K, pt, &m_temperaturekelvin, false);
	ok = ok && CalculatePartitionFunctions();
	if (crosssectionschanged != NULL) *crosssectionschanged = true;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_RayleighDryAir_Inelastic::InternalClimatology_UpdateCache		2020-03-09*/
 /** The Rayleigh cross-section has no dependency on climatology parameters
  *	such as pressure and temperature. Thus this method does nothing.
  **/
  /*---------------------------------------------------------------------------*/

bool skOpticalProperties_RayleighDryAir_Inelastic::InternalClimatology_UpdateCache(const GEODETIC_INSTANT& /*pt*/)
{
	return true;
}

bool skOpticalProperties_RayleighDryAir_Inelastic::LookupUpThreadData(RayleighWavelength_TLS** data)
{
	size_t					threadnum;
	iterator				iter;
	bool					ok;
	static std::mutex		lock;							// Add a lock to make sure map.insert is thread safe.

	threadnum = nxWorkerThreadManager::GetCurrentThreadIdCode();
	lock.lock();
	iter = m_threadstate.find(threadnum);
	ok = (iter != m_threadstate.end());
	if (!ok)
	{
		std::pair<iterator, bool>							result;
		std::pair<size_t, RayleighWavelength_TLS >	newentry(threadnum, RayleighWavelength_TLS());

		result = m_threadstate.insert(newentry);
		ok = result.second;
		iter = result.first;
		(*iter).second.m_wavenumber = 0.0;
		(*iter).second.m_xsection = 0.0;
		(*iter).second.m_delta = 0.0;
		(*iter).second.m_deltaprime = 0.0;

	}
	lock.unlock();
	if (!ok)
	{
		*data = nullptr;
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_RayleighDryAir_Inelastic::LookupUpThreadData, error looking/creating thread data for thread %d", (int)threadnum);
	}
	else
	{
		*data = &(*iter).second;
	}
	return ok;

}

bool skOpticalProperties_RayleighDryAir_Inelastic::CalculateCrossSections(double wavenum, double* absxs, double* extxs, double* scattxs, skOpticalProperties_RayleighDryAir_Inelastic::RayleighWavelength_TLS* state, bool cabannes) const
{
	static const double PiCubed = 31.006276680299820175476315067101;		// This is PI**3
	static const double lochsmidt = 1.0 / 2.686763E19;						// Particles per cm3 in Ideal gas at 0C, 1 atmosphere.
	bool				ok;

	ok = (state->m_wavenumber == wavenum && state->m_cabannes == cabannes);
	if (!ok)
	{

		double  lamda = 1.0E4 / wavenum;		// get the wavelength in microns
		double  sigma = wavenum * 1.0E-4;		// Get the wavenumber in micron-1
		double  sigma2 = sigma * sigma;			// Get the terms for the 
		double  sigma4 = sigma2 * sigma2;		// polynomial expansions
		double  sigma6 = sigma2 * sigma4;
		double  sigma8 = sigma4 * sigma4;
		double	nr;
		double  xsection;
		double  O2xsect;
		double  N2xsect;
		double  Arxsect;
		double  CO2xsect;
		double  XXxsect;
		double	O2Fk;				// King Correction factor for O2
		double	N2Fk;				// King Correction factor for N2
		double	CO2Fk;				// King Correction factor for CO2
		double	ArFk;				// King Correction factor for Argon.
		double  XXFk;
		double	nO2;				// Refractivity of O2 at STP
		double	nN2;				// Refractivity of N2 at STP
		double	nAr;				// Refractivity of Argon at STP
		double	nCO2;				// Refractivity of CO2 at STP
		double	nXX;				// Refractivity of miscellaneous components at STP.

		// ---- Calculate Refractivities of individual components at STP
		// ---- formulae for refractivities are taken from Ref 1) Bates 1984
		// ---- and Ref 2)and Ref 3)
		//
		// ---- Note that more upto date values for the refractivity of air are
		// ---- available but not for the refractivity of the individual components.
		// ---- The biggest error in all of this is the lack of account for moist air.

		if (lamda < 0.254) nN2 = 1.0E-8*(6998.749 + 3233582.0 / (144 - sigma2));
		else if (lamda < 0.468) nN2 = 1.0E-8*(5989.242 + 3363266.3 / (144 - sigma2));
		else                       nN2 = 1.0E-8*(6855.200 + 3243157.0 / (144 - sigma2));

		if (lamda < 0.221)  nO2 = 1.0E-8*(23796.7 + 168988.4 / (40.9 - sigma2));
		else if (lamda < 0.288)  nO2 = 1.0E-8*(22120.4 + 203187.6 / (40.9 - sigma2));
		else if (lamda < 0.546)  nO2 = 1.0E-8*(20564.8 + 248089.9 / (40.9 - sigma2));
		else                       nO2 = 1.0E-8*(21351.1 + 218567.0 / (40.9 - sigma2));

		nAr = (5.547E-4*0.5)*(1.0
			+ 5.15E-3*sigma2
			+ 4.19E-5*sigma4
			+ 4.09E-7*sigma6
			+ 4.32E-9*sigma8);

		nCO2 = 1.0E-5*(1205.5*(5.79925 / (166.175 - sigma2))		// Refractivity of CO2 at STP
			+ 0.12005 / (79.609 - sigma2)
			+ 0.0053334 / (56.3064 - sigma2)
			+ 0.0043244 / (46.0196 - sigma2)
			+ 0.0001218145 / (0.0584738 - sigma2)
			);
		nXX = nAr;								// Refractivity of Miscellaneous constituents at STP (close)

		// ---- Calculate King Correction factors for different dry air constituents
		// ---- formulae for King Correction factors are taken from Ref 1) Bates 1984.

		O2Fk = 1.096 + 1.385E-3*sigma2 + 1.448E-4*sigma4;
		N2Fk = 1.034 + 3.17E-4*sigma2;
		CO2Fk = 1.15;
		ArFk = 1.0;
		XXFk = 1.0;

		// This change takes us from total Rayleigh to Cabannes according to Chance and Spurr (1997)
		if (cabannes)
		{
			O2Fk = 0.25 * (O2Fk + 3.0);
			N2Fk = 0.25 * (N2Fk + 3.0);
		}

		O2xsect = m_O2mix * nO2 *nO2 *O2Fk;			// (Unnormalized) fraction of cross-section from O2
		N2xsect = m_N2mix * nN2 *nN2 *N2Fk;			// (Unnormalized) fraction of cross-section from N2
		Arxsect = m_Armix * nAr *nAr *ArFk;			// (Unnormalized) fraction of cross-section from Ar
		CO2xsect = m_CO2mix * nCO2*nCO2*CO2Fk;			// (Unnormalized) fraction of cross-section from CO2
		XXxsect = m_XXmix * nXX *nXX *XXFk;			// (Unnormalized) fraction of cross-section from Miscellaneous species
		nr = O2xsect + N2xsect + Arxsect + CO2xsect + XXxsect;	// Total Unnormalized cross-section from everything

		xsection = 32.0*PiCubed / 3.0*nr*nxmath::sqr(wavenum*wavenum*lochsmidt);	// Get the Rayleigh cross-section in cm2.

		state->m_delta = 0;
		state->m_deltaprime = 0;
		AddDepolarization(O2Fk, O2xsect, state);
		AddDepolarization(N2Fk, N2xsect, state);
		AddDepolarization(ArFk, Arxsect, state);
		AddDepolarization(CO2Fk, CO2xsect, state);
		AddDepolarization(XXFk, XXxsect, state);
		state->m_delta /= nr;
		state->m_deltaprime /= nr;
		state->m_xsection = xsection;
		state->m_wavenumber = wavenum;
		state->m_cabannes = cabannes;
	}
	*absxs = 0.0;
	*extxs = state->m_xsection;
	*scattxs = state->m_xsection;
	return true;
}


bool skOpticalProperties_RayleighDryAir_Inelastic::CalculateCrossSections(double wavenum, double* absxs, double* extxs, double* scattxs)
{
	bool	ok;
	RayleighWavelength_TLS* threaddata;
	double inelxs;
	double atemp, etemp, stemp;

	ok = LookupUpThreadData(&threaddata);
	ok = ok && CalculateCrossSections(wavenum, absxs, extxs, &stemp, threaddata, false); // extinction is determined by cabannes + raman
	ok = ok && CalculateCrossSections(wavenum, absxs, &etemp, scattxs, threaddata, true); // scattering is now determined only by cabbanes
	return ok;
}


void skOpticalProperties_RayleighDryAir_Inelastic::AddDepolarization(double Fk, double fraction, skOpticalProperties_RayleighDryAir_Inelastic::RayleighWavelength_TLS* state) const
{
	double depol;

	depol = 6.0*(Fk - 1.0) / (3.0 + 7 * Fk);					// Get the depolarization term from the King correction factor
	state->m_delta += fraction * (1 - depol) / (1 + depol * 0.5);		// Get the DELTA term eqn 2.16 Hansen and Travis
	state->m_deltaprime += fraction * (1 - 2 * depol) / (1 - depol);			// Get the DELTA' term eqn 2.16 Hansen and Travis 
}


bool skOpticalProperties_RayleighDryAir_Inelastic::Interpolate_PhaseMatrixTables(double mu, skRTPhaseMatrix* P, skOpticalProperties_RayleighDryAir_Inelastic::RayleighWavelength_TLS* state) const
{
	double cosalpha = mu;
	double cosalpha2 = mu * mu;
	double sinalpha2 = 1 - cosalpha2;
	double cosalpha2p1 = 1 + cosalpha2;


	P->At(1, 1) = (SKRTFLOAT)(state->m_delta*0.75*cosalpha2p1);
	P->At(1, 2) = (SKRTFLOAT)(-state->m_delta*0.75*sinalpha2);
	P->At(1, 3) = (SKRTFLOAT)(0);
	P->At(1, 4) = (SKRTFLOAT)(0);

	P->At(2, 1) = (SKRTFLOAT)(P->At(1, 2));
	P->At(2, 2) = (SKRTFLOAT)(P->At(1, 1));
	P->At(2, 3) = (SKRTFLOAT)(0);
	P->At(2, 4) = (SKRTFLOAT)(0);

	P->At(3, 1) = (SKRTFLOAT)(0);
	P->At(3, 2) = (SKRTFLOAT)(0);
	P->At(3, 3) = (SKRTFLOAT)(state->m_delta*1.5*cosalpha);
	P->At(3, 4) = (SKRTFLOAT)(0);

	P->At(4, 1) = (SKRTFLOAT)(0);
	P->At(4, 2) = (SKRTFLOAT)(0);
	P->At(4, 3) = (SKRTFLOAT)(0);
	P->At(4, 4) = (SKRTFLOAT)(state->m_deltaprime*state->m_delta*1.5*cosalpha);

	P->At(1, 1) += (SKRTFLOAT)((1 - state->m_delta));
	return true;
}


bool skOpticalProperties_RayleighDryAir_Inelastic::CalculatePhaseMatrix(double wavenumber, double cosscatterangle, skRTPhaseMatrix* phasematrix)
{
	bool	ok;
	double  absxs, extxs, scattxs;
	RayleighWavelength_TLS* threaddata;

	ok = LookupUpThreadData(&threaddata);
	ok = ok && CalculateCrossSections(wavenumber, &absxs, &extxs, &scattxs, threaddata, true); // scattering is determined by cabannes only
	ok = ok && Interpolate_PhaseMatrixTables(cosscatterangle, phasematrix, threaddata);
	return ok;
}

bool skOpticalProperties_RayleighDryAir_Inelastic::CalculateInelasticCrossSections_Classical(double wavenumin, double* inelxs)
{
	bool ok = true;
	RayleighWavelength_TLS* threaddata;
	double absxs, extxs;
	double scattxscabannes, scattxsrayleigh;

	ok = ok && LookupUpThreadData(&threaddata);
	ok = ok && CalculateCrossSections(wavenumin, &absxs, &extxs, &scattxscabannes, threaddata, true);
	ok = ok && CalculateCrossSections(wavenumin, &absxs, &extxs, &scattxsrayleigh, threaddata, false);
	*inelxs = scattxsrayleigh - scattxscabannes;

	return ok;
}

bool skOpticalProperties_RayleighDryAir_Inelastic::CalculateInelasticCrossSections_Classical(double wavenumin, size_t lineidx, double* wavenumout, double* inelxs)
{
	bool ok = true;
	double totalclassical, totalquantum;

	ok = ok && CalculateInelasticCrossSections_Classical(wavenumin, &totalclassical);
	ok = ok && CalculateInelasticCrossSections_Quantum(wavenumin, &totalquantum);
	ok = ok && CalculateInelasticCrossSections_Quantum(wavenumin, lineidx, wavenumout, inelxs);
	*inelxs *= totalclassical / totalquantum;

	return ok;
}

bool skOpticalProperties_RayleighDryAir_Inelastic::CalculateInelasticCrossSections_Quantum(double wavenumin, double* inelxs)
{
	bool ok = true;
	size_t numlines = NumInelasticLines();

	double wn, xs;
	*inelxs = 0.0;
	for (size_t lineidx = 0; lineidx < numlines; lineidx++)
	{
		ok = ok && CalculateInelasticCrossSections_Quantum(wavenumin, lineidx, &wn, &xs);
		*inelxs += xs;
	}
	return ok;
}

bool skOpticalProperties_RayleighDryAir_Inelastic::CalculateInelasticCrossSections_Quantum(double wavenumin, size_t lineidx, double * wavenumout, double * inelxs)
{
	bool ok = true;

	static const double constant = 2901.5199742604463716763620507822;								// This is 256 * PI**5 / 27
	static const double planck = 6.62607015E-34;													// Planck constant in J.s
	static const double light = 299792458.0;														// Speed of light in m/s
	static const double boltzmann = 1.38064852E-23;													// Boltzmann constant J/K
	const double cmtoJoverkT = 100.0 * planck * light / (boltzmann * m_temperaturekelvin);	// factor to convert from energy E (cm^-1) to E/kT (unitless, where k is the Boltzmann constant and T is the temperature in Kelvin)

	size_t numo2 = N_ELEMENTS(o2linedata);
	size_t numn2 = N_ELEMENTS(n2linedata);

	int J, g;
	double E, dE, cPT, f, gamma2;

	if (lineidx < numo2)
	{
		size_t o2idx = lineidx;

		J = o2linedata[o2idx].J;
		E = o2linedata[o2idx].Eterm;
		cPT = o2linedata[o2idx].cPT;
		dE = o2linedata[o2idx].deltaE;

		*wavenumout = wavenumin - dE;											// deltaE is the change in molecular energy, equal to the LOSS in photon energy, therefore outgoing = incoming - deltaE
		f = (2 * J + 1) * exp(-E * cmtoJoverkT) * m_O2partition;				// fractional population of initial state
		gamma2 = CalculatePolarizabilityAnisotropyO2(wavenumin);				// polarizability anisotropy
		gamma2 *= gamma2;														// squared
		*inelxs = m_O2mix * constant * pow(*wavenumout, 4) * gamma2 * f * cPT;	// 1 * cm^-4 * cm^6 * 1 * 1 = cm^2
	}
	else
	{
		size_t n2idx = lineidx - numo2;
		ok = ok && n2idx < numn2;

		if (ok)
		{
			J = n2linedata[n2idx].J;
			E = n2linedata[n2idx].Eterm;
			cPT = n2linedata[n2idx].cPT;
			g = n2linedata[n2idx].gN;
			dE = n2linedata[n2idx].deltaE;

			*wavenumout = wavenumin - dE;											// deltaE is the change in molecular energy, equal to the LOSS in photon energy, therefore outgoing = incoming - deltaE
			f = (2 * J + 1) * g * exp(-E * cmtoJoverkT) * m_N2partition;			// fractional population of initial state
			gamma2 = CalculatePolarizabilityAnisotropyN2(wavenumin);				// polarizability anisotropy
			gamma2 *= gamma2;														// squared
			*inelxs = m_N2mix * constant * pow(*wavenumout, 4) * gamma2 * f * cPT;	// 1 * cm^-4 * cm^6 * 1 * 1 = cm^2
		}
	}
	return ok;
}

bool skOpticalProperties_RayleighDryAir_Inelastic::CalculateInelasticCrossSections_ClassicalReverse(double wavenumout, double* inelxs)
{
	bool ok = true;
	size_t numlines = NumInelasticLines();

	double wn, xs;
	*inelxs = 0.0;
	for (size_t lineidx = 0; lineidx < numlines; lineidx++)
	{
		ok = ok && CalculateInelasticCrossSections_ClassicalReverse(wavenumout, lineidx, &wn, &xs);
		*inelxs += xs;
	}
	return ok;
}

bool skOpticalProperties_RayleighDryAir_Inelastic::CalculateInelasticCrossSections_ClassicalReverse(double wavenumout, size_t lineidx, double* wavenumin, double* inelxs)
{
	bool ok = true;

	size_t numo2 = N_ELEMENTS(o2linedata);
	size_t numn2 = N_ELEMENTS(n2linedata);

	size_t o2idx, n2idx;
	double dE;
	if (lineidx < numo2)
	{
		o2idx = lineidx;
		dE = o2linedata[o2idx].deltaE;
		*wavenumin = wavenumout + dE; // deltaE is the change in molecular energy, equal to the LOSS in photon energy, therefore incoming = outgoing + deltaE
	}
	else
	{
		n2idx = lineidx - numo2;
		ok = ok && n2idx < numn2;
		if (ok)
		{
			double dE = n2linedata[n2idx].deltaE;
			*wavenumin = wavenumout + dE; // deltaE is the change in molecular energy, equal to the LOSS in photon energy, therefore incoming = outgoing + deltaE
		}
	}
	ok = ok && CalculateInelasticCrossSections_Classical(*wavenumin, lineidx, &wavenumout, inelxs);

	return ok;
}

bool skOpticalProperties_RayleighDryAir_Inelastic::CalculateInelasticCrossSections_QuantumReverse(double wavenumout, double* inelxs)
{
	bool ok = true;
	size_t numlines = NumInelasticLines();

	double wn, xs;
	*inelxs = 0.0;
	for (size_t lineidx = 0; lineidx < numlines; lineidx++)
	{
		ok = ok && CalculateInelasticCrossSections_QuantumReverse(wavenumout, lineidx, &wn, &xs);
		*inelxs += xs;
	}
	return ok;
}

bool skOpticalProperties_RayleighDryAir_Inelastic::CalculateInelasticCrossSections_QuantumReverse(double wavenumout, size_t lineidx, double* wavenumin, double* inelxs)
{
	bool ok = true;

	size_t numo2 = N_ELEMENTS(o2linedata);
	size_t numn2 = N_ELEMENTS(n2linedata);

	size_t o2idx, n2idx;
	double dE;

	if (lineidx < numo2)
	{
		size_t o2idx = lineidx;
		dE = o2linedata[o2idx].deltaE; 
		*wavenumin = wavenumout + dE; // deltaE is the change in molecular energy, equal to the LOSS in photon energy, therefore incoming = outgoing + deltaE
	}
	else
	{
		size_t n2idx = lineidx - numo2;
		ok = ok && n2idx < numn2;
		if (ok)
		{
			double dE = n2linedata[n2idx].deltaE;
			*wavenumin = wavenumout + dE; // deltaE is the change in molecular energy, equal to the LOSS in photon energy, therefore incoming = outgoing + deltaE
		}
	}

	ok = ok && CalculateInelasticCrossSections_Quantum(*wavenumin, lineidx, &wavenumout, inelxs);
	return ok;
}

bool skOpticalProperties_RayleighDryAir_Inelastic::CalculateInelasticCrossSections(double wavenumin, double* inelxs)
{
	return CalculateInelasticCrossSections_Classical(wavenumin, inelxs); // fast (1)
	//return CalculateInelasticCrossSections_Quantum(wavenumin, inelxs); // slow (N)
}


bool skOpticalProperties_RayleighDryAir_Inelastic::CalculateInelasticCrossSections(double wavenumin, size_t lineidx, double* wavenumout, double* inelxs)
{
	return CalculateInelasticCrossSections_Classical(wavenumin, lineidx, wavenumout, inelxs); // slow (N)
	//return CalculateInelasticCrossSections_Quantum(wavenumin, lineidx, wavenumout, inelxs); // fast (1)
}


bool skOpticalProperties_RayleighDryAir_Inelastic::CalculateInelasticCrossSections_Reverse(double wavenumout, double* inelxs)
{
	return CalculateInelasticCrossSections_ClassicalReverse(wavenumout, inelxs); // very slow (N^3)
	//return CalculateInelasticCrossSections_QuantumReverse(wavenumout, inelxs); // fast (1)
}

bool skOpticalProperties_RayleighDryAir_Inelastic::CalculateInelasticCrossSections_Reverse(double wavenumout, size_t lineidx, double* wavenumin, double* inelxs)
{
	/* for complete consistency this should call CalculateInelasticCrossSections_ClassicalReverse,
	but this would be significantly slower, and as of now these	values are only ever required relative 
	to each other -- currently trying it with classical */
	return CalculateInelasticCrossSections_ClassicalReverse(wavenumout, lineidx, wavenumin, inelxs); // very slow (N^2)
	//return CalculateInelasticCrossSections_QuantumReverse(wavenumout, lineidx, wavenumin, inelxs); // fast (1)
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_RayleighDryAir_Inelastic::CalculatePhaseMatrix		2020-03-09*/
 /** **/
 /*---------------------------------------------------------------------------*/

bool skOpticalProperties_RayleighDryAir_Inelastic::CalculateInelasticPhaseMatrix(double wavenumber, double cosscatterangle, skRTPhaseMatrix* phasematrix)
{
	double cosalpha = cosscatterangle;
	double cosalpha2 = cosscatterangle * cosscatterangle;
	double sinalpha2 = 1 - cosalpha2;
	double cosalpha2p1 = 1 + cosalpha2;

	// depolarization factor is 6/7 for Raman scattering
	double delta = 0.1;				// Get the DELTA term eqn 2.16 Hansen and Travis | delta = (1 - depol) / (1 + depol * 0.5) = (1/7) / (1 + 3/7) = 1/10
	double deltaprime = -5.0;		// Get the DELTA' term eqn 2.16 Hansen and Travis | delta = (1 - 2 * depol) / (1 - depol) = (-5/7) / (1/7) = -5


	phasematrix->At(1, 1) = (SKRTFLOAT)(delta*0.75*cosalpha2p1);
	phasematrix->At(1, 2) = (SKRTFLOAT)(-delta * 0.75*sinalpha2);
	phasematrix->At(1, 3) = (SKRTFLOAT)(0);
	phasematrix->At(1, 4) = (SKRTFLOAT)(0);

	phasematrix->At(2, 1) = (SKRTFLOAT)(phasematrix->At(1, 2));
	phasematrix->At(2, 2) = (SKRTFLOAT)(phasematrix->At(1, 1));
	phasematrix->At(2, 3) = (SKRTFLOAT)(0);
	phasematrix->At(2, 4) = (SKRTFLOAT)(0);

	phasematrix->At(3, 1) = (SKRTFLOAT)(0);
	phasematrix->At(3, 2) = (SKRTFLOAT)(0);
	phasematrix->At(3, 3) = (SKRTFLOAT)(delta*1.5*cosalpha);
	phasematrix->At(3, 4) = (SKRTFLOAT)(0);

	phasematrix->At(4, 1) = (SKRTFLOAT)(0);
	phasematrix->At(4, 2) = (SKRTFLOAT)(0);
	phasematrix->At(4, 3) = (SKRTFLOAT)(0);
	phasematrix->At(4, 4) = (SKRTFLOAT)(deltaprime*delta*1.5*cosalpha);

	phasematrix->At(1, 1) += (SKRTFLOAT)((1 - delta));
	return true;
}

size_t skOpticalProperties_RayleighDryAir_Inelastic::NumInelasticLines() const
{
	return N_ELEMENTS(o2linedata) + N_ELEMENTS(n2linedata);
}
 
/*-----------------------------------------------------------------------------
 *					skOpticalProperties_RayleighDryAir_Inelastic::CalculatePartitionFunction		2020-03-09*/
 /** **/
 /*---------------------------------------------------------------------------*/
bool skOpticalProperties_RayleighDryAir_Inelastic::CalculatePartitionFunctions()
{
	static const double planck = 6.62607015E-34;							// Panck constant in J.s
	static const double light = 299792458.0;								// Speed of light in m/s
	static const double boltzmann = 1.38064852E-23;							// Boltzmann constant J/K
																
	double J;																// total angular momentum quantum number
	double g;																// nuclear spin statistical weight factor
	double E;																// state energy in cm^-1

	// factor to convert from energy E (cm^-1) to E/kT (unitless, where k is the Boltzmann constant and T is the temperature in Kelvin)
	double cmtoJoverkT = 100.0 * planck * light / (boltzmann * m_temperaturekelvin);		

	m_O2partition = 0;
	for (size_t idx = 0; idx < N_ELEMENTS(o2partitiondata); idx++)
	{
		J = o2partitiondata[idx].J;
		E = o2partitiondata[idx].Eterm;
		m_O2partition += (2 * J + 1) * exp(-E * cmtoJoverkT);
	}

	m_N2partition = 0;
	for (size_t idx = 0; idx < N_ELEMENTS(n2partitiondata); idx++)
	{
		J = n2partitiondata[idx].J;
		g = n2partitiondata[idx].gN;
		E = n2partitiondata[idx].Eterm;
		m_N2partition += g * (2 * J + 1) * exp(-E * cmtoJoverkT);
	}

	m_O2partition = 1.0 / m_O2partition;
	m_N2partition = 1.0 / m_N2partition;

	return true;
}

double skOpticalProperties_RayleighDryAir_Inelastic::CalculatePolarizabilityAnisotropyO2(double excitationwavenumber)
{
	double excitationwavenumber_um_squared = excitationwavenumber * excitationwavenumber * 1E-8;
	return 1E-24 * (0.07149 + 45.9364 / (48.2716 - excitationwavenumber_um_squared));
}

double skOpticalProperties_RayleighDryAir_Inelastic::CalculatePolarizabilityAnisotropyN2(double excitationwavenumber)
{
	double excitationwavenumber_um_squared = excitationwavenumber * excitationwavenumber * 1E-8;
	return 1e-25 * (-6.01466 + 2385.57 / (186.099 - excitationwavenumber_um_squared));
}

