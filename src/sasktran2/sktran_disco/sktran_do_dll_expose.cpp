#include <sktran_disco/sktran_do_lowlevelinterface.h>

extern "C" void calculate(const sasktran_disco_lowlevel::Atmosphere* atmosphere, const sasktran_disco_lowlevel::Config* config,
                          const sasktran_disco_lowlevel::WeightingFunctions* weightingfunctions,
                          const sasktran_disco_lowlevel::ViewingGeometry* geometry, const sasktran_disco_lowlevel::Output* output) {
    sasktran_disco_lowlevel::calculate(atmosphere, config, weightingfunctions, geometry, output);
}