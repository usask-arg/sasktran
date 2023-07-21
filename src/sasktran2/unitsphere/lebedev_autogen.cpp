#include <sasktran2/internal_common.h>
#include <sasktran2/math/unitsphere/lebedev.h>

#include <sasktran2/math/unitsphere/lebedev/sphere_6.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_14.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_26.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_38.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_50.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_74.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_86.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_110.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_146.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_170.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_194.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_230.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_266.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_302.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_350.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_434.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_590.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_770.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_974.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_1202.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_1454.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_1730.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_2030.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_2354.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_2702.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_3074.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_3470.h>
#include <sasktran2/math/unitsphere/lebedev/sphere_3890.h>

namespace sasktran2::math::unitsphere::lebedev {
    void get_lebedev_data(int npoints, Eigen::MatrixXd& result) {
        if(npoints == 6) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_6, 4, 6);
        } else if(npoints == 14) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_14, 4, 14);
        } else if(npoints == 26) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_26, 4, 26);
        } else if(npoints == 38) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_38, 4, 38);
        } else if(npoints == 50) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_50, 4, 50);
        } else if(npoints == 74) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_74, 4, 74);
        } else if(npoints == 86) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_86, 4, 86);
        } else if(npoints == 110) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_110, 4, 110);
        } else if(npoints == 146) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_146, 4, 146);
        } else if(npoints == 170) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_170, 4, 170);
        } else if(npoints == 194) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_194, 4, 194);
        } else if(npoints == 230) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_230, 4, 230);
        } else if(npoints == 266) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_266, 4, 266);
        } else if(npoints == 302) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_302, 4, 302);
        } else if(npoints == 350) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_350, 4, 350);
        } else if(npoints == 434) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_434, 4, 434);
        } else if(npoints == 590) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_590, 4, 590);
        } else if(npoints == 770) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_770, 4, 770);
        } else if(npoints == 974) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_974, 4, 974);
        } else if(npoints == 1202) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_1202, 4, 1202);
        } else if(npoints == 1454) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_1454, 4, 1454);
        } else if(npoints == 1730) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_1730, 4, 1730);
        } else if(npoints == 2030) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_2030, 4, 2030);
        } else if(npoints == 2354) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_2354, 4, 2354);
        } else if(npoints == 2702) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_2702, 4, 2702);
        } else if(npoints == 3074) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_3074, 4, 3074);
        } else if(npoints == 3470) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_3470, 4, 3470);
        } else if(npoints == 3890) {
            result = Eigen::Map<Eigen::MatrixXd>(g_lebedev_xyzw_3890, 4, 3890);
        } else {
            BOOST_LOG_TRIVIAL(error) << "Requested number of Lebedev quadrature points does not exist";
            throw std::runtime_error("Requested number of Lebedev quadrature points does not exist");
        }
    };
}
