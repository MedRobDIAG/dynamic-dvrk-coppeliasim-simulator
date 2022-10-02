#pragma once

// Standard header files
#include <math.h>

// Eigen header files
#include <Eigen/Dense>


#ifndef PSM_ACTIVE_JOINTS
#define PSM_ACTIVE_JOINTS 7
#endif

#ifndef PSM_FULL_JOINTS
#define PSM_FULL_JOINTS 14
#endif

#ifndef PSM_DYN_PARAMS
#define PSM_DYN_PARAMS 115
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef Eigen::Matrix<double, PSM_ACTIVE_JOINTS, 1> psmActiveJointsVectord;
typedef Eigen::Matrix<int, PSM_ACTIVE_JOINTS, 1> psmActiveJointsVectori;
typedef Eigen::Matrix<float, PSM_ACTIVE_JOINTS, 1> psmActiveJointsVectorf;
typedef Eigen::Matrix<double, PSM_FULL_JOINTS, 1> psmFullJointsVectord;
typedef Eigen::Matrix<float, PSM_FULL_JOINTS, 1> psmFullJointsVectorf;
typedef Eigen::Matrix<int, PSM_FULL_JOINTS, 1> psmFullJointsVectori;
typedef Eigen::Matrix<float, PSM_DYN_PARAMS, 1> psmDynParamsVectorf;
typedef Eigen::Matrix<double, PSM_DYN_PARAMS, 1> psmDynParamsVectord;
typedef Eigen::Matrix<float, PSM_ACTIVE_JOINTS, PSM_DYN_PARAMS> psmRegressorMatrixf;
typedef Eigen::Matrix<double, PSM_ACTIVE_JOINTS, PSM_DYN_PARAMS> psmRegressorMatrixd;
typedef Eigen::Matrix<float, PSM_ACTIVE_JOINTS, PSM_ACTIVE_JOINTS> psmMassMatrixf;
typedef Eigen::Matrix<double, PSM_ACTIVE_JOINTS, PSM_ACTIVE_JOINTS> psmMassMatrixd;

namespace Eigen {

	typedef Eigen::Matrix<float, 6, 1> Vector6f;
	typedef Eigen::Matrix<double, 6, 1> Vector6d;
	typedef Eigen::Matrix<int, 6, 1> Vector6i;

}
