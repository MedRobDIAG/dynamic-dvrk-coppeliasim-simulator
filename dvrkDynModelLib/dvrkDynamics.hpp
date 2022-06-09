/*
 *
 *  Created on: 02 feb 2021
 *     Authors: Marco Ferro, Alessandro Mirante
 *
 *  M. Ferro, A. Mirante, F. Ficuciello, M. Vendittelli, 'A CoppeliaSim 
 *  dynamic simulator for the daVinci Research Kit'. IEEE RA-L, 2022.
 *
 */
#pragma once
#include "dvrkDefines.hpp"

//static const double  l9x = 0.0091043, l9y = -0.0064588, l9z = -0.0348129;

static const double  L1xx = 1.318395, L1xy = 4.47e-5, L1xz = 0.0, L1yy = 1.3192812, L1yz = 0.0, L1zz = 0.00093;
static const double l1x = 0.0735011, l1y = -0.0036854, l1z = 0.0;

static const double  L2xx = 0.0769277, L2xy = -0.0156765, L2xz = 0.0069185, L2yy = 0.0316068, L2yz = -0.0039539, L2zz = 0.0603379;
static const double l2x = -0.1951434, l2y = -0.7218708, l2z = -0.1108571;

static const double  L4xx = 0.0021128, L4xy = -0.0103894, L4xz = -0.0090333, L4yy = 90.1865526, L4yz = -0.001028, L4zz = 90.1868409;
static const double  l4x = 0.5180607, l4y = 0.0591317, l4z = 0.0514193;

static const double  L5xx = 0.0009005, L5xy = 0.0042873, L5xz = -0.004286, L5yy = 37.088943, L5yz = 0.0004319, L5zz = 37.0889432;
static const double  l5x = 0.1428667, l5y = -0.0143942, l5z = 0.01439;

static const double  L6xx = 0.0003618, L6xy = 2.07e-5, L6xz = 3.59e-5, L6yy = 0.0001666, L6yz = 0.0001348, L6zz = 0.0003194;
static const double  l6x = -1e-7, l6y = 0.0050002, l6z = -0.0029356;

static const double  L7xx = 0.0100262, L7xy = -0.0031696, L7xz = 0.0004199, L7yy = 0.0013188, L7yz = 0.0014119, L7zz = 0.0109082;
static const double  l7x = -0.022787, l7y = -0.070001, l7z = 0.0100001;

static const double  L8xx = 0.0037648, L8xy = -0.0001809, L8xz = -0.0007013, L8yy = 0.0038083, L8yz = -0.0007102, L8zz = 0.0003704;
static const double  l8x = -0.0072563, l8y = -0.0072582, l8z = -0.0361225;

//static const double  L9xx = 0.0036083, L9xy = 0.000204, L9xz = 0.0009268, L9yy = 0.0037699, L9yz = -0.0006548, L9zz = 0.0004368;

static const double m1 = 6.0667823;
static const double m2 = 10.0001371;
static const double m4 = 2.956747;
static const double m5 = 0.4796773;
static const double m6 = 0.1000014;
static const double m7 = 0.5000069;
static const double m8 = 0.3630268;
//static const double m9 = 0.3499395;


/**
* @brief Mass matrix function
* Compute the mass matrix of the DVRK PSM arm dynamic model
* @param [out] M: the mass matrix of the PSM 
* @param [in] q: the current joint positions of the PSM
* @param [in] withCounterWeight: if the counterweight is modeled and accounted in the computation or not
* @param [in] withFriction: if the joint friction is modeled and accounted in the computation or not
*
*/
psmMassMatrixf computeDVRKPSMMassMatrix(const psmActiveJointsVectorf& q, const bool& withCounterWeight = true, bool withFriction = false);

/**
* @brief Gravity vector function
* Compute the gravity vector of the DVRK PSM arm dynamic model
* @param [out] g: the gravity of the PSM
* @param [in] q: the current joint positions of the PSM
* @param [in] withCounterWeight: if the counterweight is modeled and accounted in the computation or not
* @param [in] withFriction: if the joint friction is modeled and accounted in the computation or not
*
*/
psmActiveJointsVectorf computeDVRKPSMGravityVector(const psmActiveJointsVectorf& q, const bool& withCounterWeight = true, bool withFriction = false);

/**
* @brief Coriolis vector function
* Compute the Coriolis and centrifugal term vector of the DVRK PSM arm dynamic model
* @param [out] C: the Coriolis and centrifugal term vector of the PSM
* @param [in] q: the current joint positions of the PSM
* @param [in] dq: the current joint velocities of the PSM
* @param [in] withCounterWeight: if the counterweight is modeled and accounted in the computation or not
* @param [in] withFriction: if the joint friction is modeled and accounted in the computation or not
*
*/
psmActiveJointsVectorf computeDVRKPSMCoriolisVector(const psmActiveJointsVectorf& q, const psmActiveJointsVectorf& dq, const bool& withCounterWeight = true, bool withFriction = false);

/**
* @brief Predicted torque vector function
* Compute the model-based torque from the full regression matrix
* @param [out] tau: the predicted torque vector of the PSM
* @param [in] q: the current joint positions of the PSM
* @param [in] dq: the current joint velocities of the PSM
* @param [in] ddq: the current joint accelerations of the PSM
* @param [in] withCounterWeight: if the counterweight is modeled and accounted in the computation or not
* @param [in] withFriction: if the joint friction is modeled and accounted in the computation or not
*
*/
psmActiveJointsVectorf computePredictedTorque(const psmActiveJointsVectorf& q, const psmActiveJointsVectorf& dq, const psmActiveJointsVectorf& ddq, const bool& withCounterWeight = true, bool withFriction = false);