#pragma once

// Eigen Header files
#include <unsupported/Eigen/CXX11/Tensor>

// dvrk Header files
#include "dvrkDefines.hpp"

/**
* @brief Denavit-Hartenberg function
* Fills the DH parameters of the dVRK PSM arm according to a modified DH convention, with kinematic simplification
* @param q: the current full-dimension (i.e., including passive joints) joint configuration of the dVRK PSM arm
* @return the matrix with DH parameters
*/
Eigen::Matrix<float, PSM_ACTIVE_JOINTS, 5> psmSimplifiedModDH(const psmFullJointsVectorf& q);

/**
* @brief chain of transformation function
* Build the chain of homogeneous transformation matrices from the input DH parameters
* @param DH: the set of DH parameters
* @return the chain of resulting homogeneous transformations
*
*/
std::vector<Eigen::Matrix4f> chainOfTransformations(const Eigen::Matrix<float, PSM_ACTIVE_JOINTS, 5>& DH);

/**
* @brief kinematics function
* Compute the direct kinematics function and the Jacobian matrix from the input DH parameters (multiplying the chain transformations)
* @param q: the full joint configuration vector
* @param T0: the initial transformation world-base (default is identity matrix)
* @return a pair containing the chain of the homogeneous transformation of the direct kinematics and the Jacobian matrix
*
*/
std::pair < std::vector<Eigen::Matrix4f>, Eigen::Matrix<float, 6, PSM_ACTIVE_JOINTS> > psmKinematics(const psmFullJointsVectorf& q, const Eigen::Matrix4f& T0);
