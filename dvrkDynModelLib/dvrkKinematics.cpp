// Project Header files
#include "dvrkKinematics.hpp"
#include <iostream>

/**
* @brief Denavit-Hartenberg function
* Fills the DH parameters of the dVRK PSM arm according to a modified DH convention, with kinematic simplification
* @param q: the current full-dimension (i.e., including passive joints) joint configuration of the dVRK PSM arm
* @return the matrix with DH parameters
*/
Eigen::Matrix<float, PSM_ACTIVE_JOINTS, 5> psmSimplifiedModDH(const psmFullJointsVectorf& q){

    //   parameters are all expressed in meters (m)
    double l_RCC = 0.4318;
    double l_2L2 = 0.516;
    double l_p2y = 0.0091;
    double l_tool = 0.4162;
    double l_clip = 0.0091;
    double pi = double(M_PI);


    // joint type | d_i | theta_i | a_i  | alpha_i
    Eigen::Matrix<float, PSM_ACTIVE_JOINTS, 5> DH;

    DH << 1,        -l_2L2,  q(0) + pi / 2,     0,  pi / 2,
          1,             0,  q(1) - pi / 2,     0, -pi / 2,
          2, -l_RCC + q(7),             pi,     0,  pi / 2,
          1,        l_tool,           q(9),     0,       0,
          1,             0, q(10) + pi / 2,     0,  pi / 2,
          1,             0, q(11) + pi / 2, l_p2y, -pi / 2,
          1,       -l_clip,             pi,     0, -pi / 2;

    return DH;
}


/**
* @brief direct kinematics function
* Build the chain of homogeneous transformation matrices from the input DH parameters
* @param DH: the set of DH parameters
* @return the chain of resulting homogeneous transformations
*
*/
std::vector<Eigen::Matrix4f> chainOfTransformations(const Eigen::Matrix<float,PSM_ACTIVE_JOINTS,5>& DH) {

    std::vector<Eigen::Matrix4f> chain;
    Eigen::Matrix4f Ti;

    for (int k = 0; k < DH.rows(); k++) {

        Ti << cos(DH(k, 2)), -sin(DH(k, 2)), 0, DH(k, 3),
            cos(DH(k, 4))* sin(DH(k, 2)), cos(DH(k, 4))* cos(DH(k, 2)), -sin(DH(k, 4)), -sin(DH(k, 4)) * DH(k, 1),
            sin(DH(k, 4))* sin(DH(k, 2)), sin(DH(k, 4))* cos(DH(k, 2)), cos(DH(k, 4)), cos(DH(k, 4))* DH(k, 1),
            0, 0, 0, 1;

        chain.push_back(Ti);

    }


    return chain;
}


/**
* @brief kinematics function
* Compute the direct kinematics function and the Jacobian matrix from the input DH parameters (multiplying the chain transformations)
* @param q: the full joint configuration vector
* @param T0: the initial transformation world-base (default is identity matrix)
* @return a pair containing the chain of the homogeneous transformation of the direct kinematics and the Jacobian matrix
*
*/
std::pair < std::vector<Eigen::Matrix4f>, Eigen::Matrix<float, 6, PSM_ACTIVE_JOINTS> > psmKinematics(const psmFullJointsVectorf& q, const Eigen::Matrix4f& T0){

    std::pair < std::vector<Eigen::Matrix4f>, Eigen::Matrix<float, 6, PSM_ACTIVE_JOINTS> > kin;
    std::vector<Eigen::Matrix4f> chain;
    std::vector<Eigen::Matrix4f> dirKinMats;
    Eigen::Matrix<float, 6, PSM_ACTIVE_JOINTS> J;
    J.setZero();

    // Build the direct kinematics of the PSM arm, according to the modified DH convention with simplifications
    Eigen::Matrix<float, PSM_ACTIVE_JOINTS, 5> DH = psmSimplifiedModDH(q);

    // Compute the chain of relative transformations from the DH parameters
    chain = chainOfTransformations(DH);

    // Build the direct kinematics matrices
    dirKinMats.resize(chain.size());
    dirKinMats[0] = T0 * chain[0];
    for (int k = 1; k < PSM_ACTIVE_JOINTS; k++) {
    
        dirKinMats[k] = dirKinMats[k - 1] * chain[k];
    
    }

    // Build the Jacobian matrix
    Eigen::Matrix4f Twee = dirKinMats[PSM_ACTIVE_JOINTS - 1];
    Eigen::Vector3f pee = Twee.block<3, 1>(0, 3);
    for (int i = 0; i < PSM_ACTIVE_JOINTS; i++) {
    
        Eigen::Matrix4f Twi = dirKinMats[i];

        Eigen::Vector3f z_i = Twi.block<3, 1>(0, 2);
        if (DH(i, 0) == 1) {
            Eigen::Vector3f p_i = Twi.block<3, 1>(0, 3);
            J.block<3, 1>(0, i) = z_i.cross(pee - p_i);
            J.block<3, 1>(3, i) = z_i;
        }
        else {
            J.block<3, 1>(0, i) = z_i;
            J.block<3, 1>(3, i).setZero();
        }

    }

    // Fill the output pair
    kin.first.resize(dirKinMats.size());
    kin.first = dirKinMats;
    kin.second = J;

    return kin; // substitute with dirKinMats
}

