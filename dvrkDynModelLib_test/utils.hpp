#pragma once

// Eigen Header files
#include <Eigen/Dense>

// Standard Header files
#include <vector>
#include <iostream>
#include <fstream>
#include <conio.h>
#include <string>

#define _USE_MATH_DEFINES
#include <math.h>

/**
* @brief Load function
* Load and extract data from the input log
* @param file: the log file
* @return a structure containing the loaded data
*/
std::vector < Eigen::VectorXf > loadCSVLog(const std::string& file);

/**
* @brief Read function
* Read the conent of the file, assuming it's writte in the CSV format
* @ return a vector of string containing the set of values of the read CSV structure
*/
std::vector < std::string > readCSVFile(const char* filename);

/**
* @brief Utility function
* parse the input string line into an array of doubled-precision floating values
* @param the input string line
* @return the corresponding vector of double-precision floating values
*/
std::vector < double > parseCSVLine(const std::string& line);

/**
* @brief Filter function
* Filter the input velocity vector with a low-pass filter
* @param velocity: the input velocity vector to be filtered
* @param velocity_state: the latest filteret velocity
* @return the filtered velocity
*/
Eigen::VectorXf filter_velocity(Eigen::VectorXf velocity, Eigen::Ref<Eigen::VectorXf> velocity_state);

/**
* @brief Filter function
* Filter the input velocity vector with a low-pass filter
* @param velocity: the input velocity vector to be filtered
* @param velocity_state: the latest filteret velocity
* @return the filtered velocity
*/
Eigen::VectorXf filter_accelleration(Eigen::VectorXf accelleration, Eigen::Ref<Eigen::VectorXf> accelleration_state);

/**
* @brief Filter function
* Filter the input velocity vector with a low-pass filter
* @param velocity: the input velocity vector to be filtered
* @param velocity_state: the latest filteret velocity
* @return the filtered velocity
*/
Eigen::VectorXf filter_Torques(Eigen::VectorXf torques, Eigen::Ref<Eigen::VectorXf> torques_state);


/**
* @brief Rotation matrix 2 Roll Pitch Yaw angles
* Convert the input rotation matrix in the corresponding triple of roll pitch yaw angles
*/
Eigen::Vector3f rot2rpy(const Eigen::Matrix3f& R);
