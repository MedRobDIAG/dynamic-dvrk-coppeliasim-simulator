// Standard header files
#include <utility>
#include <limits.h>

// Library header files
#include "dvrkDynamics.hpp"

#include <iostream>

template <typename T> int sign(T val) {
	return (T(0) < val) - (val < T(0));
}

/**
* @brief Mass matrix function
* Compute the mass matrix of the DVRK PSM arm dynamic model
* @param [out] M: the mass matrix of the PSM
* @param [in] q: the current joint positions of the PSM
* @param [in] withCounterWeight: if the counterweight is modeled and accounted in the computation or not
* @param [in] withFriction: if the joint friction is modeled and accounted in the computation or not
*
*/
psmMassMatrixf computeDVRKPSMMassMatrix(const psmActiveJointsVectorf& q, const bool& withCounterWeight, bool withFriction){

	psmMassMatrixf M;
	M.setZero();

	double m9 = 0.0, l9x = 0.0, l9y = 0.0, l9z = 0.0, L9xx = 0.0, L9xy = 0.0, L9xz = 0.0, L9yy = 0.0, L9yz = 0.0, L9zz = 0.0;
	double Ia10 = 0.0, Ia11 = 0.0, Ia14 = 0.0, Ia15 = 0.0;

	double q1 = q(0);
	double q2 = q(1);
	double q3 = q(2);
	double q4 = q(3);

	double s1 = sin(q1);
	double c1 = cos(q1);
	double s2 = sin(q2);
	double c2 = cos(q2);
	double c12 = pow(c1, 2);
	double s12 = pow(s1, 2);
	double s22 = pow(s2, 2);
	double s23 = pow(s2, 3);
	double c22 = pow(c2, 2);
	double c23 = pow(c2, 3);
	double q32 = pow(q3, 2);

	if (withCounterWeight) {	
		L9xx = 0.0036083;
		L9xy = 0.000204;
		L9xz = 0.0009268;
		L9yy = 0.0037699;
		L9yz = -0.0006548;
		L9zz = 0.0004368;
		m9 = 0.3499395;
		l9x = 0.0091043;
		l9y = -0.0064588;
		l9z = -0.0348129;
	}

	if (withFriction) {
		Ia10 = -3.1e-6;
		Ia11 = -0.0006152;
		Ia14 = -0.0002932;
		Ia15 = -0.0003634;
	}

	M(0,0) = L1zz + L2xx * c22 - 2 * L2xy * s2 * c2 + L2yy * s22 + 1.0 * L4xx + 1.0 * L5xx + 1.0 * L6xx * c22 - 2.0 * L6xy * s2 * c2 + 1.0 * L6yy * s22 + 1.0 * L7xx * c22 - 2.0 * L7xy * s2 * c2 + 1.0 * L7yy * s22 + 1.0 * L8xx * c22 - 2.0 * L8xz * s2 * c2 + 1.0 * L8zz * s22 + 1.0 * L9xx * c22 - 2.0 * L9xz * s2 * c2 + 1.0 * L9zz * s22 - 0.08018 * l4y * s12 * s2 - 0.28908 * l4y * s12 * c2 - 0.08018 * l4y * s2 * c12 - 0.28908 * l4y * c12 * c2 - 0.08018 * l5y * s12 * s2 - 0.36524 * l5y * s12 * c2 - 0.08018 * l5y * s2 * c12 - 0.36524 * l5y * c12 * c2 - 0.08018 * l6x * s12 * s22 - 0.28908 * l6x * s12 * s2 * c2 - 0.08018 * l6x * s22 * c12 - 0.28908 * l6x * s2 * c12 * c2 - 0.08018 * l6y * s12 * s2 * c2 - 0.28908 * l6y * s12 * c22 - 0.08018 * l6y * s2 * c12 * c2 - 0.28908 * l6y * c12 * c22 + 0.08018 * l7x * s12 * s22 + 0.28908 * l7x * s12 * s2 * c2 + 0.08018 * l7x * s22 * c12 + 0.28908 * l7x * s2 * c12 * c2 + 0.08018 * l7y * s12 * s2 * c2 + 0.28908 * l7y * s12 * c22 + 0.08018 * l7y * s2 * c12 * c2 + 0.28908 * l7y * c12 * c22 + 2.0 * l8x * q3 * s12 * s2 * c2 + 2.0 * l8x * q3 * s2 * c12 * c2 - 0.8636 * l8x * s12 * s2 * c2 - 0.8636 * l8x * s2 * c12 * c2 + 2.0 * l8z * q3 * s12 * c22 + 2.0 * l8z * q3 * c12 * c22 - 0.8636 * l8z * s12 * c22 - 0.8636 * l8z * c12 * c22 + 2.0 * l9x * q3 * s12 * s2 * c2 + 2.0 * l9x * q3 * s2 * c12 * c2 + 0.08018 * l9x * s12 * s22 + 0.08018 * l9x * s22 * c12 + 2.0 * l9z * q3 * s12 * c22 + 2.0 * l9z * q3 * c12 * c22 + 0.08018 * l9z * s12 * s2 * c2 + 0.08018 * l9z * s2 * c12 * c2 + 0.0016072081 * m4 * s12 * s22 + 0.0115892172 * m4 * s12 * s2 * c2 + 0.0208918116 * m4 * s12 * c22 + 0.0016072081 * m4 * s22 * c12 + 0.0115892172 * m4 * s2 * c12 * c2 + 0.0208918116 * m4 * c12 * c22 + 0.0016072081 * m5 * s12 * s22 + 0.0146424716 * m5 * s12 * s2 * c2 + 0.0333500644 * m5 * s12 * c22 + 0.0016072081 * m5 * s22 * c12 + 0.0146424716 * m5 * s2 * c12 * c2 + 0.0333500644 * m5 * c12 * c22 + 0.0016072081 * m6 * s12 * s22 + 0.0115892172 * m6 * s12 * s2 * c2 + 0.0208918116 * m6 * s12 * c22 + 0.0016072081 * m6 * s22 * c12 + 0.0115892172 * m6 * s2 * c12 * c2 + 0.0208918116 * m6 * c12 * c22 + 0.0016072081 * m7 * s12 * s22 + 0.0115892172 * m7 * s12 * s2 * c2 + 0.0208918116 * m7 * s12 * c22 + 0.0016072081 * m7 * s22 * c12 + 0.0115892172 * m7 * s2 * c12 * c2 + 0.0208918116 * m7 * c12 * c22 + 1.0 * m8 * q32 * s12 * c22 + 1.0 * m8 * q32 * c12 * c22 - 0.8636 * m8 * q3 * s12 * c22 - 0.8636 * m8 * q3 * c12 * c22 + 0.18645124 * m8 * s12 * c22 + 0.18645124 * m8 * c12 * c22 + 1.0 * m9 * q32 * s12 * c22 + 1.0 * m9 * q32 * c12 * c22 + 0.08018 * m9 * q3 * s12 * s2 * c2 + 0.08018 * m9 * q3 * s2 * c12 * c2 + 0.0016072081 * m9 * s12 * s22 + 0.0016072081 * m9 * s22 * c12;
	M(0,1) = L2xz * c2 - L2yz * s2 - 1.0 * L6xz * c2 + 1.0 * L6yz * s2 + 1.0 * L7xz * c2 - 1.0 * L7yz * s2 + 1.0 * L8xy * c2 - 1.0 * L8yz * s2 - 1.0 * L9xy * c2 + 1.0 * L9yz * s2 + 0.14454 * l4z * s12 * s2 - 0.04009 * l4z * s12 * c2 + 0.14454 * l4z * s2 * c12 - 0.04009 * l4z * c12 * c2 + 0.18262 * l5z * s12 * s2 - 0.04009 * l5z * s12 * c2 + 0.18262 * l5z * s2 * c12 - 0.04009 * l5z * c12 * c2 + 0.14454 * l6z * s12 * s2 - 0.04009 * l6z * s12 * c2 + 0.14454 * l6z * s2 * c12 - 0.04009 * l6z * c12 * c2 + 0.14454 * l7z * s12 * s2 - 0.04009 * l7z * s12 * c2 + 0.14454 * l7z * s2 * c12 - 0.04009 * l7z * c12 * c2 + 1.0 * l8y * q3 * s12 * s2 + 1.0 * l8y * q3 * s2 * c12 - 0.4318 * l8y * s12 * s2 - 0.4318 * l8y * s2 * c12 - 1.0 * l9y * q3 * s12 * s2 - 1.0 * l9y * q3 * s2 * c12 + 0.04009 * l9y * s12 * c2 + 0.04009 * l9y * c12 * c2;
	M(0,2) = -1.0 * l8y * s12 * c2 - 1.0 * l8y * c12 * c2 + 1.0 * l9y * s12 * c2 + 1.0 * l9y * c12 * c2;
	M(1,0) = L2xz * c2 - L2yz * s2 - 1.0 * L6xz * c2 + 1.0 * L6yz * s2 + 1.0 * L7xz * c2 - 1.0 * L7yz * s2 + 1.0 * L8xy * c2 - 1.0 * L8yz * s2 - 1.0 * L9xy * c2 + 1.0 * L9yz * s2 + 0.14454 * l4z * s12 * s2 - 0.04009 * l4z * s12 * c2 + 0.14454 * l4z * s2 * c12 - 0.04009 * l4z * c12 * c2 + 0.18262 * l5z * s12 * s2 - 0.04009 * l5z * s12 * c2 + 0.18262 * l5z * s2 * c12 - 0.04009 * l5z * c12 * c2 + 0.14454 * l6z * s12 * s2 - 0.04009 * l6z * s12 * c2 + 0.14454 * l6z * s2 * c12 - 0.04009 * l6z * c12 * c2 + 0.14454 * l7z * s12 * s2 - 0.04009 * l7z * s12 * c2 + 0.14454 * l7z * s2 * c12 - 0.04009 * l7z * c12 * c2 + 1.0 * l8y * q3 * s12 * s2 + 1.0 * l8y * q3 * s2 * c12 - 0.4318 * l8y * s12 * s2 - 0.4318 * l8y * s2 * c12 - 1.0 * l9y * q3 * s12 * s2 - 1.0 * l9y * q3 * s2 * c12 + 0.04009 * l9y * s12 * c2 + 0.04009 * l9y * c12 * c2;
	M(1,1) = L2zz + 1.0 * L6zz + 1.0 * L7zz + 1.0 * L8yy + 1.0 * L9yy + 0.28908 * l6x * s12 * s2 * c2 - 0.08018 * l6x * s12 * c22 - 0.08018 * l6x * s22 + 0.28908 * l6x * s2 * c12 * c2 - 0.28908 * l6x * s2 * c2 - 0.08018 * l6x * c12 * c22 - 0.28908 * l6y * s12 * s22 + 0.08018 * l6y * s12 * s2 * c2 - 0.28908 * l6y * s22 * c12 + 0.08018 * l6y * s2 * c12 * c2 - 0.08018 * l6y * s2 * c2 - 0.28908 * l6y * c22 - 0.28908 * l7x * s12 * s2 * c2 + 0.08018 * l7x * s12 * c22 + 0.08018 * l7x * s22 - 0.28908 * l7x * s2 * c12 * c2 + 0.28908 * l7x * s2 * c2 + 0.08018 * l7x * c12 * c22 + 0.28908 * l7y * s12 * s22 - 0.08018 * l7y * s12 * s2 * c2 + 0.28908 * l7y * s22 * c12 - 0.08018 * l7y * s2 * c12 * c2 + 0.08018 * l7y * s2 * c2 + 0.28908 * l7y * c22 - 2.0 * l8x * q3 * s12 * s2 * c2 - 2.0 * l8x * q3 * s2 * c12 * c2 + 2.0 * l8x * q3 * s2 * c2 + 0.8636 * l8x * s12 * s2 * c2 + 0.8636 * l8x * s2 * c12 * c2 - 0.8636 * l8x * s2 * c2 + 2.0 * l8z * q3 * s12 * s22 + 2.0 * l8z * q3 * s22 * c12 + 2.0 * l8z * q3 * c22 - 0.8636 * l8z * s12 * s22 - 0.8636 * l8z * s22 * c12 - 0.8636 * l8z * c22 - 2.0 * l9x * q3 * s12 * s2 * c2 - 2.0 * l9x * q3 * s2 * c12 * c2 + 2.0 * l9x * q3 * s2 * c2 + 0.08018 * l9x * s12 * c22 + 0.08018 * l9x * s22 + 0.08018 * l9x * c12 * c22 + 2.0 * l9z * q3 * s12 * s22 + 2.0 * l9z * q3 * s22 * c12 + 2.0 * l9z * q3 * c22 - 0.08018 * l9z * s12 * s2 * c2 - 0.08018 * l9z * s2 * c12 * c2 + 0.08018 * l9z * s2 * c2 + 0.0208918116 * m4 * s12 * s22 - 0.0115892172 * m4 * s12 * s2 * c2 + 0.0016072081 * m4 * s12 * c22 + 0.0208918116 * m4 * s22 * c12 + 0.0016072081 * m4 * s22 - 0.0115892172 * m4 * s2 * c12 * c2 + 0.0115892172 * m4 * s2 * c2 + 0.0016072081 * m4 * c12 * c22 + 0.0208918116 * m4 * c22 + 0.0333500644 * m5 * s12 * s22 - 0.0146424716 * m5 * s12 * s2 * c2 + 0.0016072081 * m5 * s12 * c22 + 0.0333500644 * m5 * s22 * c12 + 0.0016072081 * m5 * s22 - 0.0146424716 * m5 * s2 * c12 * c2 + 0.0146424716 * m5 * s2 * c2 + 0.0016072081 * m5 * c12 * c22 + 0.0333500644 * m5 * c22 + 0.0208918116 * m6 * s12 * s22 - 0.0115892172 * m6 * s12 * s2 * c2 + 0.0016072081 * m6 * s12 * c22 + 0.0208918116 * m6 * s22 * c12 + 0.0016072081 * m6 * s22 - 0.0115892172 * m6 * s2 * c12 * c2 + 0.0115892172 * m6 * s2 * c2 + 0.0016072081 * m6 * c12 * c22 + 0.0208918116 * m6 * c22 + 0.0208918116 * m7 * s12 * s22 - 0.0115892172 * m7 * s12 * s2 * c2 + 0.0016072081 * m7 * s12 * c22 + 0.0208918116 * m7 * s22 * c12 + 0.0016072081 * m7 * s22 - 0.0115892172 * m7 * s2 * c12 * c2 + 0.0115892172 * m7 * s2 * c2 + 0.0016072081 * m7 * c12 * c22 + 0.0208918116 * m7 * c22 + 1.0 * m8 * q32 * s12 * s22 + 1.0 * m8 * q32 * s22 * c12 + 1.0 * m8 * q32 * c22 - 0.8636 * m8 * q3 * s12 * s22 - 0.8636 * m8 * q3 * s22 * c12 - 0.8636 * m8 * q3 * c22 + 0.18645124 * m8 * s12 * s22 + 0.18645124 * m8 * s22 * c12 + 0.18645124 * m8 * c22 + 1.0 * m9 * q32 * s12 * s22 + 1.0 * m9 * q32 * s22 * c12 + 1.0 * m9 * q32 * c22 - 0.08018 * m9 * q3 * s12 * s2 * c2 - 0.08018 * m9 * q3 * s2 * c12 * c2 + 0.08018 * m9 * q3 * s2 * c2 + 0.0016072081 * m9 * s12 * c22 + 0.0016072081 * m9 * s22 + 0.0016072081 * m9 * c12 * c22;
	M(1,2) = 1.0 * l8x * s12 * c22 + 1.0 * l8x * s22 + 1.0 * l8x * c12 * c22 - 1.0 * l8z * s12 * s2 * c2 - 1.0 * l8z * s2 * c12 * c2 + 1.0 * l8z * s2 * c2 + 1.0 * l9x * s12 * c22 + 1.0 * l9x * s22 + 1.0 * l9x * c12 * c22 - 1.0 * l9z * s12 * s2 * c2 - 1.0 * l9z * s2 * c12 * c2 + 1.0 * l9z * s2 * c2 - 1.0 * m8 * q3 * s12 * s2 * c2 - 1.0 * m8 * q3 * s2 * c12 * c2 + 1.0 * m8 * q3 * s2 * c2 + 0.4318 * m8 * s12 * s2 * c2 + 0.4318 * m8 * s2 * c12 * c2 - 0.4318 * m8 * s2 * c2 - 1.0 * m9 * q3 * s12 * s2 * c2 - 1.0 * m9 * q3 * s2 * c12 * c2 + 1.0 * m9 * q3 * s2 * c2 + 0.04009 * m9 * s12 * c22 + 0.04009 * m9 * s22 + 0.04009 * m9 * c12 * c22;
	M(2,0) = -1.0 * l8y * s12 * c2 - 1.0 * l8y * c12 * c2 + 1.0 * l9y * s12 * c2 + 1.0 * l9y * c12 * c2;
	M(2,1) = 1.0 * l8x * s12 * c22 + 1.0 * l8x * s22 + 1.0 * l8x * c12 * c22 - 1.0 * l8z * s12 * s2 * c2 - 1.0 * l8z * s2 * c12 * c2 + 1.0 * l8z * s2 * c2 + 1.0 * l9x * s12 * c22 + 1.0 * l9x * s22 + 1.0 * l9x * c12 * c22 - 1.0 * l9z * s12 * s2 * c2 - 1.0 * l9z * s2 * c12 * c2 + 1.0 * l9z * s2 * c2 - 1.0 * m8 * q3 * s12 * s2 * c2 - 1.0 * m8 * q3 * s2 * c12 * c2 + 1.0 * m8 * q3 * s2 * c2 + 0.4318 * m8 * s12 * s2 * c2 + 0.4318 * m8 * s2 * c12 * c2 - 0.4318 * m8 * s2 * c2 - 1.0 * m9 * q3 * s12 * s2 * c2 - 1.0 * m9 * q3 * s2 * c12 * c2 + 1.0 * m9 * q3 * s2 * c2 + 0.04009 * m9 * s12 * c22 + 0.04009 * m9 * s22 + 0.04009 * m9 * c12 * c22;
	M(2,2) = 1.0 * m8 * s12 * c22 + 1.0 * m8 * s22 + 1.0 * m8 * c12 * c22 + 1.0 * m9 * s12 * c22 + 1.0 * m9 * s22 + 1.0 * m9 * c12 * c22;
	M(3,3) = Ia10;
	M(4,4) = 1.0186 * Ia11;
	M(5,5) = Ia14;
	M(6,6) = Ia15;

	return M;
}

/**
* @brief Gravity vector function
* Compute the gravity vector of the DVRK PSM arm dynamic model
* @param [out] g: the gravity of the PSM
* @param [in] q: the current joint positions of the PSM
* @param [in] withCounterWeight: if the counterweight is modeled and accounted in the computation or not
* @param [in] withFriction: if the joint friction is modeled and accounted in the computation or not
*
*/
psmActiveJointsVectorf computeDVRKPSMGravityVector(const psmActiveJointsVectorf& q, const bool& withCounterWeight, bool withFriction){

	psmActiveJointsVectorf g;
	g.setZero();

	double m9 = 0.0, l9x = 0.0, l9y = 0.0, l9z = 0.0, L9xx = 0.0, L9xy = 0.0, L9xz = 0.0, L9yy = 0.0, L9yz = 0.0, L9zz = 0.0;
	double Fo1 = 0.0, Fo2 = 0.0, Fo8 = 0.0, Fo10 = 0.0, Fo11 = 0.0, Fo12 = 0.0, Fo13 = 0.0, Fo14 = 0.0, Fo15 = 0.0, Fo16 = 0.0;

	double q1 = q(0);
	double q2 = q(1);
	double q3 = q(2);
	double q4 = q(3);

	double s1 = sin(q1);
	double c1 = cos(q1);
	double s2 = sin(q2);
	double c2 = cos(q2);
	double s22 = pow(s2, 2);
	double s23 = pow(s2, 3);
	double c22 = pow(c2, 2);
	double c23 = pow(c2, 3);

	double K10 = 0.0;
	//double K10 = 0.0001873;

	if (withCounterWeight) {
		L9xx = 0.0036083;
		L9xy = 0.000204;
		L9xz = 0.0009268;
		L9yy = 0.0037699;
		L9yz = -0.0006548;
		L9zz = 0.0004368;
		m9 = 0.3499395;
		l9x = 0.0091043;
		l9y = -0.0064588;
		l9z = -0.0348129;
	}

	if (withFriction) {
	
		Fo1 = -0.3000027;
		Fo2 = 0.1419199;
		Fo8 = 0.0038872;
		Fo10 = -0.0048904;
		Fo11 = 0.1329238;
		Fo12 = 0.0777171;
		Fo13 = 0.0561429;
		Fo14 = -0.1661489;
		Fo15 = -0.0169157;
		Fo16 = -0.0458516;

	}

	g(0) = Fo1 - 9.81 * l1x * s1 - 9.81 * l1y * c1 - 9.81 * l2x * s1 * s2 - 9.81 * l2y * s1 * c2 - 9.81 * l2z * c1 + 9.81 * l4y * s1 * s22 + 9.81 * l4y * s1 * c22 - 9.81 * l4z * c1 + 9.81 * l5y * s1 * s22 + 9.81 * l5y * s1 * c22 - 9.81 * l5z * c1 + 9.81 * l6x * s1 * s23 + 9.81 * l6x * s1 * s2 * c22 + 9.81 * l6y * s1 * s22 * c2 + 9.81 * l6y * s1 * c23 - 9.81 * l6z * c1 - 9.81 * l7x * s1 * s23 - 9.81 * l7x * s1 * s2 * c22 - 9.81 * l7y * s1 * s22 * c2 - 9.81 * l7y * s1 * c23 - 9.81 * l7z * c1 + 9.81 * l8x * s1 * s23 + 9.81 * l8x * s1 * s2 * c22 + 9.81 * l8y * c1 + 9.81 * l8z * s1 * s22 * c2 + 9.81 * l8z * s1 * c23 - 9.81 * l9x * s1 * s2 + 9.81 * l9y * c1 - 9.81 * l9z * s1 * c2 - 0.3932829 * m4 * s1 * s2 - 1.4179374 * m4 * s1 * c2 - 0.3932829 * m5 * s1 * s2 - 1.7915022 * m5 * s1 * c2 - 0.3932829 * m6 * s1 * s2 - 1.4179374 * m6 * s1 * c2 - 0.3932829 * m7 * s1 * s2 - 1.4179374 * m7 * s1 * c2 + 9.81 * m8 * q3 * s1 * s22 * c2 + 9.81 * m8 * q3 * s1 * c23 + 0.3932829 * m8 * s1 * s23 - 2.8180206 * m8 * s1 * s22 * c2 + 0.3932829 * m8 * s1 * s2 * c22 - 0.3932829 * m8 * s1 * s2 - 2.8180206 * m8 * s1 * c23 - 1.4179374 * m8 * s1 * c2 - 9.81 * m9 * q3 * s1 * c2 - 0.3932829 * m9 * s1 * s2;
	g(1) = Fo2 + 9.81 * l2x * c1 * c2 - 9.81 * l2y * s2 * c1 - 9.81 * l6x * s22 * c1 * c2 - 9.81 * l6x * c1 * c23 + 9.81 * l6y * s23 * c1 + 9.81 * l6y * s2 * c1 * c22 + 9.81 * l7x * s22 * c1 * c2 + 9.81 * l7x * c1 * c23 - 9.81 * l7y * s23 * c1 - 9.81 * l7y * s2 * c1 * c22 - 9.81 * l8x * s22 * c1 * c2 - 9.81 * l8x * c1 * c23 + 9.81 * l8z * s23 * c1 + 9.81 * l8z * s2 * c1 * c22 + 9.81 * l9x * c1 * c2 - 9.81 * l9z * s2 * c1 - 1.4179374 * m4 * s2 * c1 + 0.3932829 * m4 * c1 * c2 - 1.7915022 * m5 * s2 * c1 + 0.3932829 * m5 * c1 * c2 - 1.4179374 * m6 * s2 * c1 + 0.3932829 * m6 * c1 * c2 - 1.4179374 * m7 * s2 * c1 + 0.3932829 * m7 * c1 * c2 + 9.81 * m8 * q3 * s23 * c1 + 9.81 * m8 * q3 * s2 * c1 * c22 - 2.8180206 * m8 * s23 * c1 - 0.3932829 * m8 * s22 * c1 * c2 - 2.8180206 * m8 * s2 * c1 * c22 - 1.4179374 * m8 * s2 * c1 - 0.3932829 * m8 * c1 * c23 + 0.3932829 * m8 * c1 * c2 - 9.81 * m9 * q3 * s2 * c1 + 0.3932829 * m9 * c1 * c2;
	g(2) = Fo8 - 9.81 * m8 * s22 * c1 * c2 - 9.81 * m8 * c1 * pow(c2,3) + 9.81 * m9 * c1 * c2;
	g(3) = Fo10 + K10 * q4;
	g(4) = 1.0186 * Fo11 - 0.8306 * Fo12 - 0.8306 * Fo13;
	g(5) = 1.2178 * Fo12 + Fo14 - 1.2178 * Fo16;
	g(6) = 1.2178 * Fo13 + Fo15 + 1.2178 * Fo16;

	return g;

}

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
psmActiveJointsVectorf computeDVRKPSMCoriolisVector(const psmActiveJointsVectorf& q, const psmActiveJointsVectorf& dq, const bool& withCounterWeight, bool withFriction){

	psmActiveJointsVectorf C;
	C.setZero();

	double m9 = 0.0, l9x = 0.0, l9y = 0.0, l9z = 0.0, L9xx = 0.0, L9xy = 0.0, L9xz = 0.0, L9yy = 0.0, L9yz = 0.0, L9zz = 0.0;
	double Fc1 = 0.0, Fc2 = 0.0, Fc8 = 0.0, Fc10 = 0.0, Fc11 = 0.0, Fc12 = 0.0, Fc13 = 0.0, Fc14 = 0.0, Fc15 = 0.0, Fc16 = 0.0;
	double Fv1 = 0.0, Fv2 = 0.0, Fv8 = 0.0, Fv10 = 0.0, Fv11 = 0.0, Fv12 = 0.0, Fv13 = 0.0, Fv14 = 0.0, Fv15 = 0.0, Fv16 = 0.0;

	if (withCounterWeight) {
		L9xx = 0.0036083;
		L9xy = 0.000204;
		L9xz = 0.0009268;
		L9yy = 0.0037699;
		L9yz = -0.0006548;
		L9zz = 0.0004368;
		m9 = 0.3499395;
		l9x = 0.0091043;
		l9y = -0.0064588;
		l9z = -0.0348129;
	}

	if (withFriction) {

		Fc1 = 0.05003;
		Fc2 = 0.0668202;
		Fc8 = 0.4475953;
		Fc10 = 0.002465;
		Fc11 = 0.0069935;
		Fc12 = 0.0007803;
		Fc13 = 0.0021698;
		Fc14 = 0.0062966;
		Fc15 = 0.0061915;
		Fc16 = 0.0001049;

		Fv1 = 0.0905326;
		Fv2 = 0.2180617;
		Fv8 = 1.2035645;
		Fv10 = 0.0011114;
		Fv11 = 0.012171;
		Fv12 = 0.0020153;
		Fv13 = 0.0022465;
		Fv14 = 0.006165;
		Fv15 = 0.0073884;
		Fv16 = 0.0034451;

	}

	double q1 = q(0);
	double q2 = q(1);
	double q3 = q(2);
	double q4 = q(3);

	double dq1 = dq(0);
	double dq2 = dq(1);
	double dq3 = dq(2);
	double dq4 = dq(3);
	double dq5 = dq(4);
	double dq6 = dq(5);
	double dq7 = dq(6);

	double s1 = sin(q1);
	double c1 = cos(q1);
	double s2 = sin(q2);
	double c2 = cos(q2);
	double c12 = pow(c1, 2);
	double s12 = pow(s1, 2);
	double s22 = pow(s2, 2);
	double c22 = pow(c2, 2);
	double s23 = pow(s2, 3);
	double c23 = pow(c2, 3);
	double q32 = pow(q3, 2);
	double dq12 = pow(dq1, 2);
	double dq22 = pow(dq2, 2);


	C(0) = Fc1 * sign(dq1) + Fv1 * dq1 - 2 * L2xx * dq1 * dq2 * s2 * c2 + 2 * L2xy * dq1 * dq2 * s22 - 2 * L2xy * dq1 * dq2 * c22 - L2xz * dq22 * s2 + 2 * L2yy * dq1 * dq2 * s2 * c2 - L2yz * dq22 * c2 - 2.0 * L6xx * dq1 * dq2 * s2 * c2 + 2.0 * L6xy * dq1 * dq2 * s22 - 2.0 * L6xy * dq1 * dq2 * c22 + 1.0 * L6xz * dq22 * s2 + 2.0 * L6yy * dq1 * dq2 * s2 * c2 + 1.0 * L6yz * dq22 * c2 - 2.0 * L7xx * dq1 * dq2 * s2 * c2 + 2.0 * L7xy * dq1 * dq2 * s22 - 2.0 * L7xy * dq1 * dq2 * c22 - 1.0 * L7xz * dq22 * s2 + 2.0 * L7yy * dq1 * dq2 * s2 * c2 - 1.0 * L7yz * dq22 * c2 - 2.0 * L8xx * dq1 * dq2 * s2 * c2 - 1.0 * L8xy * dq22 * s2 + 2.0 * L8xz * dq1 * dq2 * s22 - 2.0 * L8xz * dq1 * dq2 * c22 - 1.0 * L8yz * dq22 * c2 + 2.0 * L8zz * dq1 * dq2 * s2 * c2 - 2.0 * L9xx * dq1 * dq2 * s2 * c2 + 1.0 * L9xy * dq22 * s2 + 2.0 * L9xz * dq1 * dq2 * s22 - 2.0 * L9xz * dq1 * dq2 * c22 + 1.0 * L9yz * dq22 * c2 + 2.0 * L9zz * dq1 * dq2 * s2 * c2 + 0.28908 * dq1 * dq2 * l4y * s12 * s2 - 0.08018 * dq1 * dq2 * l4y * s12 * c2 + 0.28908 * dq1 * dq2 * l4y * s2 * c12 - 0.08018 * dq1 * dq2 * l4y * c12 * c2 + 0.36524 * dq1 * dq2 * l5y * s12 * s2 - 0.08018 * dq1 * dq2 * l5y * s12 * c2 + 0.36524 * dq1 * dq2 * l5y * s2 * c12 - 0.08018 * dq1 * dq2 * l5y * c12 * c2 + 0.28908 * dq1 * dq2 * l6x * s12 * s22 - 0.16036 * dq1 * dq2 * l6x * s12 * s2 * c2 - 0.28908 * dq1 * dq2 * l6x * s12 * c22 + 0.28908 * dq1 * dq2 * l6x * s22 * c12 - 0.16036 * dq1 * dq2 * l6x * s2 * c12 * c2 - 0.28908 * dq1 * dq2 * l6x * c12 * c22 + 0.08018 * dq1 * dq2 * l6y * s12 * s22 + 0.57816 * dq1 * dq2 * l6y * s12 * s2 * c2 - 0.08018 * dq1 * dq2 * l6y * s12 * c22 + 0.08018 * dq1 * dq2 * l6y * s22 * c12 + 0.57816 * dq1 * dq2 * l6y * s2 * c12 * c2 - 0.08018 * dq1 * dq2 * l6y * c12 * c22 - 0.28908 * dq1 * dq2 * l7x * s12 * s22 + 0.16036 * dq1 * dq2 * l7x * s12 * s2 * c2 + 0.28908 * dq1 * dq2 * l7x * s12 * c22 - 0.28908 * dq1 * dq2 * l7x * s22 * c12 + 0.16036 * dq1 * dq2 * l7x * s2 * c12 * c2 + 0.28908 * dq1 * dq2 * l7x * c12 * c22 - 0.08018 * dq1 * dq2 * l7y * s12 * s22 - 0.57816 * dq1 * dq2 * l7y * s12 * s2 * c2 + 0.08018 * dq1 * dq2 * l7y * s12 * c22 - 0.08018 * dq1 * dq2 * l7y * s22 * c12 - 0.57816 * dq1 * dq2 * l7y * s2 * c12 * c2 + 0.08018 * dq1 * dq2 * l7y * c12 * c22 - 2.0 * dq1 * dq2 * l8x * q3 * s12 * s22 + 2.0 * dq1 * dq2 * l8x * q3 * s12 * c22 - 2.0 * dq1 * dq2 * l8x * q3 * s22 * c12 + 2.0 * dq1 * dq2 * l8x * q3 * c12 * c22 + 0.8636 * dq1 * dq2 * l8x * s12 * s22 - 0.8636 * dq1 * dq2 * l8x * s12 * c22 + 0.8636 * dq1 * dq2 * l8x * s22 * c12 - 0.8636 * dq1 * dq2 * l8x * c12 * c22 - 4.0 * dq1 * dq2 * l8z * q3 * s12 * s2 * c2 - 4.0 * dq1 * dq2 * l8z * q3 * s2 * c12 * c2 + 1.7272 * dq1 * dq2 * l8z * s12 * s2 * c2 + 1.7272 * dq1 * dq2 * l8z * s2 * c12 * c2 - 2.0 * dq1 * dq2 * l9x * q3 * s12 * s22 + 2.0 * dq1 * dq2 * l9x * q3 * s12 * c22 - 2.0 * dq1 * dq2 * l9x * q3 * s22 * c12 + 2.0 * dq1 * dq2 * l9x * q3 * c12 * c22 + 0.16036 * dq1 * dq2 * l9x * s12 * s2 * c2 + 0.16036 * dq1 * dq2 * l9x * s2 * c12 * c2 - 4.0 * dq1 * dq2 * l9z * q3 * s12 * s2 * c2 - 4.0 * dq1 * dq2 * l9z * q3 * s2 * c12 * c2 - 0.08018 * dq1 * dq2 * l9z * s12 * s22 + 0.08018 * dq1 * dq2 * l9z * s12 * c22 - 0.08018 * dq1 * dq2 * l9z * s22 * c12 + 0.08018 * dq1 * dq2 * l9z * c12 * c22 - 0.0115892172 * dq1 * dq2 * m4 * s12 * s22 - 0.038569207 * dq1 * dq2 * m4 * s12 * s2 * c2 + 0.0115892172 * dq1 * dq2 * m4 * s12 * c22 - 0.0115892172 * dq1 * dq2 * m4 * s22 * c12 - 0.038569207 * dq1 * dq2 * m4 * s2 * c12 * c2 + 0.0115892172 * dq1 * dq2 * m4 * c12 * c22 - 0.0146424716 * dq1 * dq2 * m5 * s12 * s22 - 0.0634857126 * dq1 * dq2 * m5 * s12 * s2 * c2 + 0.0146424716 * dq1 * dq2 * m5 * s12 * c22 - 0.0146424716 * dq1 * dq2 * m5 * s22 * c12 - 0.0634857126 * dq1 * dq2 * m5 * s2 * c12 * c2 + 0.0146424716 * dq1 * dq2 * m5 * c12 * c22 - 0.0115892172 * dq1 * dq2 * m6 * s12 * s22 - 0.038569207 * dq1 * dq2 * m6 * s12 * s2 * c2 + 0.0115892172 * dq1 * dq2 * m6 * s12 * c22 - 0.0115892172 * dq1 * dq2 * m6 * s22 * c12 - 0.038569207 * dq1 * dq2 * m6 * s2 * c12 * c2 + 0.0115892172 * dq1 * dq2 * m6 * c12 * c22 - 0.0115892172 * dq1 * dq2 * m7 * s12 * s22 - 0.038569207 * dq1 * dq2 * m7 * s12 * s2 * c2 + 0.0115892172 * dq1 * dq2 * m7 * s12 * c22 - 0.0115892172 * dq1 * dq2 * m7 * s22 * c12 - 0.038569207 * dq1 * dq2 * m7 * s2 * c12 * c2 + 0.0115892172 * dq1 * dq2 * m7 * c12 * c22 - 2.0 * dq1 * dq2 * m8 * q32 * s12 * s2 * c2 - 2.0 * dq1 * dq2 * m8 * q32 * s2 * c12 * c2 + 1.7272 * dq1 * dq2 * m8 * q3 * s12 * s2 * c2 + 1.7272 * dq1 * dq2 * m8 * q3 * s2 * c12 * c2 - 0.37290248 * dq1 * dq2 * m8 * s12 * s2 * c2 - 0.37290248 * dq1 * dq2 * m8 * s2 * c12 * c2 - 2.0 * dq1 * dq2 * m9 * q32 * s12 * s2 * c2 - 2.0 * dq1 * dq2 * m9 * q32 * s2 * c12 * c2 - 0.08018 * dq1 * dq2 * m9 * q3 * s12 * s22 + 0.08018 * dq1 * dq2 * m9 * q3 * s12 * c22 - 0.08018 * dq1 * dq2 * m9 * q3 * s22 * c12 + 0.08018 * dq1 * dq2 * m9 * q3 * c12 * c22 + 0.0032144162 * dq1 * dq2 * m9 * s12 * s2 * c2 + 0.0032144162 * dq1 * dq2 * m9 * s2 * c12 * c2 + 2.0 * dq1 * dq3 * l8x * s12 * s2 * c2 + 2.0 * dq1 * dq3 * l8x * s2 * c12 * c2 + 2.0 * dq1 * dq3 * l8z * s12 * c22 + 2.0 * dq1 * dq3 * l8z * c12 * c22 + 2.0 * dq1 * dq3 * l9x * s12 * s2 * c2 + 2.0 * dq1 * dq3 * l9x * s2 * c12 * c2 + 2.0 * dq1 * dq3 * l9z * s12 * c22 + 2.0 * dq1 * dq3 * l9z * c12 * c22 + 2.0 * dq1 * dq3 * m8 * q3 * s12 * c22 + 2.0 * dq1 * dq3 * m8 * q3 * c12 * c22 - 0.8636 * dq1 * dq3 * m8 * s12 * c22 - 0.8636 * dq1 * dq3 * m8 * c12 * c22 + 2.0 * dq1 * dq3 * m9 * q3 * s12 * c22 + 2.0 * dq1 * dq3 * m9 * q3 * c12 * c22 + 0.08018 * dq1 * dq3 * m9 * s12 * s2 * c2 + 0.08018 * dq1 * dq3 * m9 * s2 * c12 * c2 + 0.04009 * dq22 * l4z * s12 * s2 + 0.14454 * dq22 * l4z * s12 * c2 + 0.04009 * dq22 * l4z * s2 * c12 + 0.14454 * dq22 * l4z * c12 * c2 + 0.04009 * dq22 * l5z * s12 * s2 + 0.18262 * dq22 * l5z * s12 * c2 + 0.04009 * dq22 * l5z * s2 * c12 + 0.18262 * dq22 * l5z * c12 * c2 + 0.04009 * dq22 * l6z * s12 * s2 + 0.14454 * dq22 * l6z * s12 * c2 + 0.04009 * dq22 * l6z * s2 * c12 + 0.14454 * dq22 * l6z * c12 * c2 + 0.04009 * dq22 * l7z * s12 * s2 + 0.14454 * dq22 * l7z * s12 * c2 + 0.04009 * dq22 * l7z * s2 * c12 + 0.14454 * dq22 * l7z * c12 * c2 + 1.0 * dq22 * l8y * q3 * s12 * c2 + 1.0 * dq22 * l8y * q3 * c12 * c2 - 0.4318 * dq22 * l8y * s12 * c2 - 0.4318 * dq22 * l8y * c12 * c2 - 1.0 * dq22 * l9y * q3 * s12 * c2 - 1.0 * dq22 * l9y * q3 * c12 * c2 - 0.04009 * dq22 * l9y * s12 * s2 - 0.04009 * dq22 * l9y * s2 * c12 + 2.0 * dq2 * dq3 * l8y * s12 * s2 + 2.0 * dq2 * dq3 * l8y * s2 * c12 - 2.0 * dq2 * dq3 * l9y * s12 * s2 - 2.0 * dq2 * dq3 * l9y * s2 * c12;
	C(1) = Fc2 * sign(dq2) + Fv2 * dq2 + L2xx * dq12 * s2 * c2 - L2xy * dq12 * s22 + L2xy * dq12 * c22 - L2yy * dq12 * s2 * c2 + 1.0 * L6xx * dq12 * s2 * c2 - 1.0 * L6xy * dq12 * s22 + 1.0 * L6xy * dq12 * c22 - 1.0 * L6yy * dq12 * s2 * c2 + 1.0 * L7xx * dq12 * s2 * c2 - 1.0 * L7xy * dq12 * s22 + 1.0 * L7xy * dq12 * c22 - 1.0 * L7yy * dq12 * s2 * c2 + 1.0 * L8xx * dq12 * s2 * c2 - 1.0 * L8xz * dq12 * s22 + 1.0 * L8xz * dq12 * c22 - 1.0 * L8zz * dq12 * s2 * c2 + 1.0 * L9xx * dq12 * s2 * c2 - 1.0 * L9xz * dq12 * s22 + 1.0 * L9xz * dq12 * c22 - 1.0 * L9zz * dq12 * s2 * c2 - 0.14454 * dq12 * l4y * s12 * s2 + 0.04009 * dq12 * l4y * s12 * c2 - 0.14454 * dq12 * l4y * s2 * c12 + 0.04009 * dq12 * l4y * c12 * c2 - 0.18262 * dq12 * l5y * s12 * s2 + 0.04009 * dq12 * l5y * s12 * c2 - 0.18262 * dq12 * l5y * s2 * c12 + 0.04009 * dq12 * l5y * c12 * c2 - 0.14454 * dq12 * l6x * s12 * s22 + 0.08018 * dq12 * l6x * s12 * s2 * c2 + 0.14454 * dq12 * l6x * s12 * c22 - 0.14454 * dq12 * l6x * s22 * c12 + 0.08018 * dq12 * l6x * s2 * c12 * c2 + 0.14454 * dq12 * l6x * c12 * c22 - 0.04009 * dq12 * l6y * s12 * s22 - 0.28908 * dq12 * l6y * s12 * s2 * c2 + 0.04009 * dq12 * l6y * s12 * c22 - 0.04009 * dq12 * l6y * s22 * c12 - 0.28908 * dq12 * l6y * s2 * c12 * c2 + 0.04009 * dq12 * l6y * c12 * c22 + 0.14454 * dq12 * l7x * s12 * s22 - 0.08018 * dq12 * l7x * s12 * s2 * c2 - 0.14454 * dq12 * l7x * s12 * c22 + 0.14454 * dq12 * l7x * s22 * c12 - 0.08018 * dq12 * l7x * s2 * c12 * c2 - 0.14454 * dq12 * l7x * c12 * c22 + 0.04009 * dq12 * l7y * s12 * s22 + 0.28908 * dq12 * l7y * s12 * s2 * c2 - 0.04009 * dq12 * l7y * s12 * c22 + 0.04009 * dq12 * l7y * s22 * c12 + 0.28908 * dq12 * l7y * s2 * c12 * c2 - 0.04009 * dq12 * l7y * c12 * c22 + 1.0 * dq12 * l8x * q3 * s12 * s22 - 1.0 * dq12 * l8x * q3 * s12 * c22 + 1.0 * dq12 * l8x * q3 * s22 * c12 - 1.0 * dq12 * l8x * q3 * c12 * c22 - 0.4318 * dq12 * l8x * s12 * s22 + 0.4318 * dq12 * l8x * s12 * c22 - 0.4318 * dq12 * l8x * s22 * c12 + 0.4318 * dq12 * l8x * c12 * c22 + 2.0 * dq12 * l8z * q3 * s12 * s2 * c2 + 2.0 * dq12 * l8z * q3 * s2 * c12 * c2 - 0.8636 * dq12 * l8z * s12 * s2 * c2 - 0.8636 * dq12 * l8z * s2 * c12 * c2 + 1.0 * dq12 * l9x * q3 * s12 * s22 - 1.0 * dq12 * l9x * q3 * s12 * c22 + 1.0 * dq12 * l9x * q3 * s22 * c12 - 1.0 * dq12 * l9x * q3 * c12 * c22 - 0.08018 * dq12 * l9x * s12 * s2 * c2 - 0.08018 * dq12 * l9x * s2 * c12 * c2 + 2.0 * dq12 * l9z * q3 * s12 * s2 * c2 + 2.0 * dq12 * l9z * q3 * s2 * c12 * c2 + 0.04009 * dq12 * l9z * s12 * s22 - 0.04009 * dq12 * l9z * s12 * c22 + 0.04009 * dq12 * l9z * s22 * c12 - 0.04009 * dq12 * l9z * c12 * c22 + 0.0057946086 * dq12 * m4 * s12 * s22 + 0.0192846035 * dq12 * m4 * s12 * s2 * c2 - 0.0057946086 * dq12 * m4 * s12 * c22 + 0.0057946086 * dq12 * m4 * s22 * c12 + 0.0192846035 * dq12 * m4 * s2 * c12 * c2 - 0.0057946086 * dq12 * m4 * c12 * c22 + 0.0073212358 * dq12 * m5 * s12 * s22 + 0.0317428563 * dq12 * m5 * s12 * s2 * c2 - 0.0073212358 * dq12 * m5 * s12 * c22 + 0.0073212358 * dq12 * m5 * s22 * c12 + 0.0317428563 * dq12 * m5 * s2 * c12 * c2 - 0.0073212358 * dq12 * m5 * c12 * c22 + 0.0057946086 * dq12 * m6 * s12 * s22 + 0.0192846035 * dq12 * m6 * s12 * s2 * c2 - 0.0057946086 * dq12 * m6 * s12 * c22 + 0.0057946086 * dq12 * m6 * s22 * c12 + 0.0192846035 * dq12 * m6 * s2 * c12 * c2 - 0.0057946086 * dq12 * m6 * c12 * c22 + 0.0057946086 * dq12 * m7 * s12 * s22 + 0.0192846035 * dq12 * m7 * s12 * s2 * c2 - 0.0057946086 * dq12 * m7 * s12 * c22 + 0.0057946086 * dq12 * m7 * s22 * c12 + 0.0192846035 * dq12 * m7 * s2 * c12 * c2 - 0.0057946086 * dq12 * m7 * c12 * c22 + 1.0 * dq12 * m8 * q32 * s12 * s2 * c2 + 1.0 * dq12 * m8 * q32 * s2 * c12 * c2 - 0.8636 * dq12 * m8 * q3 * s12 * s2 * c2 - 0.8636 * dq12 * m8 * q3 * s2 * c12 * c2 + 0.18645124 * dq12 * m8 * s12 * s2 * c2 + 0.18645124 * dq12 * m8 * s2 * c12 * c2 + 1.0 * dq12 * m9 * q32 * s12 * s2 * c2 + 1.0 * dq12 * m9 * q32 * s2 * c12 * c2 + 0.04009 * dq12 * m9 * q3 * s12 * s22 - 0.04009 * dq12 * m9 * q3 * s12 * c22 + 0.04009 * dq12 * m9 * q3 * s22 * c12 - 0.04009 * dq12 * m9 * q3 * c12 * c22 - 0.0016072081 * dq12 * m9 * s12 * s2 * c2 - 0.0016072081 * dq12 * m9 * s2 * c12 * c2 - 0.14454 * dq22 * l6x * s12 * s22 + 0.08018 * dq22 * l6x * s12 * s2 * c2 + 0.14454 * dq22 * l6x * s12 * c22 - 0.14454 * dq22 * l6x * s22 * c12 + 0.14454 * dq22 * l6x * s22 + 0.08018 * dq22 * l6x * s2 * c12 * c2 - 0.08018 * dq22 * l6x * s2 * c2 + 0.14454 * dq22 * l6x * c12 * c22 - 0.14454 * dq22 * l6x * c22 - 0.04009 * dq22 * l6y * s12 * s22 - 0.28908 * dq22 * l6y * s12 * s2 * c2 + 0.04009 * dq22 * l6y * s12 * c22 - 0.04009 * dq22 * l6y * s22 * c12 + 0.04009 * dq22 * l6y * s22 - 0.28908 * dq22 * l6y * s2 * c12 * c2 + 0.28908 * dq22 * l6y * s2 * c2 + 0.04009 * dq22 * l6y * c12 * c22 - 0.04009 * dq22 * l6y * c22 + 0.14454 * dq22 * l7x * s12 * s22 - 0.08018 * dq22 * l7x * s12 * s2 * c2 - 0.14454 * dq22 * l7x * s12 * c22 + 0.14454 * dq22 * l7x * s22 * c12 - 0.14454 * dq22 * l7x * s22 - 0.08018 * dq22 * l7x * s2 * c12 * c2 + 0.08018 * dq22 * l7x * s2 * c2 - 0.14454 * dq22 * l7x * c12 * c22 + 0.14454 * dq22 * l7x * c22 + 0.04009 * dq22 * l7y * s12 * s22 + 0.28908 * dq22 * l7y * s12 * s2 * c2 - 0.04009 * dq22 * l7y * s12 * c22 + 0.04009 * dq22 * l7y * s22 * c12 - 0.04009 * dq22 * l7y * s22 + 0.28908 * dq22 * l7y * s2 * c12 * c2 - 0.28908 * dq22 * l7y * s2 * c2 - 0.04009 * dq22 * l7y * c12 * c22 + 0.04009 * dq22 * l7y * c22 + 1.0 * dq22 * l8x * q3 * s12 * s22 - 1.0 * dq22 * l8x * q3 * s12 * c22 + 1.0 * dq22 * l8x * q3 * s22 * c12 - 1.0 * dq22 * l8x * q3 * s22 - 1.0 * dq22 * l8x * q3 * c12 * c22 + 1.0 * dq22 * l8x * q3 * c22 - 0.4318 * dq22 * l8x * s12 * s22 + 0.4318 * dq22 * l8x * s12 * c22 - 0.4318 * dq22 * l8x * s22 * c12 + 0.4318 * dq22 * l8x * s22 + 0.4318 * dq22 * l8x * c12 * c22 - 0.4318 * dq22 * l8x * c22 + 2.0 * dq22 * l8z * q3 * s12 * s2 * c2 + 2.0 * dq22 * l8z * q3 * s2 * c12 * c2 - 2.0 * dq22 * l8z * q3 * s2 * c2 - 0.8636 * dq22 * l8z * s12 * s2 * c2 - 0.8636 * dq22 * l8z * s2 * c12 * c2 + 0.8636 * dq22 * l8z * s2 * c2 + 1.0 * dq22 * l9x * q3 * s12 * s22 - 1.0 * dq22 * l9x * q3 * s12 * c22 + 1.0 * dq22 * l9x * q3 * s22 * c12 - 1.0 * dq22 * l9x * q3 * s22 - 1.0 * dq22 * l9x * q3 * c12 * c22 + 1.0 * dq22 * l9x * q3 * c22 - 0.08018 * dq22 * l9x * s12 * s2 * c2 - 0.08018 * dq22 * l9x * s2 * c12 * c2 + 0.08018 * dq22 * l9x * s2 * c2 + 2.0 * dq22 * l9z * q3 * s12 * s2 * c2 + 2.0 * dq22 * l9z * q3 * s2 * c12 * c2 - 2.0 * dq22 * l9z * q3 * s2 * c2 + 0.04009 * dq22 * l9z * s12 * s22 - 0.04009 * dq22 * l9z * s12 * c22 + 0.04009 * dq22 * l9z * s22 * c12 - 0.04009 * dq22 * l9z * s22 - 0.04009 * dq22 * l9z * c12 * c22 + 0.04009 * dq22 * l9z * c22 + 0.0057946086 * dq22 * m4 * s12 * s22 + 0.0192846035 * dq22 * m4 * s12 * s2 * c2 - 0.0057946086 * dq22 * m4 * s12 * c22 + 0.0057946086 * dq22 * m4 * s22 * c12 - 0.0057946086 * dq22 * m4 * s22 + 0.0192846035 * dq22 * m4 * s2 * c12 * c2 - 0.0192846035 * dq22 * m4 * s2 * c2 - 0.0057946086 * dq22 * m4 * c12 * c22 + 0.0057946086 * dq22 * m4 * c22 + 0.0073212358 * dq22 * m5 * s12 * s22 + 0.0317428563 * dq22 * m5 * s12 * s2 * c2 - 0.0073212358 * dq22 * m5 * s12 * c22 + 0.0073212358 * dq22 * m5 * s22 * c12 - 0.0073212358 * dq22 * m5 * s22 + 0.0317428563 * dq22 * m5 * s2 * c12 * c2 - 0.0317428563 * dq22 * m5 * s2 * c2 - 0.0073212358 * dq22 * m5 * c12 * c22 + 0.0073212358 * dq22 * m5 * c22 + 0.0057946086 * dq22 * m6 * s12 * s22 + 0.0192846035 * dq22 * m6 * s12 * s2 * c2 - 0.0057946086 * dq22 * m6 * s12 * c22 + 0.0057946086 * dq22 * m6 * s22 * c12 - 0.0057946086 * dq22 * m6 * s22 + 0.0192846035 * dq22 * m6 * s2 * c12 * c2 - 0.0192846035 * dq22 * m6 * s2 * c2 - 0.0057946086 * dq22 * m6 * c12 * c22 + 0.0057946086 * dq22 * m6 * c22 + 0.0057946086 * dq22 * m7 * s12 * s22 + 0.0192846035 * dq22 * m7 * s12 * s2 * c2 - 0.0057946086 * dq22 * m7 * s12 * c22 + 0.0057946086 * dq22 * m7 * s22 * c12 - 0.0057946086 * dq22 * m7 * s22 + 0.0192846035 * dq22 * m7 * s2 * c12 * c2 - 0.0192846035 * dq22 * m7 * s2 * c2 - 0.0057946086 * dq22 * m7 * c12 * c22 + 0.0057946086 * dq22 * m7 * c22 + 1.0 * dq22 * m8 * q32 * s12 * s2 * c2 + 1.0 * dq22 * m8 * q32 * s2 * c12 * c2 - 1.0 * dq22 * m8 * q32 * s2 * c2 - 0.8636 * dq22 * m8 * q3 * s12 * s2 * c2 - 0.8636 * dq22 * m8 * q3 * s2 * c12 * c2 + 0.8636 * dq22 * m8 * q3 * s2 * c2 + 0.18645124 * dq22 * m8 * s12 * s2 * c2 + 0.18645124 * dq22 * m8 * s2 * c12 * c2 - 0.18645124 * dq22 * m8 * s2 * c2 + 1.0 * dq22 * m9 * q32 * s12 * s2 * c2 + 1.0 * dq22 * m9 * q32 * s2 * c12 * c2 - 1.0 * dq22 * m9 * q32 * s2 * c2 + 0.04009 * dq22 * m9 * q3 * s12 * s22 - 0.04009 * dq22 * m9 * q3 * s12 * c22 + 0.04009 * dq22 * m9 * q3 * s22 * c12 - 0.04009 * dq22 * m9 * q3 * s22 - 0.04009 * dq22 * m9 * q3 * c12 * c22 + 0.04009 * dq22 * m9 * q3 * c22 - 0.0016072081 * dq22 * m9 * s12 * s2 * c2 - 0.0016072081 * dq22 * m9 * s2 * c12 * c2 + 0.0016072081 * dq22 * m9 * s2 * c2 - 2.0 * dq2 * dq3 * l8x * s12 * s2 * c2 - 2.0 * dq2 * dq3 * l8x * s2 * c12 * c2 + 2.0 * dq2 * dq3 * l8x * s2 * c2 + 2.0 * dq2 * dq3 * l8z * s12 * s22 + 2.0 * dq2 * dq3 * l8z * s22 * c12 + 2.0 * dq2 * dq3 * l8z * c22 - 2.0 * dq2 * dq3 * l9x * s12 * s2 * c2 - 2.0 * dq2 * dq3 * l9x * s2 * c12 * c2 + 2.0 * dq2 * dq3 * l9x * s2 * c2 + 2.0 * dq2 * dq3 * l9z * s12 * s22 + 2.0 * dq2 * dq3 * l9z * s22 * c12 + 2.0 * dq2 * dq3 * l9z * c22 + 2.0 * dq2 * dq3 * m8 * q3 * s12 * s22 + 2.0 * dq2 * dq3 * m8 * q3 * s22 * c12 + 2.0 * dq2 * dq3 * m8 * q3 * c22 - 0.8636 * dq2 * dq3 * m8 * s12 * s22 - 0.8636 * dq2 * dq3 * m8 * s22 * c12 - 0.8636 * dq2 * dq3 * m8 * c22 + 2.0 * dq2 * dq3 * m9 * q3 * s12 * s22 + 2.0 * dq2 * dq3 * m9 * q3 * s22 * c12 + 2.0 * dq2 * dq3 * m9 * q3 * c22 - 0.08018 * dq2 * dq3 * m9 * s12 * s2 * c2 - 0.08018 * dq2 * dq3 * m9 * s2 * c12 * c2 + 0.08018 * dq2 * dq3 * m9 * s2 * c2;
	C(2) = Fc8 * sign(dq3) + Fv8 * dq3 - 1.0 * dq12 * l8x * s12 * s2 * c2 - 1.0 * dq12 * l8x * s2 * c12 * c2 - 1.0 * dq12 * l8z * s12 * c22 - 1.0 * dq12 * l8z * c12 * c22 - 1.0 * dq12 * l9x * s12 * s2 * c2 - 1.0 * dq12 * l9x * s2 * c12 * c2 - 1.0 * dq12 * l9z * s12 * c22 - 1.0 * dq12 * l9z * c12 * c22 - 1.0 * dq12 * m8 * q3 * s12 * c22 - 1.0 * dq12 * m8 * q3 * c12 * c22 + 0.4318 * dq12 * m8 * s12 * c22 + 0.4318 * dq12 * m8 * c12 * c22 - 1.0 * dq12 * m9 * q3 * s12 * c22 - 1.0 * dq12 * m9 * q3 * c12 * c22 - 0.04009 * dq12 * m9 * s12 * s2 * c2 - 0.04009 * dq12 * m9 * s2 * c12 * c2 - 1.0 * dq22 * l8x * s12 * s2 * c2 - 1.0 * dq22 * l8x * s2 * c12 * c2 + 1.0 * dq22 * l8x * s2 * c2 - 1.0 * dq22 * l8z * s12 * c22 - 1.0 * dq22 * l8z * s22 - 1.0 * dq22 * l8z * c12 * c22 - 1.0 * dq22 * l9x * s12 * s2 * c2 - 1.0 * dq22 * l9x * s2 * c12 * c2 + 1.0 * dq22 * l9x * s2 * c2 - 1.0 * dq22 * l9z * s12 * c22 - 1.0 * dq22 * l9z * s22 - 1.0 * dq22 * l9z * c12 * c22 - 1.0 * dq22 * m8 * q3 * s12 * c22 - 1.0 * dq22 * m8 * q3 * s22 - 1.0 * dq22 * m8 * q3 * c12 * c22 + 0.4318 * dq22 * m8 * s12 * c22 + 0.4318 * dq22 * m8 * s22 + 0.4318 * dq22 * m8 * c12 * c22 - 1.0 * dq22 * m9 * q3 * s12 * c22 - 1.0 * dq22 * m9 * q3 * s22 - 1.0 * dq22 * m9 * q3 * c12 * c22 - 0.04009 * dq22 * m9 * s12 * s2 * c2 - 0.04009 * dq22 * m9 * s2 * c12 * c2 + 0.04009 * dq22 * m9 * s2 * c2 - 2.0 * dq2 * dq3 * m8 * s12 * s2 * c2 - 2.0 * dq2 * dq3 * m8 * s2 * c12 * c2 + 2.0 * dq2 * dq3 * m8 * s2 * c2 - 2.0 * dq2 * dq3 * m9 * s12 * s2 * c2 - 2.0 * dq2 * dq3 * m9 * s2 * c12 * c2 + 2.0 * dq2 * dq3 * m9 * s2 * c2;
	C(3) = Fc10 * sign(dq4) + Fv10 * dq4;
	C(4) = 1.0186 * Fc11 * sign(dq5) - 0.8306 * Fc12 * sign(-0.8306 * dq5 + 1.2178 * dq6) - 0.8306 * Fc13 * sign(-0.8306 * dq5 + 1.2178 * dq7) + 1.03754596 * Fv11 * dq5 - 0.8306 * Fv12 * (-0.8306 * dq5 + 1.2178 * dq6) - 0.8306 * Fv13 * (-0.8306 * dq5 + 1.2178 * dq7);
	C(5) = 1.2178 * Fc12 * sign(-0.8306 * dq5 + 1.2178 * dq6) + Fc14 * sign(dq6) - 1.2178 * Fc16 * sign(-1.2178 * dq6 + 1.2178 * dq7) + 1.2178 * Fv12 * (-0.8306 * dq5 + 1.2178 * dq6) + Fv14 * dq6 - 1.2178 * Fv16 * (-1.2178 * dq6 + 1.2178 * dq7);
	C(6) = 1.2178 * Fc13 * sign(-0.8306 * dq5 + 1.2178 * dq7) + Fc15 * sign(dq7) + 1.2178 * Fc16 * sign(-1.2178 * dq6 + 1.2178 * dq7) + 1.2178 * Fv13 * (-0.8306 * dq5 + 1.2178 * dq7) + Fv15 * dq7 + 1.2178 * Fv16 * (-1.2178 * dq6 + 1.2178 * dq7);
	
	return C;

}


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
psmActiveJointsVectorf computePredictedTorque(const psmActiveJointsVectorf& q, const psmActiveJointsVectorf& dq, const psmActiveJointsVectorf& ddq, const bool& withCounterWeight, bool withFriction){
	psmActiveJointsVectorf tau;
	psmDynParamsVectorf delta;
	tau.setZero();

	double m9 = 0.0, l9x = 0.0, l9y = 0.0, l9z = 0.0, L9xx = 0.0, L9xy = 0.0, L9xz = 0.0, L9yy = 0.0, L9yz = 0.0, L9zz = 0.0;
	double Fo1 = 0.0, Fo2 = 0.0, Fo8 = 0.0, Fo10 = 0.0, Fo11 = 0.0, Fo12 = 0.0, Fo13 = 0.0, Fo14 = 0.0, Fo15 = 0.0, Fo16 = 0.0;
	double Fc1 = 0.0, Fc2 = 0.0, Fc8 = 0.0, Fc10 = 0.0, Fc11 = 0.0, Fc12 = 0.0, Fc13 = 0.0, Fc14 = 0.0, Fc15 = 0.0, Fc16 = 0.0;
	double Fv1 = 0.0, Fv2 = 0.0, Fv8 = 0.0, Fv10 = 0.0, Fv11 = 0.0, Fv12 = 0.0, Fv13 = 0.0, Fv14 = 0.0, Fv15 = 0.0, Fv16 = 0.0;
	double Ia10 = 0.0, Ia11 = 0.0, Ia14 = 0.0, Ia15 = 0.0;
	double K10 = 0.0;

	if (withCounterWeight) {
		L9xx = 0.0036083; 
		L9xy = 0.000204;
		L9xz = 0.0009268;
		L9yy = 0.0037699; 
		L9yz = -0.0006548; 
		L9zz = 0.0004368;
		m9 = 0.3499395;
		l9x = 0.0091043;
		l9y = -0.0064588; 
		l9z = -0.0348129;
	}

	if (withFriction) {

		Fo1 = -0.3000027;
		Fo2 = 0.1419199;
		Fo8 = 0.0038872;
		Fo10 = -0.0048904;
		Fo11 = 0.1329238;
		Fo12 = 0.0777171;
		Fo13 = 0.0561429;
		Fo14 = -0.1661489;
		Fo15 = -0.0169157;
		Fo16 = -0.0458516;

		Fc1 = 0.05003;
		Fc2 = 0.0668202;
		Fc8 = 0.4475953;
		Fc10 = 0.002465;
		Fc11 = 0.0069935;
		Fc12 = 0.0007803;
		Fc13 = 0.0021698;
		Fc14 = 0.0062966;
		Fc15 = 0.0061915;
		Fc16 = 0.0001049;

		Fv1 = 0.0905326;
		Fv2 = 0.2180617;
		Fv8 = 1.2035645;
		Fv10 = 0.0011114;
		Fv11 = 0.012171;
		Fv12 = 0.0020153;
		Fv13 = 0.0022465;
		Fv14 = 0.006165;
		Fv15 = 0.0073884;
		Fv16 = 0.0034451;

		Ia10 = -3.1e-6;
		Ia11 = -0.0006152;
		Ia14 = -0.0002932;
		Ia15 = -0.0003634;

		K10 = 0.0001873;

	}

	delta <<
		L1xx, L1xy, L1xz, L1yy, L1yz, L1zz, l1x, l1y, l1z, m1, Fc1, Fv1, Fo1,
		L2xx, L2xy, L2xz, L2yy, L2yz, L2zz, l2x, l2y, l2z, m2, Fc2, Fv2, Fo2,
		L4xx, L4xy, L4xz, L4yy, L4yz, L4zz, l4x, l4y, l4z, m4,
		L5xx, L5xy, L5xz, L5yy, L5yz, L5zz, l5x, l5y, l5z, m5,
		L6xx, L6xy, L6xz, L6yy, L6yz, L6zz, l6x, l6y, l6z, m6,
		L7xx, L7xy, L7xz, L7yy, L7yz, L7zz, l7x, l7y, l7z, m7,
		L8xx, L8xy, L8xz, L8yy, L8yz, L8zz, l8x, l8y, l8z, m8, Fc8, Fv8, Fo8,
		L9xx, L9xy, L9xz, L9yy, L9yz, L9zz, l9x, l9y, l9z, m9, 
		Fc10, Fv10, Fo10, Ia10, K10,
		Fc11, Fv11, Fo11, Ia11, 
		Fc12, Fv12, Fo12, 
		Fc13, Fv13, Fo13, 
		Fc14, Fv14, Fo14, Ia14,
		Fc15, Fv15, Fo15, Ia15, 
		Fc16, Fv16, Fo16;

	double q1 = q(0);
	double q2 = q(1);
	double q3 = q(2);
	double q4 = q(3);
	double q5 = q(4);
	double q6 = q(5);
	double q7 = q(6);

	double s1 = sin(q1);
	double c1 = cos(q1);
	double s2 = sin(q2);
	double c2 = cos(q2);

	double dq1 = dq(0);
	double dq2 = dq(1);
	double dq3 = dq(2);
	double dq4 = dq(3);
	double dq5 = dq(4);
	double dq6 = dq(5);
	double dq7 = dq(6);

	double ddq1 = ddq(0);
	double ddq2 = ddq(1);
	double ddq3 = ddq(2);
	double ddq4 = ddq(3);
	double ddq5 = ddq(4);
	double ddq6 = ddq(5);
	double ddq7 = ddq(6);

	psmRegressorMatrixf H;

	// 79:88 (0-index) : counterweight contribution indexes
	H << 0, 0, 0, 0, 0, ddq1, -9.81 * s1, -9.81 * c1, 0, 0, sign(dq1), dq1, 1, ddq1* pow(c2, 2) - 2 * dq1 * dq2 * s2 * c2, -2 * ddq1 * s2 * c2 + 2 * dq1 * dq2 * pow(s2, 2) - 2 * dq1 * dq2 * pow(c2, 2), ddq2* c2 - pow(dq2, 2) * s2, ddq1* pow(s2, 2) + 2 * dq1 * dq2 * s2 * c2, -ddq2 * s2 - pow(dq2, 2) * c2, 0, -9.81 * s1 * s2, -9.81 * s1 * c2, -9.81 * c1, 0, 0, 0, 0, 1.0 * ddq1, 0, 0, 0, 0, 0, 0, -0.08018 * ddq1 * pow(s1, 2) * s2 - 0.28908 * ddq1 * pow(s1, 2) * c2 - 0.08018 * ddq1 * s2 * pow(c1, 2) - 0.28908 * ddq1 * pow(c1, 2) * c2 + 0.28908 * dq1 * dq2 * pow(s1, 2) * s2 - 0.08018 * dq1 * dq2 * pow(s1, 2) * c2 + 0.28908 * dq1 * dq2 * s2 * pow(c1, 2) - 0.08018 * dq1 * dq2 * pow(c1, 2) * c2 + 9.81 * s1 * pow(s2, 2) + 9.81 * s1 * pow(c2, 2), 0.14454 * ddq2 * pow(s1, 2) * s2 - 0.04009 * ddq2 * pow(s1, 2) * c2 + 0.14454 * ddq2 * s2 * pow(c1, 2) - 0.04009 * ddq2 * pow(c1, 2) * c2 + 0.04009 * pow(dq2, 2) * pow(s1, 2) * s2 + 0.14454 * pow(dq2, 2) * pow(s1, 2) * c2 + 0.04009 * pow(dq2, 2) * s2 * pow(c1, 2) + 0.14454 * pow(dq2, 2) * pow(c1, 2) * c2 - 9.81 * c1, 0.0016072081 * ddq1 * pow(s1, 2) * pow(s2, 2) + 0.0115892172 * ddq1 * pow(s1, 2) * s2 * c2 + 0.0208918116 * ddq1 * pow(s1, 2) * pow(c2, 2) + 0.0016072081 * ddq1 * pow(s2, 2) * pow(c1, 2) + 0.0115892172 * ddq1 * s2 * pow(c1, 2) * c2 + 0.0208918116 * ddq1 * pow(c1, 2) * pow(c2, 2) - 0.0115892172 * dq1 * dq2 * pow(s1, 2) * pow(s2, 2) - 0.038569207 * dq1 * dq2 * pow(s1, 2) * s2 * c2 + 0.0115892172 * dq1 * dq2 * pow(s1, 2) * pow(c2, 2) - 0.0115892172 * dq1 * dq2 * pow(s2, 2) * pow(c1, 2) - 0.038569207 * dq1 * dq2 * s2 * pow(c1, 2) * c2 + 0.0115892172 * dq1 * dq2 * pow(c1, 2) * pow(c2, 2) - 0.3932829 * s1 * s2 - 1.4179374 * s1 * c2, 1.0 * ddq1, 0, 0, 0, 0, 0, 0, -0.08018 * ddq1 * pow(s1, 2) * s2 - 0.36524 * ddq1 * pow(s1, 2) * c2 - 0.08018 * ddq1 * s2 * pow(c1, 2) - 0.36524 * ddq1 * pow(c1, 2) * c2 + 0.36524 * dq1 * dq2 * pow(s1, 2) * s2 - 0.08018 * dq1 * dq2 * pow(s1, 2) * c2 + 0.36524 * dq1 * dq2 * s2 * pow(c1, 2) - 0.08018 * dq1 * dq2 * pow(c1, 2) * c2 + 9.81 * s1 * pow(s2, 2) + 9.81 * s1 * pow(c2, 2), 0.18262 * ddq2 * pow(s1, 2) * s2 - 0.04009 * ddq2 * pow(s1, 2) * c2 + 0.18262 * ddq2 * s2 * pow(c1, 2) - 0.04009 * ddq2 * pow(c1, 2) * c2 + 0.04009 * pow(dq2, 2) * pow(s1, 2) * s2 + 0.18262 * pow(dq2, 2) * pow(s1, 2) * c2 + 0.04009 * pow(dq2, 2) * s2 * pow(c1, 2) + 0.18262 * pow(dq2, 2) * pow(c1, 2) * c2 - 9.81 * c1, 0.0016072081 * ddq1 * pow(s1, 2) * pow(s2, 2) + 0.0146424716 * ddq1 * pow(s1, 2) * s2 * c2 + 0.0333500644 * ddq1 * pow(s1, 2) * pow(c2, 2) + 0.0016072081 * ddq1 * pow(s2, 2) * pow(c1, 2) + 0.0146424716 * ddq1 * s2 * pow(c1, 2) * c2 + 0.0333500644 * ddq1 * pow(c1, 2) * pow(c2, 2) - 0.0146424716 * dq1 * dq2 * pow(s1, 2) * pow(s2, 2) - 0.0634857126 * dq1 * dq2 * pow(s1, 2) * s2 * c2 + 0.0146424716 * dq1 * dq2 * pow(s1, 2) * pow(c2, 2) - 0.0146424716 * dq1 * dq2 * pow(s2, 2) * pow(c1, 2) - 0.0634857126 * dq1 * dq2 * s2 * pow(c1, 2) * c2 + 0.0146424716 * dq1 * dq2 * pow(c1, 2) * pow(c2, 2) - 0.3932829 * s1 * s2 - 1.7915022 * s1 * c2, 1.0 * ddq1 * pow(c2, 2) - 2.0 * dq1 * dq2 * s2 * c2, -2.0 * ddq1 * s2 * c2 + 2.0 * dq1 * dq2 * pow(s2, 2) - 2.0 * dq1 * dq2 * pow(c2, 2), -1.0 * ddq2 * c2 + 1.0 * pow(dq2, 2) * s2, 1.0 * ddq1 * pow(s2, 2) + 2.0 * dq1 * dq2 * s2 * c2, 1.0 * ddq2 * s2 + 1.0 * pow(dq2, 2) * c2, 0, -0.08018 * ddq1 * pow(s1, 2) * pow(s2, 2) - 0.28908 * ddq1 * pow(s1, 2) * s2 * c2 - 0.08018 * ddq1 * pow(s2, 2) * pow(c1, 2) - 0.28908 * ddq1 * s2 * pow(c1, 2) * c2 + 0.28908 * dq1 * dq2 * pow(s1, 2) * pow(s2, 2) - 0.16036 * dq1 * dq2 * pow(s1, 2) * s2 * c2 - 0.28908 * dq1 * dq2 * pow(s1, 2) * pow(c2, 2) + 0.28908 * dq1 * dq2 * pow(s2, 2) * pow(c1, 2) - 0.16036 * dq1 * dq2 * s2 * pow(c1, 2) * c2 - 0.28908 * dq1 * dq2 * pow(c1, 2) * pow(c2, 2) + 9.81 * s1 * pow(s2, 3) + 9.81 * s1 * s2 * pow(c2, 2), -0.08018 * ddq1 * pow(s1, 2) * s2 * c2 - 0.28908 * ddq1 * pow(s1, 2) * pow(c2, 2) - 0.08018 * ddq1 * s2 * pow(c1, 2) * c2 - 0.28908 * ddq1 * pow(c1, 2) * pow(c2, 2) + 0.08018 * dq1 * dq2 * pow(s1, 2) * pow(s2, 2) + 0.57816 * dq1 * dq2 * pow(s1, 2) * s2 * c2 - 0.08018 * dq1 * dq2 * pow(s1, 2) * pow(c2, 2) + 0.08018 * dq1 * dq2 * pow(s2, 2) * pow(c1, 2) + 0.57816 * dq1 * dq2 * s2 * pow(c1, 2) * c2 - 0.08018 * dq1 * dq2 * pow(c1, 2) * pow(c2, 2) + 9.81 * s1 * pow(s2, 2) * c2 + 9.81 * s1 * pow(c2, 3), 0.14454 * ddq2 * pow(s1, 2) * s2 - 0.04009 * ddq2 * pow(s1, 2) * c2 + 0.14454 * ddq2 * s2 * pow(c1, 2) - 0.04009 * ddq2 * pow(c1, 2) * c2 + 0.04009 * pow(dq2, 2) * pow(s1, 2) * s2 + 0.14454 * pow(dq2, 2) * pow(s1, 2) * c2 + 0.04009 * pow(dq2, 2) * s2 * pow(c1, 2) + 0.14454 * pow(dq2, 2) * pow(c1, 2) * c2 - 9.81 * c1, 0.0016072081 * ddq1 * pow(s1, 2) * pow(s2, 2) + 0.0115892172 * ddq1 * pow(s1, 2) * s2 * c2 + 0.0208918116 * ddq1 * pow(s1, 2) * pow(c2, 2) + 0.0016072081 * ddq1 * pow(s2, 2) * pow(c1, 2) + 0.0115892172 * ddq1 * s2 * pow(c1, 2) * c2 + 0.0208918116 * ddq1 * pow(c1, 2) * pow(c2, 2) - 0.0115892172 * dq1 * dq2 * pow(s1, 2) * pow(s2, 2) - 0.038569207 * dq1 * dq2 * pow(s1, 2) * s2 * c2 + 0.0115892172 * dq1 * dq2 * pow(s1, 2) * pow(c2, 2) - 0.0115892172 * dq1 * dq2 * pow(s2, 2) * pow(c1, 2) - 0.038569207 * dq1 * dq2 * s2 * pow(c1, 2) * c2 + 0.0115892172 * dq1 * dq2 * pow(c1, 2) * pow(c2, 2) - 0.3932829 * s1 * s2 - 1.4179374 * s1 * c2, 1.0 * ddq1 * pow(c2, 2) - 2.0 * dq1 * dq2 * s2 * c2, -2.0 * ddq1 * s2 * c2 + 2.0 * dq1 * dq2 * pow(s2, 2) - 2.0 * dq1 * dq2 * pow(c2, 2), 1.0 * ddq2 * c2 - 1.0 * pow(dq2, 2) * s2, 1.0 * ddq1 * pow(s2, 2) + 2.0 * dq1 * dq2 * s2 * c2, -1.0 * ddq2 * s2 - 1.0 * pow(dq2, 2) * c2, 0, 0.08018 * ddq1 * pow(s1, 2) * pow(s2, 2) + 0.28908 * ddq1 * pow(s1, 2) * s2 * c2 + 0.08018 * ddq1 * pow(s2, 2) * pow(c1, 2) + 0.28908 * ddq1 * s2 * pow(c1, 2) * c2 - 0.28908 * dq1 * dq2 * pow(s1, 2) * pow(s2, 2) + 0.16036 * dq1 * dq2 * pow(s1, 2) * s2 * c2 + 0.28908 * dq1 * dq2 * pow(s1, 2) * pow(c2, 2) - 0.28908 * dq1 * dq2 * pow(s2, 2) * pow(c1, 2) + 0.16036 * dq1 * dq2 * s2 * pow(c1, 2) * c2 + 0.28908 * dq1 * dq2 * pow(c1, 2) * pow(c2, 2) - 9.81 * s1 * pow(s2, 3) - 9.81 * s1 * s2 * pow(c2, 2), 0.08018 * ddq1 * pow(s1, 2) * s2 * c2 + 0.28908 * ddq1 * pow(s1, 2) * pow(c2, 2) + 0.08018 * ddq1 * s2 * pow(c1, 2) * c2 + 0.28908 * ddq1 * pow(c1, 2) * pow(c2, 2) - 0.08018 * dq1 * dq2 * pow(s1, 2) * pow(s2, 2) - 0.57816 * dq1 * dq2 * pow(s1, 2) * s2 * c2 + 0.08018 * dq1 * dq2 * pow(s1, 2) * pow(c2, 2) - 0.08018 * dq1 * dq2 * pow(s2, 2) * pow(c1, 2) - 0.57816 * dq1 * dq2 * s2 * pow(c1, 2) * c2 + 0.08018 * dq1 * dq2 * pow(c1, 2) * pow(c2, 2) - 9.81 * s1 * pow(s2, 2) * c2 - 9.81 * s1 * pow(c2, 3), 0.14454 * ddq2 * pow(s1, 2) * s2 - 0.04009 * ddq2 * pow(s1, 2) * c2 + 0.14454 * ddq2 * s2 * pow(c1, 2) - 0.04009 * ddq2 * pow(c1, 2) * c2 + 0.04009 * pow(dq2, 2) * pow(s1, 2) * s2 + 0.14454 * pow(dq2, 2) * pow(s1, 2) * c2 + 0.04009 * pow(dq2, 2) * s2 * pow(c1, 2) + 0.14454 * pow(dq2, 2) * pow(c1, 2) * c2 - 9.81 * c1, 0.0016072081 * ddq1 * pow(s1, 2) * pow(s2, 2) + 0.0115892172 * ddq1 * pow(s1, 2) * s2 * c2 + 0.0208918116 * ddq1 * pow(s1, 2) * pow(c2, 2) + 0.0016072081 * ddq1 * pow(s2, 2) * pow(c1, 2) + 0.0115892172 * ddq1 * s2 * pow(c1, 2) * c2 + 0.0208918116 * ddq1 * pow(c1, 2) * pow(c2, 2) - 0.0115892172 * dq1 * dq2 * pow(s1, 2) * pow(s2, 2) - 0.038569207 * dq1 * dq2 * pow(s1, 2) * s2 * c2 + 0.0115892172 * dq1 * dq2 * pow(s1, 2) * pow(c2, 2) - 0.0115892172 * dq1 * dq2 * pow(s2, 2) * pow(c1, 2) - 0.038569207 * dq1 * dq2 * s2 * pow(c1, 2) * c2 + 0.0115892172 * dq1 * dq2 * pow(c1, 2) * pow(c2, 2) - 0.3932829 * s1 * s2 - 1.4179374 * s1 * c2, 1.0 * ddq1 * pow(c2, 2) - 2.0 * dq1 * dq2 * s2 * c2, 1.0 * ddq2 * c2 - 1.0 * pow(dq2, 2) * s2, -2.0 * ddq1 * s2 * c2 + 2.0 * dq1 * dq2 * pow(s2, 2) - 2.0 * dq1 * dq2 * pow(c2, 2), 0, -1.0 * ddq2 * s2 - 1.0 * pow(dq2, 2) * c2, 1.0 * ddq1 * pow(s2, 2) + 2.0 * dq1 * dq2 * s2 * c2, 2.0 * ddq1 * q3 * pow(s1, 2) * s2 * c2 + 2.0 * ddq1 * q3 * s2 * pow(c1, 2) * c2 - 0.8636 * ddq1 * pow(s1, 2) * s2 * c2 - 0.8636 * ddq1 * s2 * pow(c1, 2) * c2 - 2.0 * dq1 * dq2 * q3 * pow(s1, 2) * pow(s2, 2) + 2.0 * dq1 * dq2 * q3 * pow(s1, 2) * pow(c2, 2) - 2.0 * dq1 * dq2 * q3 * pow(s2, 2) * pow(c1, 2) + 2.0 * dq1 * dq2 * q3 * pow(c1, 2) * pow(c2, 2) + 0.8636 * dq1 * dq2 * pow(s1, 2) * pow(s2, 2) - 0.8636 * dq1 * dq2 * pow(s1, 2) * pow(c2, 2) + 0.8636 * dq1 * dq2 * pow(s2, 2) * pow(c1, 2) - 0.8636 * dq1 * dq2 * pow(c1, 2) * pow(c2, 2) + 2.0 * dq1 * dq3 * pow(s1, 2) * s2 * c2 + 2.0 * dq1 * dq3 * s2 * pow(c1, 2) * c2 + 9.81 * s1 * pow(s2, 3) + 9.81 * s1 * s2 * pow(c2, 2), 1.0 * ddq2 * q3 * pow(s1, 2) * s2 + 1.0 * ddq2 * q3 * s2 * pow(c1, 2) - 0.4318 * ddq2 * pow(s1, 2) * s2 - 0.4318 * ddq2 * s2 * pow(c1, 2) - 1.0 * ddq3 * pow(s1, 2) * c2 - 1.0 * ddq3 * pow(c1, 2) * c2 + 1.0 * pow(dq2, 2) * q3 * pow(s1, 2) * c2 + 1.0 * pow(dq2, 2) * q3 * pow(c1, 2) * c2 - 0.4318 * pow(dq2, 2) * pow(s1, 2) * c2 - 0.4318 * pow(dq2, 2) * pow(c1, 2) * c2 + 2.0 * dq2 * dq3 * pow(s1, 2) * s2 + 2.0 * dq2 * dq3 * s2 * pow(c1, 2) + 9.81 * c1, 2.0 * ddq1 * q3 * pow(s1, 2) * pow(c2, 2) + 2.0 * ddq1 * q3 * pow(c1, 2) * pow(c2, 2) - 0.8636 * ddq1 * pow(s1, 2) * pow(c2, 2) - 0.8636 * ddq1 * pow(c1, 2) * pow(c2, 2) - 4.0 * dq1 * dq2 * q3 * pow(s1, 2) * s2 * c2 - 4.0 * dq1 * dq2 * q3 * s2 * pow(c1, 2) * c2 + 1.7272 * dq1 * dq2 * pow(s1, 2) * s2 * c2 + 1.7272 * dq1 * dq2 * s2 * pow(c1, 2) * c2 + 2.0 * dq1 * dq3 * pow(s1, 2) * pow(c2, 2) + 2.0 * dq1 * dq3 * pow(c1, 2) * pow(c2, 2) + 9.81 * s1 * pow(s2, 2) * c2 + 9.81 * s1 * pow(c2, 3), 1.0 * ddq1 * pow(q3, 2) * pow(s1, 2) * pow(c2, 2) + 1.0 * ddq1 * pow(q3, 2) * pow(c1, 2) * pow(c2, 2) - 0.8636 * ddq1 * q3 * pow(s1, 2) * pow(c2, 2) - 0.8636 * ddq1 * q3 * pow(c1, 2) * pow(c2, 2) + 0.18645124 * ddq1 * pow(s1, 2) * pow(c2, 2) + 0.18645124 * ddq1 * pow(c1, 2) * pow(c2, 2) - 2.0 * dq1 * dq2 * pow(q3, 2) * pow(s1, 2) * s2 * c2 - 2.0 * dq1 * dq2 * pow(q3, 2) * s2 * pow(c1, 2) * c2 + 1.7272 * dq1 * dq2 * q3 * pow(s1, 2) * s2 * c2 + 1.7272 * dq1 * dq2 * q3 * s2 * pow(c1, 2) * c2 - 0.37290248 * dq1 * dq2 * pow(s1, 2) * s2 * c2 - 0.37290248 * dq1 * dq2 * s2 * pow(c1, 2) * c2 + 2.0 * dq1 * dq3 * q3 * pow(s1, 2) * pow(c2, 2) + 2.0 * dq1 * dq3 * q3 * pow(c1, 2) * pow(c2, 2) - 0.8636 * dq1 * dq3 * pow(s1, 2) * pow(c2, 2) - 0.8636 * dq1 * dq3 * pow(c1, 2) * pow(c2, 2) + 9.81 * q3 * s1 * pow(s2, 2) * c2 + 9.81 * q3 * s1 * pow(c2, 3) + 0.3932829 * s1 * pow(s2, 3) - 2.8180206 * s1 * pow(s2, 2) * c2 + 0.3932829 * s1 * s2 * pow(c2, 2) - 0.3932829 * s1 * s2 - 2.8180206 * s1 * pow(c2, 3) - 1.4179374 * s1 * c2, 0, 0, 0, 1.0 * ddq1 * pow(c2, 2) - 2.0 * dq1 * dq2 * s2 * c2, -1.0 * ddq2 * c2 + 1.0 * pow(dq2, 2) * s2, -2.0 * ddq1 * s2 * c2 + 2.0 * dq1 * dq2 * pow(s2, 2) - 2.0 * dq1 * dq2 * pow(c2, 2), 0, 1.0 * ddq2 * s2 + 1.0 * pow(dq2, 2) * c2, 1.0 * ddq1 * pow(s2, 2) + 2.0 * dq1 * dq2 * s2 * c2, 2.0 * ddq1 * q3 * pow(s1, 2) * s2 * c2 + 2.0 * ddq1 * q3 * s2 * pow(c1, 2) * c2 + 0.08018 * ddq1 * pow(s1, 2) * pow(s2, 2) + 0.08018 * ddq1 * pow(s2, 2) * pow(c1, 2) - 2.0 * dq1 * dq2 * q3 * pow(s1, 2) * pow(s2, 2) + 2.0 * dq1 * dq2 * q3 * pow(s1, 2) * pow(c2, 2) - 2.0 * dq1 * dq2 * q3 * pow(s2, 2) * pow(c1, 2) + 2.0 * dq1 * dq2 * q3 * pow(c1, 2) * pow(c2, 2) + 0.16036 * dq1 * dq2 * pow(s1, 2) * s2 * c2 + 0.16036 * dq1 * dq2 * s2 * pow(c1, 2) * c2 + 2.0 * dq1 * dq3 * pow(s1, 2) * s2 * c2 + 2.0 * dq1 * dq3 * s2 * pow(c1, 2) * c2 - 9.81 * s1 * s2, -1.0 * ddq2 * q3 * pow(s1, 2) * s2 - 1.0 * ddq2 * q3 * s2 * pow(c1, 2) + 0.04009 * ddq2 * pow(s1, 2) * c2 + 0.04009 * ddq2 * pow(c1, 2) * c2 + 1.0 * ddq3 * pow(s1, 2) * c2 + 1.0 * ddq3 * pow(c1, 2) * c2 - 1.0 * pow(dq2, 2) * q3 * pow(s1, 2) * c2 - 1.0 * pow(dq2, 2) * q3 * pow(c1, 2) * c2 - 0.04009 * pow(dq2, 2) * pow(s1, 2) * s2 - 0.04009 * pow(dq2, 2) * s2 * pow(c1, 2) - 2.0 * dq2 * dq3 * pow(s1, 2) * s2 - 2.0 * dq2 * dq3 * s2 * pow(c1, 2) + 9.81 * c1, 2.0 * ddq1 * q3 * pow(s1, 2) * pow(c2, 2) + 2.0 * ddq1 * q3 * pow(c1, 2) * pow(c2, 2) + 0.08018 * ddq1 * pow(s1, 2) * s2 * c2 + 0.08018 * ddq1 * s2 * pow(c1, 2) * c2 - 4.0 * dq1 * dq2 * q3 * pow(s1, 2) * s2 * c2 - 4.0 * dq1 * dq2 * q3 * s2 * pow(c1, 2) * c2 - 0.08018 * dq1 * dq2 * pow(s1, 2) * pow(s2, 2) + 0.08018 * dq1 * dq2 * pow(s1, 2) * pow(c2, 2) - 0.08018 * dq1 * dq2 * pow(s2, 2) * pow(c1, 2) + 0.08018 * dq1 * dq2 * pow(c1, 2) * pow(c2, 2) + 2.0 * dq1 * dq3 * pow(s1, 2) * pow(c2, 2) + 2.0 * dq1 * dq3 * pow(c1, 2) * pow(c2, 2) - 9.81 * s1 * c2, 1.0 * ddq1 * pow(q3, 2) * pow(s1, 2) * pow(c2, 2) + 1.0 * ddq1 * pow(q3, 2) * pow(c1, 2) * pow(c2, 2) + 0.08018 * ddq1 * q3 * pow(s1, 2) * s2 * c2 + 0.08018 * ddq1 * q3 * s2 * pow(c1, 2) * c2 + 0.0016072081 * ddq1 * pow(s1, 2) * pow(s2, 2) + 0.0016072081 * ddq1 * pow(s2, 2) * pow(c1, 2) - 2.0 * dq1 * dq2 * pow(q3, 2) * pow(s1, 2) * s2 * c2 - 2.0 * dq1 * dq2 * pow(q3, 2) * s2 * pow(c1, 2) * c2 - 0.08018 * dq1 * dq2 * q3 * pow(s1, 2) * pow(s2, 2) + 0.08018 * dq1 * dq2 * q3 * pow(s1, 2) * pow(c2, 2) - 0.08018 * dq1 * dq2 * q3 * pow(s2, 2) * pow(c1, 2) + 0.08018 * dq1 * dq2 * q3 * pow(c1, 2) * pow(c2, 2) + 0.0032144162 * dq1 * dq2 * pow(s1, 2) * s2 * c2 + 0.0032144162 * dq1 * dq2 * s2 * pow(c1, 2) * c2 + 2.0 * dq1 * dq3 * q3 * pow(s1, 2) * pow(c2, 2) + 2.0 * dq1 * dq3 * q3 * pow(c1, 2) * pow(c2, 2) + 0.08018 * dq1 * dq3 * pow(s1, 2) * s2 * c2 + 0.08018 * dq1 * dq3 * s2 * pow(c1, 2) * c2 - 9.81 * q3 * s1 * c2 - 0.3932829 * s1 * s2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pow(dq1, 2)* s2* c2, -pow(dq1, 2) * pow(s2, 2) + pow(dq1, 2) * pow(c2, 2), ddq1* c2, -pow(dq1, 2) * s2 * c2, -ddq1 * s2, ddq2, 9.81 * c1 * c2, -9.81 * s2 * c1, 0, 0, sign(dq2), dq2, 1, 0, 0, 0, 0, 0, 0, 0, -0.14454 * pow(dq1, 2) * pow(s1, 2) * s2 + 0.04009 * pow(dq1, 2) * pow(s1, 2) * c2 - 0.14454 * pow(dq1, 2) * s2 * pow(c1, 2) + 0.04009 * pow(dq1, 2) * pow(c1, 2) * c2, 0.14454 * ddq1 * pow(s1, 2) * s2 - 0.04009 * ddq1 * pow(s1, 2) * c2 + 0.14454 * ddq1 * s2 * pow(c1, 2) - 0.04009 * ddq1 * pow(c1, 2) * c2, 0.0208918116 * ddq2 * pow(s1, 2) * pow(s2, 2) - 0.0115892172 * ddq2 * pow(s1, 2) * s2 * c2 + 0.0016072081 * ddq2 * pow(s1, 2) * pow(c2, 2) + 0.0208918116 * ddq2 * pow(s2, 2) * pow(c1, 2) + 0.0016072081 * ddq2 * pow(s2, 2) - 0.0115892172 * ddq2 * s2 * pow(c1, 2) * c2 + 0.0115892172 * ddq2 * s2 * c2 + 0.0016072081 * ddq2 * pow(c1, 2) * pow(c2, 2) + 0.0208918116 * ddq2 * pow(c2, 2) + 0.0057946086 * pow(dq1, 2) * pow(s1, 2) * pow(s2, 2) + 0.0192846035 * pow(dq1, 2) * pow(s1, 2) * s2 * c2 - 0.0057946086 * pow(dq1, 2) * pow(s1, 2) * pow(c2, 2) + 0.0057946086 * pow(dq1, 2) * pow(s2, 2) * pow(c1, 2) + 0.0192846035 * pow(dq1, 2) * s2 * pow(c1, 2) * c2 - 0.0057946086 * pow(dq1, 2) * pow(c1, 2) * pow(c2, 2) + 0.0057946086 * pow(dq2, 2) * pow(s1, 2) * pow(s2, 2) + 0.0192846035 * pow(dq2, 2) * pow(s1, 2) * s2 * c2 - 0.0057946086 * pow(dq2, 2) * pow(s1, 2) * pow(c2, 2) + 0.0057946086 * pow(dq2, 2) * pow(s2, 2) * pow(c1, 2) - 0.0057946086 * pow(dq2, 2) * pow(s2, 2) + 0.0192846035 * pow(dq2, 2) * s2 * pow(c1, 2) * c2 - 0.0192846035 * pow(dq2, 2) * s2 * c2 - 0.0057946086 * pow(dq2, 2) * pow(c1, 2) * pow(c2, 2) + 0.0057946086 * pow(dq2, 2) * pow(c2, 2) - 1.4179374 * s2 * c1 + 0.3932829 * c1 * c2, 0, 0, 0, 0, 0, 0, 0, -0.18262 * pow(dq1, 2) * pow(s1, 2) * s2 + 0.04009 * pow(dq1, 2) * pow(s1, 2) * c2 - 0.18262 * pow(dq1, 2) * s2 * pow(c1, 2) + 0.04009 * pow(dq1, 2) * pow(c1, 2) * c2, 0.18262 * ddq1 * pow(s1, 2) * s2 - 0.04009 * ddq1 * pow(s1, 2) * c2 + 0.18262 * ddq1 * s2 * pow(c1, 2) - 0.04009 * ddq1 * pow(c1, 2) * c2, 0.0333500644 * ddq2 * pow(s1, 2) * pow(s2, 2) - 0.0146424716 * ddq2 * pow(s1, 2) * s2 * c2 + 0.0016072081 * ddq2 * pow(s1, 2) * pow(c2, 2) + 0.0333500644 * ddq2 * pow(s2, 2) * pow(c1, 2) + 0.0016072081 * ddq2 * pow(s2, 2) - 0.0146424716 * ddq2 * s2 * pow(c1, 2) * c2 + 0.0146424716 * ddq2 * s2 * c2 + 0.0016072081 * ddq2 * pow(c1, 2) * pow(c2, 2) + 0.0333500644 * ddq2 * pow(c2, 2) + 0.0073212358 * pow(dq1, 2) * pow(s1, 2) * pow(s2, 2) + 0.0317428563 * pow(dq1, 2) * pow(s1, 2) * s2 * c2 - 0.0073212358 * pow(dq1, 2) * pow(s1, 2) * pow(c2, 2) + 0.0073212358 * pow(dq1, 2) * pow(s2, 2) * pow(c1, 2) + 0.0317428563 * pow(dq1, 2) * s2 * pow(c1, 2) * c2 - 0.0073212358 * pow(dq1, 2) * pow(c1, 2) * pow(c2, 2) + 0.0073212358 * pow(dq2, 2) * pow(s1, 2) * pow(s2, 2) + 0.0317428563 * pow(dq2, 2) * pow(s1, 2) * s2 * c2 - 0.0073212358 * pow(dq2, 2) * pow(s1, 2) * pow(c2, 2) + 0.0073212358 * pow(dq2, 2) * pow(s2, 2) * pow(c1, 2) - 0.0073212358 * pow(dq2, 2) * pow(s2, 2) + 0.0317428563 * pow(dq2, 2) * s2 * pow(c1, 2) * c2 - 0.0317428563 * pow(dq2, 2) * s2 * c2 - 0.0073212358 * pow(dq2, 2) * pow(c1, 2) * pow(c2, 2) + 0.0073212358 * pow(dq2, 2) * pow(c2, 2) - 1.7915022 * s2 * c1 + 0.3932829 * c1 * c2, 1.0 * pow(dq1, 2) * s2 * c2, -1.0 * pow(dq1, 2) * pow(s2, 2) + 1.0 * pow(dq1, 2) * pow(c2, 2), -1.0 * ddq1 * c2, -1.0 * pow(dq1, 2) * s2 * c2, 1.0 * ddq1 * s2, 1.0 * ddq2, 0.28908 * ddq2 * pow(s1, 2) * s2 * c2 - 0.08018 * ddq2 * pow(s1, 2) * pow(c2, 2) - 0.08018 * ddq2 * pow(s2, 2) + 0.28908 * ddq2 * s2 * pow(c1, 2) * c2 - 0.28908 * ddq2 * s2 * c2 - 0.08018 * ddq2 * pow(c1, 2) * pow(c2, 2) - 0.14454 * pow(dq1, 2) * pow(s1, 2) * pow(s2, 2) + 0.08018 * pow(dq1, 2) * pow(s1, 2) * s2 * c2 + 0.14454 * pow(dq1, 2) * pow(s1, 2) * pow(c2, 2) - 0.14454 * pow(dq1, 2) * pow(s2, 2) * pow(c1, 2) + 0.08018 * pow(dq1, 2) * s2 * pow(c1, 2) * c2 + 0.14454 * pow(dq1, 2) * pow(c1, 2) * pow(c2, 2) - 0.14454 * pow(dq2, 2) * pow(s1, 2) * pow(s2, 2) + 0.08018 * pow(dq2, 2) * pow(s1, 2) * s2 * c2 + 0.14454 * pow(dq2, 2) * pow(s1, 2) * pow(c2, 2) - 0.14454 * pow(dq2, 2) * pow(s2, 2) * pow(c1, 2) + 0.14454 * pow(dq2, 2) * pow(s2, 2) + 0.08018 * pow(dq2, 2) * s2 * pow(c1, 2) * c2 - 0.08018 * pow(dq2, 2) * s2 * c2 + 0.14454 * pow(dq2, 2) * pow(c1, 2) * pow(c2, 2) - 0.14454 * pow(dq2, 2) * pow(c2, 2) - 9.81 * pow(s2, 2) * c1 * c2 - 9.81 * c1 * pow(c2, 3), -0.28908 * ddq2 * pow(s1, 2) * pow(s2, 2) + 0.08018 * ddq2 * pow(s1, 2) * s2 * c2 - 0.28908 * ddq2 * pow(s2, 2) * pow(c1, 2) + 0.08018 * ddq2 * s2 * pow(c1, 2) * c2 - 0.08018 * ddq2 * s2 * c2 - 0.28908 * ddq2 * pow(c2, 2) - 0.04009 * pow(dq1, 2) * pow(s1, 2) * pow(s2, 2) - 0.28908 * pow(dq1, 2) * pow(s1, 2) * s2 * c2 + 0.04009 * pow(dq1, 2) * pow(s1, 2) * pow(c2, 2) - 0.04009 * pow(dq1, 2) * pow(s2, 2) * pow(c1, 2) - 0.28908 * pow(dq1, 2) * s2 * pow(c1, 2) * c2 + 0.04009 * pow(dq1, 2) * pow(c1, 2) * pow(c2, 2) - 0.04009 * pow(dq2, 2) * pow(s1, 2) * pow(s2, 2) - 0.28908 * pow(dq2, 2) * pow(s1, 2) * s2 * c2 + 0.04009 * pow(dq2, 2) * pow(s1, 2) * pow(c2, 2) - 0.04009 * pow(dq2, 2) * pow(s2, 2) * pow(c1, 2) + 0.04009 * pow(dq2, 2) * pow(s2, 2) - 0.28908 * pow(dq2, 2) * s2 * pow(c1, 2) * c2 + 0.28908 * pow(dq2, 2) * s2 * c2 + 0.04009 * pow(dq2, 2) * pow(c1, 2) * pow(c2, 2) - 0.04009 * pow(dq2, 2) * pow(c2, 2) + 9.81 * pow(s2, 3) * c1 + 9.81 * s2 * c1 * pow(c2, 2), 0.14454 * ddq1 * pow(s1, 2) * s2 - 0.04009 * ddq1 * pow(s1, 2) * c2 + 0.14454 * ddq1 * s2 * pow(c1, 2) - 0.04009 * ddq1 * pow(c1, 2) * c2, 0.0208918116 * ddq2 * pow(s1, 2) * pow(s2, 2) - 0.0115892172 * ddq2 * pow(s1, 2) * s2 * c2 + 0.0016072081 * ddq2 * pow(s1, 2) * pow(c2, 2) + 0.0208918116 * ddq2 * pow(s2, 2) * pow(c1, 2) + 0.0016072081 * ddq2 * pow(s2, 2) - 0.0115892172 * ddq2 * s2 * pow(c1, 2) * c2 + 0.0115892172 * ddq2 * s2 * c2 + 0.0016072081 * ddq2 * pow(c1, 2) * pow(c2, 2) + 0.0208918116 * ddq2 * pow(c2, 2) + 0.0057946086 * pow(dq1, 2) * pow(s1, 2) * pow(s2, 2) + 0.0192846035 * pow(dq1, 2) * pow(s1, 2) * s2 * c2 - 0.0057946086 * pow(dq1, 2) * pow(s1, 2) * pow(c2, 2) + 0.0057946086 * pow(dq1, 2) * pow(s2, 2) * pow(c1, 2) + 0.0192846035 * pow(dq1, 2) * s2 * pow(c1, 2) * c2 - 0.0057946086 * pow(dq1, 2) * pow(c1, 2) * pow(c2, 2) + 0.0057946086 * pow(dq2, 2) * pow(s1, 2) * pow(s2, 2) + 0.0192846035 * pow(dq2, 2) * pow(s1, 2) * s2 * c2 - 0.0057946086 * pow(dq2, 2) * pow(s1, 2) * pow(c2, 2) + 0.0057946086 * pow(dq2, 2) * pow(s2, 2) * pow(c1, 2) - 0.0057946086 * pow(dq2, 2) * pow(s2, 2) + 0.0192846035 * pow(dq2, 2) * s2 * pow(c1, 2) * c2 - 0.0192846035 * pow(dq2, 2) * s2 * c2 - 0.0057946086 * pow(dq2, 2) * pow(c1, 2) * pow(c2, 2) + 0.0057946086 * pow(dq2, 2) * pow(c2, 2) - 1.4179374 * s2 * c1 + 0.3932829 * c1 * c2, 1.0 * pow(dq1, 2) * s2 * c2, -1.0 * pow(dq1, 2) * pow(s2, 2) + 1.0 * pow(dq1, 2) * pow(c2, 2), 1.0 * ddq1 * c2, -1.0 * pow(dq1, 2) * s2 * c2, -1.0 * ddq1 * s2, 1.0 * ddq2, -0.28908 * ddq2 * pow(s1, 2) * s2 * c2 + 0.08018 * ddq2 * pow(s1, 2) * pow(c2, 2) + 0.08018 * ddq2 * pow(s2, 2) - 0.28908 * ddq2 * s2 * pow(c1, 2) * c2 + 0.28908 * ddq2 * s2 * c2 + 0.08018 * ddq2 * pow(c1, 2) * pow(c2, 2) + 0.14454 * pow(dq1, 2) * pow(s1, 2) * pow(s2, 2) - 0.08018 * pow(dq1, 2) * pow(s1, 2) * s2 * c2 - 0.14454 * pow(dq1, 2) * pow(s1, 2) * pow(c2, 2) + 0.14454 * pow(dq1, 2) * pow(s2, 2) * pow(c1, 2) - 0.08018 * pow(dq1, 2) * s2 * pow(c1, 2) * c2 - 0.14454 * pow(dq1, 2) * pow(c1, 2) * pow(c2, 2) + 0.14454 * pow(dq2, 2) * pow(s1, 2) * pow(s2, 2) - 0.08018 * pow(dq2, 2) * pow(s1, 2) * s2 * c2 - 0.14454 * pow(dq2, 2) * pow(s1, 2) * pow(c2, 2) + 0.14454 * pow(dq2, 2) * pow(s2, 2) * pow(c1, 2) - 0.14454 * pow(dq2, 2) * pow(s2, 2) - 0.08018 * pow(dq2, 2) * s2 * pow(c1, 2) * c2 + 0.08018 * pow(dq2, 2) * s2 * c2 - 0.14454 * pow(dq2, 2) * pow(c1, 2) * pow(c2, 2) + 0.14454 * pow(dq2, 2) * pow(c2, 2) + 9.81 * pow(s2, 2) * c1 * c2 + 9.81 * c1 * pow(c2, 3), 0.28908 * ddq2 * pow(s1, 2) * pow(s2, 2) - 0.08018 * ddq2 * pow(s1, 2) * s2 * c2 + 0.28908 * ddq2 * pow(s2, 2) * pow(c1, 2) - 0.08018 * ddq2 * s2 * pow(c1, 2) * c2 + 0.08018 * ddq2 * s2 * c2 + 0.28908 * ddq2 * pow(c2, 2) + 0.04009 * pow(dq1, 2) * pow(s1, 2) * pow(s2, 2) + 0.28908 * pow(dq1, 2) * pow(s1, 2) * s2 * c2 - 0.04009 * pow(dq1, 2) * pow(s1, 2) * pow(c2, 2) + 0.04009 * pow(dq1, 2) * pow(s2, 2) * pow(c1, 2) + 0.28908 * pow(dq1, 2) * s2 * pow(c1, 2) * c2 - 0.04009 * pow(dq1, 2) * pow(c1, 2) * pow(c2, 2) + 0.04009 * pow(dq2, 2) * pow(s1, 2) * pow(s2, 2) + 0.28908 * pow(dq2, 2) * pow(s1, 2) * s2 * c2 - 0.04009 * pow(dq2, 2) * pow(s1, 2) * pow(c2, 2) + 0.04009 * pow(dq2, 2) * pow(s2, 2) * pow(c1, 2) - 0.04009 * pow(dq2, 2) * pow(s2, 2) + 0.28908 * pow(dq2, 2) * s2 * pow(c1, 2) * c2 - 0.28908 * pow(dq2, 2) * s2 * c2 - 0.04009 * pow(dq2, 2) * pow(c1, 2) * pow(c2, 2) + 0.04009 * pow(dq2, 2) * pow(c2, 2) - 9.81 * pow(s2, 3) * c1 - 9.81 * s2 * c1 * pow(c2, 2), 0.14454 * ddq1 * pow(s1, 2) * s2 - 0.04009 * ddq1 * pow(s1, 2) * c2 + 0.14454 * ddq1 * s2 * pow(c1, 2) - 0.04009 * ddq1 * pow(c1, 2) * c2, 0.0208918116 * ddq2 * pow(s1, 2) * pow(s2, 2) - 0.0115892172 * ddq2 * pow(s1, 2) * s2 * c2 + 0.0016072081 * ddq2 * pow(s1, 2) * pow(c2, 2) + 0.0208918116 * ddq2 * pow(s2, 2) * pow(c1, 2) + 0.0016072081 * ddq2 * pow(s2, 2) - 0.0115892172 * ddq2 * s2 * pow(c1, 2) * c2 + 0.0115892172 * ddq2 * s2 * c2 + 0.0016072081 * ddq2 * pow(c1, 2) * pow(c2, 2) + 0.0208918116 * ddq2 * pow(c2, 2) + 0.0057946086 * pow(dq1, 2) * pow(s1, 2) * pow(s2, 2) + 0.0192846035 * pow(dq1, 2) * pow(s1, 2) * s2 * c2 - 0.0057946086 * pow(dq1, 2) * pow(s1, 2) * pow(c2, 2) + 0.0057946086 * pow(dq1, 2) * pow(s2, 2) * pow(c1, 2) + 0.0192846035 * pow(dq1, 2) * s2 * pow(c1, 2) * c2 - 0.0057946086 * pow(dq1, 2) * pow(c1, 2) * pow(c2, 2) + 0.0057946086 * pow(dq2, 2) * pow(s1, 2) * pow(s2, 2) + 0.0192846035 * pow(dq2, 2) * pow(s1, 2) * s2 * c2 - 0.0057946086 * pow(dq2, 2) * pow(s1, 2) * pow(c2, 2) + 0.0057946086 * pow(dq2, 2) * pow(s2, 2) * pow(c1, 2) - 0.0057946086 * pow(dq2, 2) * pow(s2, 2) + 0.0192846035 * pow(dq2, 2) * s2 * pow(c1, 2) * c2 - 0.0192846035 * pow(dq2, 2) * s2 * c2 - 0.0057946086 * pow(dq2, 2) * pow(c1, 2) * pow(c2, 2) + 0.0057946086 * pow(dq2, 2) * pow(c2, 2) - 1.4179374 * s2 * c1 + 0.3932829 * c1 * c2, 1.0 * pow(dq1, 2) * s2 * c2, 1.0 * ddq1 * c2, -1.0 * pow(dq1, 2) * pow(s2, 2) + 1.0 * pow(dq1, 2) * pow(c2, 2), 1.0 * ddq2, -1.0 * ddq1 * s2, -1.0 * pow(dq1, 2) * s2 * c2, -2.0 * ddq2 * q3 * pow(s1, 2) * s2 * c2 - 2.0 * ddq2 * q3 * s2 * pow(c1, 2) * c2 + 2.0 * ddq2 * q3 * s2 * c2 + 0.8636 * ddq2 * pow(s1, 2) * s2 * c2 + 0.8636 * ddq2 * s2 * pow(c1, 2) * c2 - 0.8636 * ddq2 * s2 * c2 + 1.0 * ddq3 * pow(s1, 2) * pow(c2, 2) + 1.0 * ddq3 * pow(s2, 2) + 1.0 * ddq3 * pow(c1, 2) * pow(c2, 2) + 1.0 * pow(dq1, 2) * q3 * pow(s1, 2) * pow(s2, 2) - 1.0 * pow(dq1, 2) * q3 * pow(s1, 2) * pow(c2, 2) + 1.0 * pow(dq1, 2) * q3 * pow(s2, 2) * pow(c1, 2) - 1.0 * pow(dq1, 2) * q3 * pow(c1, 2) * pow(c2, 2) - 0.4318 * pow(dq1, 2) * pow(s1, 2) * pow(s2, 2) + 0.4318 * pow(dq1, 2) * pow(s1, 2) * pow(c2, 2) - 0.4318 * pow(dq1, 2) * pow(s2, 2) * pow(c1, 2) + 0.4318 * pow(dq1, 2) * pow(c1, 2) * pow(c2, 2) + 1.0 * pow(dq2, 2) * q3 * pow(s1, 2) * pow(s2, 2) - 1.0 * pow(dq2, 2) * q3 * pow(s1, 2) * pow(c2, 2) + 1.0 * pow(dq2, 2) * q3 * pow(s2, 2) * pow(c1, 2) - 1.0 * pow(dq2, 2) * q3 * pow(s2, 2) - 1.0 * pow(dq2, 2) * q3 * pow(c1, 2) * pow(c2, 2) + 1.0 * pow(dq2, 2) * q3 * pow(c2, 2) - 0.4318 * pow(dq2, 2) * pow(s1, 2) * pow(s2, 2) + 0.4318 * pow(dq2, 2) * pow(s1, 2) * pow(c2, 2) - 0.4318 * pow(dq2, 2) * pow(s2, 2) * pow(c1, 2) + 0.4318 * pow(dq2, 2) * pow(s2, 2) + 0.4318 * pow(dq2, 2) * pow(c1, 2) * pow(c2, 2) - 0.4318 * pow(dq2, 2) * pow(c2, 2) - 2.0 * dq2 * dq3 * pow(s1, 2) * s2 * c2 - 2.0 * dq2 * dq3 * s2 * pow(c1, 2) * c2 + 2.0 * dq2 * dq3 * s2 * c2 - 9.81 * pow(s2, 2) * c1 * c2 - 9.81 * c1 * pow(c2, 3), 1.0 * ddq1 * q3 * pow(s1, 2) * s2 + 1.0 * ddq1 * q3 * s2 * pow(c1, 2) - 0.4318 * ddq1 * pow(s1, 2) * s2 - 0.4318 * ddq1 * s2 * pow(c1, 2), 2.0 * ddq2 * q3 * pow(s1, 2) * pow(s2, 2) + 2.0 * ddq2 * q3 * pow(s2, 2) * pow(c1, 2) + 2.0 * ddq2 * q3 * pow(c2, 2) - 0.8636 * ddq2 * pow(s1, 2) * pow(s2, 2) - 0.8636 * ddq2 * pow(s2, 2) * pow(c1, 2) - 0.8636 * ddq2 * pow(c2, 2) - 1.0 * ddq3 * pow(s1, 2) * s2 * c2 - 1.0 * ddq3 * s2 * pow(c1, 2) * c2 + 1.0 * ddq3 * s2 * c2 + 2.0 * pow(dq1, 2) * q3 * pow(s1, 2) * s2 * c2 + 2.0 * pow(dq1, 2) * q3 * s2 * pow(c1, 2) * c2 - 0.8636 * pow(dq1, 2) * pow(s1, 2) * s2 * c2 - 0.8636 * pow(dq1, 2) * s2 * pow(c1, 2) * c2 + 2.0 * pow(dq2, 2) * q3 * pow(s1, 2) * s2 * c2 + 2.0 * pow(dq2, 2) * q3 * s2 * pow(c1, 2) * c2 - 2.0 * pow(dq2, 2) * q3 * s2 * c2 - 0.8636 * pow(dq2, 2) * pow(s1, 2) * s2 * c2 - 0.8636 * pow(dq2, 2) * s2 * pow(c1, 2) * c2 + 0.8636 * pow(dq2, 2) * s2 * c2 + 2.0 * dq2 * dq3 * pow(s1, 2) * pow(s2, 2) + 2.0 * dq2 * dq3 * pow(s2, 2) * pow(c1, 2) + 2.0 * dq2 * dq3 * pow(c2, 2) + 9.81 * pow(s2, 3) * c1 + 9.81 * s2 * c1 * pow(c2, 2), 1.0 * ddq2 * pow(q3, 2) * pow(s1, 2) * pow(s2, 2) + 1.0 * ddq2 * pow(q3, 2) * pow(s2, 2) * pow(c1, 2) + 1.0 * ddq2 * pow(q3, 2) * pow(c2, 2) - 0.8636 * ddq2 * q3 * pow(s1, 2) * pow(s2, 2) - 0.8636 * ddq2 * q3 * pow(s2, 2) * pow(c1, 2) - 0.8636 * ddq2 * q3 * pow(c2, 2) + 0.18645124 * ddq2 * pow(s1, 2) * pow(s2, 2) + 0.18645124 * ddq2 * pow(s2, 2) * pow(c1, 2) + 0.18645124 * ddq2 * pow(c2, 2) - 1.0 * ddq3 * q3 * pow(s1, 2) * s2 * c2 - 1.0 * ddq3 * q3 * s2 * pow(c1, 2) * c2 + 1.0 * ddq3 * q3 * s2 * c2 + 0.4318 * ddq3 * pow(s1, 2) * s2 * c2 + 0.4318 * ddq3 * s2 * pow(c1, 2) * c2 - 0.4318 * ddq3 * s2 * c2 + 1.0 * pow(dq1, 2) * pow(q3, 2) * pow(s1, 2) * s2 * c2 + 1.0 * pow(dq1, 2) * pow(q3, 2) * s2 * pow(c1, 2) * c2 - 0.8636 * pow(dq1, 2) * q3 * pow(s1, 2) * s2 * c2 - 0.8636 * pow(dq1, 2) * q3 * s2 * pow(c1, 2) * c2 + 0.18645124 * pow(dq1, 2) * pow(s1, 2) * s2 * c2 + 0.18645124 * pow(dq1, 2) * s2 * pow(c1, 2) * c2 + 1.0 * pow(dq2, 2) * pow(q3, 2) * pow(s1, 2) * s2 * c2 + 1.0 * pow(dq2, 2) * pow(q3, 2) * s2 * pow(c1, 2) * c2 - 1.0 * pow(dq2, 2) * pow(q3, 2) * s2 * c2 - 0.8636 * pow(dq2, 2) * q3 * pow(s1, 2) * s2 * c2 - 0.8636 * pow(dq2, 2) * q3 * s2 * pow(c1, 2) * c2 + 0.8636 * pow(dq2, 2) * q3 * s2 * c2 + 0.18645124 * pow(dq2, 2) * pow(s1, 2) * s2 * c2 + 0.18645124 * pow(dq2, 2) * s2 * pow(c1, 2) * c2 - 0.18645124 * pow(dq2, 2) * s2 * c2 + 2.0 * dq2 * dq3 * q3 * pow(s1, 2) * pow(s2, 2) + 2.0 * dq2 * dq3 * q3 * pow(s2, 2) * pow(c1, 2) + 2.0 * dq2 * dq3 * q3 * pow(c2, 2) - 0.8636 * dq2 * dq3 * pow(s1, 2) * pow(s2, 2) - 0.8636 * dq2 * dq3 * pow(s2, 2) * pow(c1, 2) - 0.8636 * dq2 * dq3 * pow(c2, 2) + 9.81 * q3 * pow(s2, 3) * c1 + 9.81 * q3 * s2 * c1 * pow(c2, 2) - 2.8180206 * pow(s2, 3) * c1 - 0.3932829 * pow(s2, 2) * c1 * c2 - 2.8180206 * s2 * c1 * pow(c2, 2) - 1.4179374 * s2 * c1 - 0.3932829 * c1 * pow(c2, 3) + 0.3932829 * c1 * c2, 0, 0, 0, 1.0 * pow(dq1, 2) * s2 * c2, -1.0 * ddq1 * c2, -1.0 * pow(dq1, 2) * pow(s2, 2) + 1.0 * pow(dq1, 2) * pow(c2, 2), 1.0 * ddq2, 1.0 * ddq1 * s2, -1.0 * pow(dq1, 2) * s2 * c2, -2.0 * ddq2 * q3 * pow(s1, 2) * s2 * c2 - 2.0 * ddq2 * q3 * s2 * pow(c1, 2) * c2 + 2.0 * ddq2 * q3 * s2 * c2 + 0.08018 * ddq2 * pow(s1, 2) * pow(c2, 2) + 0.08018 * ddq2 * pow(s2, 2) + 0.08018 * ddq2 * pow(c1, 2) * pow(c2, 2) + 1.0 * ddq3 * pow(s1, 2) * pow(c2, 2) + 1.0 * ddq3 * pow(s2, 2) + 1.0 * ddq3 * pow(c1, 2) * pow(c2, 2) + 1.0 * pow(dq1, 2) * q3 * pow(s1, 2) * pow(s2, 2) - 1.0 * pow(dq1, 2) * q3 * pow(s1, 2) * pow(c2, 2) + 1.0 * pow(dq1, 2) * q3 * pow(s2, 2) * pow(c1, 2) - 1.0 * pow(dq1, 2) * q3 * pow(c1, 2) * pow(c2, 2) - 0.08018 * pow(dq1, 2) * pow(s1, 2) * s2 * c2 - 0.08018 * pow(dq1, 2) * s2 * pow(c1, 2) * c2 + 1.0 * pow(dq2, 2) * q3 * pow(s1, 2) * pow(s2, 2) - 1.0 * pow(dq2, 2) * q3 * pow(s1, 2) * pow(c2, 2) + 1.0 * pow(dq2, 2) * q3 * pow(s2, 2) * pow(c1, 2) - 1.0 * pow(dq2, 2) * q3 * pow(s2, 2) - 1.0 * pow(dq2, 2) * q3 * pow(c1, 2) * pow(c2, 2) + 1.0 * pow(dq2, 2) * q3 * pow(c2, 2) - 0.08018 * pow(dq2, 2) * pow(s1, 2) * s2 * c2 - 0.08018 * pow(dq2, 2) * s2 * pow(c1, 2) * c2 + 0.08018 * pow(dq2, 2) * s2 * c2 - 2.0 * dq2 * dq3 * pow(s1, 2) * s2 * c2 - 2.0 * dq2 * dq3 * s2 * pow(c1, 2) * c2 + 2.0 * dq2 * dq3 * s2 * c2 + 9.81 * c1 * c2, -1.0 * ddq1 * q3 * pow(s1, 2) * s2 - 1.0 * ddq1 * q3 * s2 * pow(c1, 2) + 0.04009 * ddq1 * pow(s1, 2) * c2 + 0.04009 * ddq1 * pow(c1, 2) * c2, 2.0 * ddq2 * q3 * pow(s1, 2) * pow(s2, 2) + 2.0 * ddq2 * q3 * pow(s2, 2) * pow(c1, 2) + 2.0 * ddq2 * q3 * pow(c2, 2) - 0.08018 * ddq2 * pow(s1, 2) * s2 * c2 - 0.08018 * ddq2 * s2 * pow(c1, 2) * c2 + 0.08018 * ddq2 * s2 * c2 - 1.0 * ddq3 * pow(s1, 2) * s2 * c2 - 1.0 * ddq3 * s2 * pow(c1, 2) * c2 + 1.0 * ddq3 * s2 * c2 + 2.0 * pow(dq1, 2) * q3 * pow(s1, 2) * s2 * c2 + 2.0 * pow(dq1, 2) * q3 * s2 * pow(c1, 2) * c2 + 0.04009 * pow(dq1, 2) * pow(s1, 2) * pow(s2, 2) - 0.04009 * pow(dq1, 2) * pow(s1, 2) * pow(c2, 2) + 0.04009 * pow(dq1, 2) * pow(s2, 2) * pow(c1, 2) - 0.04009 * pow(dq1, 2) * pow(c1, 2) * pow(c2, 2) + 2.0 * pow(dq2, 2) * q3 * pow(s1, 2) * s2 * c2 + 2.0 * pow(dq2, 2) * q3 * s2 * pow(c1, 2) * c2 - 2.0 * pow(dq2, 2) * q3 * s2 * c2 + 0.04009 * pow(dq2, 2) * pow(s1, 2) * pow(s2, 2) - 0.04009 * pow(dq2, 2) * pow(s1, 2) * pow(c2, 2) + 0.04009 * pow(dq2, 2) * pow(s2, 2) * pow(c1, 2) - 0.04009 * pow(dq2, 2) * pow(s2, 2) - 0.04009 * pow(dq2, 2) * pow(c1, 2) * pow(c2, 2) + 0.04009 * pow(dq2, 2) * pow(c2, 2) + 2.0 * dq2 * dq3 * pow(s1, 2) * pow(s2, 2) + 2.0 * dq2 * dq3 * pow(s2, 2) * pow(c1, 2) + 2.0 * dq2 * dq3 * pow(c2, 2) - 9.81 * s2 * c1, 1.0 * ddq2 * pow(q3, 2) * pow(s1, 2) * pow(s2, 2) + 1.0 * ddq2 * pow(q3, 2) * pow(s2, 2) * pow(c1, 2) + 1.0 * ddq2 * pow(q3, 2) * pow(c2, 2) - 0.08018 * ddq2 * q3 * pow(s1, 2) * s2 * c2 - 0.08018 * ddq2 * q3 * s2 * pow(c1, 2) * c2 + 0.08018 * ddq2 * q3 * s2 * c2 + 0.0016072081 * ddq2 * pow(s1, 2) * pow(c2, 2) + 0.0016072081 * ddq2 * pow(s2, 2) + 0.0016072081 * ddq2 * pow(c1, 2) * pow(c2, 2) - 1.0 * ddq3 * q3 * pow(s1, 2) * s2 * c2 - 1.0 * ddq3 * q3 * s2 * pow(c1, 2) * c2 + 1.0 * ddq3 * q3 * s2 * c2 + 0.04009 * ddq3 * pow(s1, 2) * pow(c2, 2) + 0.04009 * ddq3 * pow(s2, 2) + 0.04009 * ddq3 * pow(c1, 2) * pow(c2, 2) + 1.0 * pow(dq1, 2) * pow(q3, 2) * pow(s1, 2) * s2 * c2 + 1.0 * pow(dq1, 2) * pow(q3, 2) * s2 * pow(c1, 2) * c2 + 0.04009 * pow(dq1, 2) * q3 * pow(s1, 2) * pow(s2, 2) - 0.04009 * pow(dq1, 2) * q3 * pow(s1, 2) * pow(c2, 2) + 0.04009 * pow(dq1, 2) * q3 * pow(s2, 2) * pow(c1, 2) - 0.04009 * pow(dq1, 2) * q3 * pow(c1, 2) * pow(c2, 2) - 0.0016072081 * pow(dq1, 2) * pow(s1, 2) * s2 * c2 - 0.0016072081 * pow(dq1, 2) * s2 * pow(c1, 2) * c2 + 1.0 * pow(dq2, 2) * pow(q3, 2) * pow(s1, 2) * s2 * c2 + 1.0 * pow(dq2, 2) * pow(q3, 2) * s2 * pow(c1, 2) * c2 - 1.0 * pow(dq2, 2) * pow(q3, 2) * s2 * c2 + 0.04009 * pow(dq2, 2) * q3 * pow(s1, 2) * pow(s2, 2) - 0.04009 * pow(dq2, 2) * q3 * pow(s1, 2) * pow(c2, 2) + 0.04009 * pow(dq2, 2) * q3 * pow(s2, 2) * pow(c1, 2) - 0.04009 * pow(dq2, 2) * q3 * pow(s2, 2) - 0.04009 * pow(dq2, 2) * q3 * pow(c1, 2) * pow(c2, 2) + 0.04009 * pow(dq2, 2) * q3 * pow(c2, 2) - 0.0016072081 * pow(dq2, 2) * pow(s1, 2) * s2 * c2 - 0.0016072081 * pow(dq2, 2) * s2 * pow(c1, 2) * c2 + 0.0016072081 * pow(dq2, 2) * s2 * c2 + 2.0 * dq2 * dq3 * q3 * pow(s1, 2) * pow(s2, 2) + 2.0 * dq2 * dq3 * q3 * pow(s2, 2) * pow(c1, 2) + 2.0 * dq2 * dq3 * q3 * pow(c2, 2) - 0.08018 * dq2 * dq3 * pow(s1, 2) * s2 * c2 - 0.08018 * dq2 * dq3 * s2 * pow(c1, 2) * c2 + 0.08018 * dq2 * dq3 * s2 * c2 - 9.81 * q3 * s2 * c1 + 0.3932829 * c1 * c2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0 * ddq2 * pow(s1, 2) * pow(c2, 2) + 1.0 * ddq2 * pow(s2, 2) + 1.0 * ddq2 * pow(c1, 2) * pow(c2, 2) - 1.0 * pow(dq1, 2) * pow(s1, 2) * s2 * c2 - 1.0 * pow(dq1, 2) * s2 * pow(c1, 2) * c2 - 1.0 * pow(dq2, 2) * pow(s1, 2) * s2 * c2 - 1.0 * pow(dq2, 2) * s2 * pow(c1, 2) * c2 + 1.0 * pow(dq2, 2) * s2 * c2, -1.0 * ddq1 * pow(s1, 2) * c2 - 1.0 * ddq1 * pow(c1, 2) * c2, -1.0 * ddq2 * pow(s1, 2) * s2 * c2 - 1.0 * ddq2 * s2 * pow(c1, 2) * c2 + 1.0 * ddq2 * s2 * c2 - 1.0 * pow(dq1, 2) * pow(s1, 2) * pow(c2, 2) - 1.0 * pow(dq1, 2) * pow(c1, 2) * pow(c2, 2) - 1.0 * pow(dq2, 2) * pow(s1, 2) * pow(c2, 2) - 1.0 * pow(dq2, 2) * pow(s2, 2) - 1.0 * pow(dq2, 2) * pow(c1, 2) * pow(c2, 2), -1.0 * ddq2 * q3 * pow(s1, 2) * s2 * c2 - 1.0 * ddq2 * q3 * s2 * pow(c1, 2) * c2 + 1.0 * ddq2 * q3 * s2 * c2 + 0.4318 * ddq2 * pow(s1, 2) * s2 * c2 + 0.4318 * ddq2 * s2 * pow(c1, 2) * c2 - 0.4318 * ddq2 * s2 * c2 + 1.0 * ddq3 * pow(s1, 2) * pow(c2, 2) + 1.0 * ddq3 * pow(s2, 2) + 1.0 * ddq3 * pow(c1, 2) * pow(c2, 2) - 1.0 * pow(dq1, 2) * q3 * pow(s1, 2) * pow(c2, 2) - 1.0 * pow(dq1, 2) * q3 * pow(c1, 2) * pow(c2, 2) + 0.4318 * pow(dq1, 2) * pow(s1, 2) * pow(c2, 2) + 0.4318 * pow(dq1, 2) * pow(c1, 2) * pow(c2, 2) - 1.0 * pow(dq2, 2) * q3 * pow(s1, 2) * pow(c2, 2) - 1.0 * pow(dq2, 2) * q3 * pow(s2, 2) - 1.0 * pow(dq2, 2) * q3 * pow(c1, 2) * pow(c2, 2) + 0.4318 * pow(dq2, 2) * pow(s1, 2) * pow(c2, 2) + 0.4318 * pow(dq2, 2) * pow(s2, 2) + 0.4318 * pow(dq2, 2) * pow(c1, 2) * pow(c2, 2) - 2.0 * dq2 * dq3 * pow(s1, 2) * s2 * c2 - 2.0 * dq2 * dq3 * s2 * pow(c1, 2) * c2 + 2.0 * dq2 * dq3 * s2 * c2 - 9.81 * pow(s2, 2) * c1 * c2 - 9.81 * c1 * pow(c2, 3), sign(dq3), dq3, 1, 0, 0, 0, 0, 0, 0, 1.0 * ddq2 * pow(s1, 2) * pow(c2, 2) + 1.0 * ddq2 * pow(s2, 2) + 1.0 * ddq2 * pow(c1, 2) * pow(c2, 2) - 1.0 * pow(dq1, 2) * pow(s1, 2) * s2 * c2 - 1.0 * pow(dq1, 2) * s2 * pow(c1, 2) * c2 - 1.0 * pow(dq2, 2) * pow(s1, 2) * s2 * c2 - 1.0 * pow(dq2, 2) * s2 * pow(c1, 2) * c2 + 1.0 * pow(dq2, 2) * s2 * c2, 1.0 * ddq1 * pow(s1, 2) * c2 + 1.0 * ddq1 * pow(c1, 2) * c2, -1.0 * ddq2 * pow(s1, 2) * s2 * c2 - 1.0 * ddq2 * s2 * pow(c1, 2) * c2 + 1.0 * ddq2 * s2 * c2 - 1.0 * pow(dq1, 2) * pow(s1, 2) * pow(c2, 2) - 1.0 * pow(dq1, 2) * pow(c1, 2) * pow(c2, 2) - 1.0 * pow(dq2, 2) * pow(s1, 2) * pow(c2, 2) - 1.0 * pow(dq2, 2) * pow(s2, 2) - 1.0 * pow(dq2, 2) * pow(c1, 2) * pow(c2, 2), -1.0 * ddq2 * q3 * pow(s1, 2) * s2 * c2 - 1.0 * ddq2 * q3 * s2 * pow(c1, 2) * c2 + 1.0 * ddq2 * q3 * s2 * c2 + 0.04009 * ddq2 * pow(s1, 2) * pow(c2, 2) + 0.04009 * ddq2 * pow(s2, 2) + 0.04009 * ddq2 * pow(c1, 2) * pow(c2, 2) + 1.0 * ddq3 * pow(s1, 2) * pow(c2, 2) + 1.0 * ddq3 * pow(s2, 2) + 1.0 * ddq3 * pow(c1, 2) * pow(c2, 2) - 1.0 * pow(dq1, 2) * q3 * pow(s1, 2) * pow(c2, 2) - 1.0 * pow(dq1, 2) * q3 * pow(c1, 2) * pow(c2, 2) - 0.04009 * pow(dq1, 2) * pow(s1, 2) * s2 * c2 - 0.04009 * pow(dq1, 2) * s2 * pow(c1, 2) * c2 - 1.0 * pow(dq2, 2) * q3 * pow(s1, 2) * pow(c2, 2) - 1.0 * pow(dq2, 2) * q3 * pow(s2, 2) - 1.0 * pow(dq2, 2) * q3 * pow(c1, 2) * pow(c2, 2) - 0.04009 * pow(dq2, 2) * pow(s1, 2) * s2 * c2 - 0.04009 * pow(dq2, 2) * s2 * pow(c1, 2) * c2 + 0.04009 * pow(dq2, 2) * s2 * c2 - 2.0 * dq2 * dq3 * pow(s1, 2) * s2 * c2 - 2.0 * dq2 * dq3 * s2 * pow(c1, 2) * c2 + 2.0 * dq2 * dq3 * s2 * c2 + 9.81 * c1 * c2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, sign(dq4), dq4, 1, ddq4, q4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0186 * sign(dq5), 1.03754596 * dq5, 1.01860000000000, 1.0186 * ddq5, -0.8306 * sign(-0.8306 * dq5 + 1.2178 * dq6), -0.8306 * (-0.8306 * dq5 + 1.2178 * dq6), -0.830600000000000, -0.8306 * sign(-0.8306 * dq5 + 1.2178 * dq7), -0.8306 * (-0.8306 * dq5 + 1.2178 * dq7), -0.830600000000000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.2178 * sign(-0.8306 * dq5 + 1.2178 * dq6), 1.2178 * (-0.8306 * dq5 + 1.2178 * dq6), 1.21780000000000, 0, 0, 0, sign(dq6), dq6, 1, ddq6, 0, 0, 0, 0, -1.2178 * sign(-1.2178 * dq6 + 1.2178 * dq7), -1.2178 * (-1.2178 * dq6 + 1.2178 * dq7), -1.21780000000000,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.2178 * sign(-0.8306 * dq5 + 1.2178 * dq7), 1.2178 * (-0.8306 * dq5 + 1.2178 * dq7), 1.21780000000000, 0, 0, 0, 0, sign(dq7), dq7, 1, ddq7, 1.2178 * sign(-1.2178 * dq6 + 1.2178 * dq7), 1.2178 * (-1.2178 * dq6 + 1.2178 * dq7), 1.21780000000000;
	
	// Uncomment if you want to exclude the counterweight contribution in the predicted torque (for debug purposes only)
	//std::cout << "tau1 cw contribution = " << H.block<1, 10>(0, 79) * delta.block<10, 1>(79, 0) << std::endl;
	//std::cout << "tau2 cw contribution = " << H.block<1, 10>(1, 79) * delta.block<10, 1>(79, 0) << std::endl;
	//std::cout << "tau3 cw contribution = " << H.block<1, 10>(2, 79) * delta.block<10, 1>(79, 0) << std::endl;
	//std::cout << std::endl;
	//H.block<7, 10>(0, 79).setZero();

	tau = H * delta;

	return tau;
}