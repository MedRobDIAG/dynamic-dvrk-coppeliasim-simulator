// dvrk dyn model library header file
#include "dvrkDynamics.hpp"
#include "dvrkKinematics.hpp"
#include "Timer.hpp"
#include "utils.hpp"
#include <chrono>

#define ESCAPE 27

// Remote CoppeliaSim API functions
extern "C" {
#include "extApi.h"
}

// Integer flags specifying the type of control used for the robot motion
// Kinematic: computes and sends joint velocity inputs
// Dynamic: computes and sends joint torque inputs
enum CTRL_TYPE {UNDEFINED = -1, KINEMATIC, DYNAMIC};
enum TASK_TYPE {REGULATION_TASK = 1, LINEAR_TRAJECTORY_TASK, SPIRAL_TRAJECTORY_TASK, EXCITING_JOINT_TRAJECTORY_TASK};

/**
* @brief Get function
* Get the initial fixed transformation Tw0 expressing the pose of the base frame of the arm specified by the input string, wrt the CoppeliaSim world frame
* @param clientID: the remote CoppeliaSim ID
* @param ref0: string with the name of the dummy object representing the base frame of the considered arm
* @return the 4x4 homogeneous transformation with the requested pose
*/
Eigen::Matrix4f getTw0(const int& clientID, const std::string& ref0);


int main(int argc, char** argv) {

	std::cout << "Test program for the dll library of the DVRK PSM dynamic model" << std::endl;

	// Setup flags
	int choose;
	float ctrlType = CTRL_TYPE::UNDEFINED;
	bool saveLogs = false; //  Set this flag to true if you want to log data from the simulation
	bool sync = true; // synchronous flag for CoppeliaSim simulation
	bool runPlayback = true; // playback flag - optin valid only in case of joint excitation trajectories
	bool withCW = true; // Set this flag to true if you want to account the counterweight in the robot model
	bool withFriction = false; // Set this flag to true if you want to account the friction in the robot model 

	// Log-related variables
	int logNum = 0;
	char curdir[256];

	// Timer variables
	double tictoc, tic, toc;
	float t_curr, t_prev;
	Timer clock;

	// CoppeliaSim object name and signal lists
	std::string jointNames[PSM_FULL_JOINTS] = { "J1_PSM1","J2_PSM1","","J23_PSM1","J22_PSM1","J24_PSM1","J25_PSM1","J3_PSM1","","J4_PSM1","J5_PSM1","J6_PSM1","J7_PSM1","" }; // Joint names in CoppeliaSim simulation
	std::string jointTorqueSigNames[3] = { "Torque_J1","Torque_J2","Force_J3" };
	std::string jointMirror3Name = "J3_mirror_PSM1";
	std::string posSig[3] = { "x_ee","y_ee","z_ee" };
	std::string qdesSig[3] = { "qLdes_1", "qLdes_2", "qLdes_3" };
	std::string errPosSig[3] = { "error_x","error_y","error_z" };
	std::string errOriSig[3] = { "error_alfa","error_beta","error_gamma" };

	// Playback variables
	std::vector < Eigen::VectorXf > offData;
	std::string offlineDatasetPath;

	// Robot-related variables
	psmFullJointsVectori jointHandles;
	psmActiveJointsVectori jointActiveHandles;
	psmFullJointsVectorf psmJointPositions, psmJointMsrTorques;
	psmActiveJointsVectorf psmJointActivePositions, psmJointActivePrevPositions, psmJointActiveVelocities, psmJointActivePrevVelocities, psmJointActiveCmdVelocities, psmJointActiveMsrVelocities, psmJointActiveAccelerations;
	psmActiveJointsVectorf tauModel, tauMsr, tauCmd, tauMsr_old;
	psmActiveJointsVectorf g;
	Eigen::VectorXf qdState, qddState, qdot_f_prev, qddot_f_prev;
	Eigen::Array<float, PSM_ACTIVE_JOINTS, 1> activeJointsIdxs; // Mask of indexes to switch from full to active only joint vector
	activeJointsIdxs << 0, 1, 7, 9, 10, 11, 12;

	// Stringstream variables
	std::stringstream qmsrSS("");
	std::stringstream qdcmdSS("");
	std::stringstream qdmsrSS("");
	std::stringstream qddmsrSS("");
	std::stringstream posSS("");
	std::stringstream posRefSS("");
	std::stringstream oriSS("");
	std::stringstream oriRefSS("");
	std::stringstream tauModelSS("");
	std::stringstream tauMsrSS("");
	std::stringstream tauCmdSS("");
	std::stringstream gSS("");

	// Reference trajectory variables
	Eigen::VectorXf refTrajectoryVisual;
	Eigen::VectorXf draw_RefTrajectory_x, draw_RefTrajectory_y, draw_RefTrajectory_z;
	Eigen::VectorXf refTrajectory_x, refTrajectory_y, refTrajectory_z;
	Eigen::VectorXf refTrajectoryDiff_x, refTrajectoryDiff_y, refTrajectoryDiff_z;
	Eigen::VectorXf refTrajectoryDiffDiff_x, refTrajectoryDiffDiff_y, refTrajectoryDiffDiff_z;
	Eigen::VectorXf refTrajectory_q1, refTrajectory_q2, refTrajectory_q3;
	Eigen::VectorXf refTrajectoryDiff_q1, refTrajectoryDiff_q2, refTrajectoryDiff_q3;
	Eigen::MatrixXf refTrajectory_q;
	Eigen::MatrixXf refTrajectoryDiff_q;
	Eigen::MatrixXf refTrajectoryDiffDiff_q;
	Eigen::Vector3f startingPoint;
	Eigen::Vector3f refOrientation;
	Eigen::Vector3f pdes;
	Eigen::VectorXf qdes, qd_des, qdd_des;
	Eigen::Vector3f pdot_des;
	Eigen::Vector3f pdotdot_des;
	Eigen::Matrix3f Rdes;
	Eigen::Quaternionf quatdes;
	Eigen::MatrixXf qExcitParams;
	double t = 0.0;
	double trajDuration = 10.0; // [s]
	double Tsim = 0.001; // [s]
	int totSteps;
	int startHandle; // Handle of dummy in simulated scene specifying the starting point of the trajectory
	int iterations;

	// Robot and control variables
	Eigen::Matrix4f Tee;
	Eigen::Matrix3f Ree;
	Eigen::Vector3f pee, rpy_ee;
	Eigen::Vector3f err_p, errI_p, errD_p;
	Eigen::VectorXf errP_q, errD_q;
	Eigen::Quaternionf quat_ee, err_quat;
	Eigen::Matrix<float, 6, PSM_ACTIVE_JOINTS - 1> J;
	Eigen::Matrix<float, 3, 3> Jl, Jo, Jlprev, Jldot;
	Eigen::Matrix<float, PSM_ACTIVE_JOINTS, 6> JT, Jpinv;
	Eigen::Matrix<float, 3, 3> JlT, Jlpinv;
	Eigen::Vector6f u_cart;
	Eigen::Vector3f u_pos, u_ori;
	Eigen::Matrix3f Kpp, Kpd, Kpi, Kop;
	Eigen::MatrixXf Kqp, Kqd;
	psmActiveJointsVectorf coriolis;
	psmMassMatrixf M;

	// Variables initialization
	Kpp.setZero();
	Kpd.setZero();
	Kpi.setZero();
	Kop.setZero();
	errI_p.setZero();
	J.setZero();
	Jl.setZero();
	Jlprev.setZero();
	JT.setZero();
	JlT.setZero();
	Jpinv.setZero();
	Jlpinv.setZero();
	Jldot.setZero();
	pdot_des.setZero();
	pdotdot_des.setZero();
	qddState.setZero(PSM_ACTIVE_JOINTS);
	qdState.setZero(PSM_ACTIVE_JOINTS);
	psmJointActiveCmdVelocities.setZero();
	psmJointActivePositions.setZero();
	psmJointActivePrevPositions.setZero();
	psmJointActiveVelocities.setZero();
	psmJointActivePrevVelocities.setZero();
	psmJointActiveAccelerations.setZero();
	psmJointMsrTorques.setZero();
	tauModel.setZero();
	tauCmd.setZero();
	jointHandles.setZero();
	psmJointPositions.setZero();
	qdot_f_prev.setZero(PSM_ACTIVE_JOINTS);
	qddot_f_prev.setZero(PSM_ACTIVE_JOINTS);
	errP_q.setZero(PSM_ACTIVE_JOINTS);
	errD_q.setZero(PSM_ACTIVE_JOINTS);
	qdes.setZero(PSM_ACTIVE_JOINTS);
	qd_des.setZero(PSM_ACTIVE_JOINTS);
	qdd_des.setZero(PSM_ACTIVE_JOINTS);
	Kqp.setZero(PSM_ACTIVE_JOINTS, PSM_ACTIVE_JOINTS);
	Kqd.setZero(PSM_ACTIVE_JOINTS, PSM_ACTIVE_JOINTS);

	// Get current directory
	GetCurrentDirectoryA(256, curdir);

	// Close previously opened remote communications with CoppeliaSim and instantiate a new one
	simxFinish(-1);

	//Create and start a connection with the CoppeliaSim at default port
	int clientID = simxStart("127.0.0.1", 19997, true, true, 5000, 5);
	if (clientID > -1) {
		std::cout << "CoppeliaSim remote connection successfully connected." << std::endl;
	}
	else {
		return -1;
	}

	// Get robot base frame CoppeliaSim handle
	simxInt psm1BaseHandle;
	simxInt dummyHandle;
	simxGetObjectHandle(clientID, "DH0_ref", &psm1BaseHandle, simx_opmode_blocking);
	simxGetObjectHandle(clientID, "eeDummy", &dummyHandle, simx_opmode_blocking);

	// Get initial joint configuration
	for (int i = 0; i < PSM_FULL_JOINTS; i++) {
		simxGetObjectHandle(clientID, jointNames[i].c_str(), &jointHandles[i], simx_opmode_blocking);
		simxGetJointPosition(clientID, jointHandles[i], &psmJointPositions(i), simx_opmode_blocking);
		simxGetJointForce(clientID, jointHandles[i], &psmJointMsrTorques(i), simx_opmode_blocking);
	}

	// Extract active joint position
	jointActiveHandles = (activeJointsIdxs.unaryExpr(jointHandles)).matrix().cast<int>();
	psmJointActivePositions = (activeJointsIdxs.unaryExpr(psmJointPositions)).matrix();
	tauMsr = (activeJointsIdxs.unaryExpr(psmJointMsrTorques)).matrix();
	psmJointActivePrevPositions = psmJointActivePositions;
	tauMsr_old = tauMsr;

	// Get the reference trajectory starting point coordinates
	simxGetObjectHandle(clientID, "Start", &startHandle, simx_opmode_blocking);
	simxGetObjectPosition(clientID, startHandle, -1, &startingPoint[0], simx_opmode_blocking);

	// Let the user choose which trajectory should be executed autonomously by the robot 
	std::cout << "Select task type: " <<
		"\n\t1. regulation task" <<
		"\n\t2. tracking task -- linear trajectory" <<
		"\n\t3. tracking task -- spiral trajectory" <<
		"\n\t4. joint excitation trajectory (work with dt=5.0ms for consistency with offline database)" << std::endl;
	std::cin >> choose;

	if (choose == TASK_TYPE::EXCITING_JOINT_TRAJECTORY_TASK) {
		// Force runPlayback flag to true
		runPlayback = true;
	}
	else {
		runPlayback = false;
	}

	// Set CoppeliaSim synchronous mode
	simxSynchronous(clientID, sync);

	// Start the CoppeliaSim simulation remotely
	simxStartSimulation(clientID, simx_opmode_blocking);

	if (sync) {
		simxSynchronousTrigger(clientID);
	}

	// Get simtulation time step
	Tsim = -1.0;
	simxFloat Tsim_f = -1.0;
	while (Tsim == -1.0) {
		simxGetFloatSignal(clientID, "simulationTimeStep", &Tsim_f, simx_opmode_blocking);
		Tsim = Tsim_f;
		std::cout << "Retrieving simulation time step ... " << std::endl;
	}
	std::cout << "Simulation Time step = " << Tsim << " [s]" << std::endl;

	// Initialize the reference trajectories for tracking and visualization
	totSteps = int(trajDuration / Tsim);
	iterations = totSteps;
	int num_microstep = 50; // Set to 1 to make this step uneffective
	double Ts = Tsim / num_microstep;
	int bufferVisualSize = trajDuration / Ts;
	int bufferSize = int(trajDuration / Tsim);
	int xyzBufferVisualSize = 3 * bufferVisualSize;

	// Initialization of Cartesian reference visual trajectory
	refTrajectoryVisual.setZero(xyzBufferVisualSize);
	draw_RefTrajectory_x.setZero(bufferVisualSize);
	draw_RefTrajectory_y.setZero(bufferVisualSize);
	draw_RefTrajectory_z.setZero(bufferVisualSize);

	// Initialization of Cartesian reference trajetory and derivatives
	refTrajectory_x.setZero(bufferSize);
	refTrajectory_y.setZero(bufferSize);
	refTrajectory_z.setZero(bufferSize);
	refTrajectoryDiff_x.setZero(bufferSize);
	refTrajectoryDiff_y.setZero(bufferSize);
	refTrajectoryDiff_z.setZero(bufferSize);
	refTrajectoryDiffDiff_x.setZero(bufferSize);
	refTrajectoryDiffDiff_y.setZero(bufferSize);
	refTrajectoryDiffDiff_z.setZero(bufferSize);
	refTrajectory_q1.setZero(bufferSize);
	refTrajectory_q2.setZero(bufferSize);
	refTrajectory_q3.setZero(bufferSize);
	refTrajectoryDiff_q1.setZero(bufferSize);
	refTrajectoryDiff_q2.setZero(bufferSize);
	refTrajectoryDiff_q3.setZero(bufferSize);
	refTrajectory_q.setZero(PSM_ACTIVE_JOINTS, bufferSize);
	refTrajectoryDiff_q.setZero(PSM_ACTIVE_JOINTS, bufferSize);
	refTrajectoryDiffDiff_q.setZero(PSM_ACTIVE_JOINTS, bufferSize);


	if (runPlayback) {

		// Get offline data if playback run is set (TODO: Check if can be kept)
		if (choose == TASK_TYPE::EXCITING_JOINT_TRAJECTORY_TASK) {
			offlineDatasetPath = std::string(curdir) + "\\one_no_weight_results.csv";
		
			offData = loadCSVLog(offlineDatasetPath);
			iterations = offData.size();
			std::cout << "Offline data loaded. Size = " << offData.size() << std::endl;

			// Get the first joint configuration at time 0 and set it on the robot in the scene
			psmJointActivePositions = offData[0].head(PSM_ACTIVE_JOINTS);
			tauMsr = offData[0].segment(PSM_ACTIVE_JOINTS, 2 * PSM_ACTIVE_JOINTS - 1);

		}

	}
	else {

		// Popoulate the samples of the reference trajectory
		if (choose == TASK_TYPE::REGULATION_TASK) { // Regulation task
			simxInt eeRef;
			simxGetObjectHandle(clientID, "eeRef_regulation", &eeRef, simx_opmode_blocking);
			simxGetObjectPosition(clientID, eeRef, -1, &pdes[0], simx_opmode_blocking);
			std::cout << "pdes = " << pdes.transpose() << std::endl;

		}
		else if (choose == TASK_TYPE::LINEAR_TRAJECTORY_TASK) { // Linear trajectory

			std::cout << "\n\nA straigth line trajectory has been selected, continue on Coppelia Sim\n";
			Eigen::Vector3d angular_coef(-0.005, -0.005, 0.005);

			// Trajectory buffer for visualization (higher sample time Ts)
			double ix = 0, iy = 0, iz = 0;
			for (double i = 0; i < trajDuration - Ts; i = i + Ts) { draw_RefTrajectory_x(ix) = angular_coef(0) * i + startingPoint(0); ix++; }
			for (double i = 0; i < trajDuration - Ts; i = i + Ts) { draw_RefTrajectory_y(iy) = angular_coef(1) * i + startingPoint(1); iy++; }
			for (double i = 0; i < trajDuration - Ts; i = i + Ts) { draw_RefTrajectory_z(iz) = angular_coef(2) * i + startingPoint(2); iz++; }

			// Trajectory buffer for control (sample time Tsim)
			std::cout << "startingPoint = " << startingPoint.transpose() << std::endl;
			ix = 0; iy = 0; iz = 0;
			for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectory_x(ix) = angular_coef(0) * i + startingPoint(0); ix++; }
			for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectory_y(iy) = angular_coef(1) * i + startingPoint(1); iy++; }
			for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectory_z(iz) = angular_coef(2) * i + startingPoint(2); iz++; }

			// Trajectory derivative buffer for control (sample time Tsim)
			ix = 0; iy = 0; iz = 0;
			for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectoryDiff_x(ix) = angular_coef(0); ix++; }
			for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectoryDiff_y(iy) = angular_coef(1); iy++; }
			for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectoryDiff_z(iz) = angular_coef(2); iz++; }

			// Trajectory second derivative buffer for control (sample time Tsim)
			ix = 0; iy = 0; iz = 0;
			for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectoryDiffDiff_x(ix) = 0.0; ix++; }
			for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectoryDiffDiff_y(iy) = 0.0; iy++; }
			for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectoryDiffDiff_z(iz) = 0.0; iz++; }

		}
		else if (choose == TASK_TYPE::SPIRAL_TRAJECTORY_TASK) { // Spiral trajectory

			std::cout << "\n\nA spyral trajectory has been selected, continue on Coppelia Sim\n";
			Eigen::Vector3d Amplitude(0.045, 0.045, 0.015);
			Eigen::Vector3d Phase(0, 0, 0);
			Eigen::Vector3d Frequency(0.75, 0.75, 0.1);

			// Trajectory buffer for visualization (higher sample time Ts)
			double ix = 0, iy = 0, iz = 0;
			for (double i = 0; i < trajDuration - Ts; i = i + Ts) { draw_RefTrajectory_x(ix) = Amplitude(0) * cos(2.0 * M_PI * Frequency(0) * i + Phase(0)) + startingPoint(0); ix++; }
			for (double i = 0; i < trajDuration - Ts; i = i + Ts) { draw_RefTrajectory_y(iy) = Amplitude(1) * sin(2.0 * M_PI * Frequency(1) * i + Phase(1)) + startingPoint(1); iy++; }
			for (double i = 0; i < trajDuration - Ts; i = i + Ts) { draw_RefTrajectory_z(iz) = Amplitude(2) * cos(2.0 * M_PI * Frequency(2) * i + Phase(2)) + startingPoint(2); iz++; }

			// Trajectory buffer for control (sample time Tsim)
			ix = 0; iy = 0; iz = 0;
			for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectory_x(ix) = Amplitude(0) * cos(2.0 * M_PI * Frequency(0) * i + Phase(0)) + startingPoint(0); ix++; }
			for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectory_y(iy) = Amplitude(1) * sin(2.0 * M_PI * Frequency(1) * i + Phase(1)) + startingPoint(1); iy++; }
			for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectory_z(iz) = Amplitude(2) * cos(2.0 * M_PI * Frequency(2) * i + Phase(2)) + startingPoint(2); iz++; }

			// Trajectory derivative buffer for control (sample time Tsim)
			ix = 0; iy = 0; iz = 0;
			for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectoryDiff_x(ix) = -2.0 * M_PI * Frequency(0) * Amplitude(0) * sin(2.0 * M_PI * Frequency(0) * i + Phase(0)); ix++; }
			for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectoryDiff_y(iy) = 2.0 * M_PI * Frequency(1) * Amplitude(1) * cos(2.0 * M_PI * Frequency(1) * i + Phase(1)); iy++; }
			for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectoryDiff_z(iz) = -2.0 * M_PI * Frequency(2) * Amplitude(2) * sin(2.0 * M_PI * Frequency(2) * i + Phase(2)); iz++; }

			// Trajectory second derivative buffer for control (sample time Tsim)
			ix = 0; iy = 0; iz = 0;
			for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectoryDiffDiff_x(ix) = -pow(2.0 * M_PI * Frequency(0), 2) * Amplitude(0) * cos(2.0 * M_PI * Frequency(0) * i + Phase(0)); ix++; }
			for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectoryDiffDiff_y(iy) = -pow(2.0 * M_PI * Frequency(1), 2) * Amplitude(1) * sin(2.0 * M_PI * Frequency(1) * i + Phase(1)); iy++; }
			for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectoryDiffDiff_z(iz) = -pow(2.0 * M_PI * Frequency(2), 2) * Amplitude(2) * cos(2.0 * M_PI * Frequency(2) * i + Phase(2)); iz++; }

		}
		else if (choose == TASK_TYPE::EXCITING_JOINT_TRAJECTORY_TASK) { // Exciting joint trajectory

			double ff = 0.18;
			// q0 a1 a2 a3 a4 a5 a6 b1 b2 b3 b4 b5 b6 (from Wang et al., 2019 - see: 
			qExcitParams.setZero(7, 13);
			qExcitParams << 0.1877, -0.0191, 0.0107, 0.226, 0.2893, 0.2369, 0.1862, 0.0109, -0.0813, -0.116, -0.0245, 0.1865, -1.0593,
				0.0864, 0.1467, 0.47, 0.6868, -0.1115, 0.233, 0.434, 0.0722, 1.2791, 0.0129, -0.0364, -0.1937, 0.2359,
				0.1486, 0.0464, 0.0323, 0.0296, 0.0491, 0.1968, 0.2014, 0.0307, 0.0485, 0.0152, 0.218, 0.0493, 0.1782,
				0.2258, 0.155, 0.3255, -0.0803, 0.3699, -0.0934, 0.4492, 0.2419, 1.0061, 0.9141, -0.7797, 0.5985, 0.0095,
				0.3635, -0.1258, 0.8912, 0.0766, -0.2848, 0.4034, 0.3895, 0.7037, 1.1342, 0.2914, 0.2584, -0.4521, 0.4273,
				-0.3571, 0.4269, 0.4307, -0.2174, -0.356, 0.7692, 0.4352, 0.4253, 1.1555, 0.8344, 0.3126, 0.0576, 0.4472,
				0.3501, 0.6519, 0.4726, -0.1491, 0.5939, 0.5882, -0.1018, 0.4729, 1.3575, 0.1049, 0.3389, -0.5517, 0.4689;

			// Trajectory buffer for control (sample time Tsim)
			//double iq[] = {0, 0, 0};
			for (int k = 0; k < 7; k++) {
				int iq = 0;
				for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) {
					refTrajectory_q(k, iq) = qExcitParams(k, 0);

					for (int j = 1; j < 6; j++) {
						refTrajectory_q(k, iq) += qExcitParams(k, j) / (2.0 * M_PI * ff * j) * sin(2.0 * M_PI * ff * j * i) - qExcitParams(k, 6 + j) / (2.0 * M_PI * ff * j) * cos(2.0 * M_PI * ff * j * i);
						refTrajectoryDiff_q(k, iq) += qExcitParams(k, j) * cos(2.0 * M_PI * ff * j * i) + qExcitParams(k, 6 + j) * sin(2.0 * M_PI * ff * j * i);
						refTrajectoryDiffDiff_q(k, iq) += -qExcitParams(k, j) * 2.0 * M_PI * ff * j * sin(2.0 * M_PI * ff * j * i) + qExcitParams(k, 6 + j) * 2.0 * M_PI * ff * j * cos(2.0 * M_PI * ff * j * i);
					}

					iq++;
				}
			}

		}
	}


	// Collect the trajectory samples for visualization
	refTrajectoryVisual << draw_RefTrajectory_x, draw_RefTrajectory_y, draw_RefTrajectory_z;

	// Define a reference desired attitude for the PSM gripper
	refOrientation << 0.0, 0.0, 0.0;
	Rdes = Eigen::Matrix3f(Eigen::AngleAxisf(refOrientation[2], Eigen::Vector3f::UnitZ())
		* Eigen::AngleAxisf(refOrientation[1], Eigen::Vector3f::UnitY())
		* Eigen::AngleAxisf(refOrientation[0], Eigen::Vector3f::UnitX()));
	quatdes = Eigen::Quaternionf(Rdes.transpose());

	// Get the fixed initial transformation (Tw0)
	Eigen::Matrix4f Tw0 = getTw0(clientID, "DH0_ref");

	// Build the direct kinematics of the PSM arm, according to the modified DH convention with simplifications
	Eigen::Matrix<float, PSM_ACTIVE_JOINTS, 5> DH = psmSimplifiedModDH(psmJointPositions);

	// Build the corresponding chain of homogeneous transformation and the resulting direct kinematics matrix
	std::pair < std::vector<Eigen::Matrix4f>, Eigen::Matrix<float, 6, PSM_ACTIVE_JOINTS> > psmKin = psmKinematics(psmJointPositions, Tw0);

	// Draw the reference trajectory on CoppeliaSim through remote call of script function
	simxCallScriptFunction(clientID, "ApiRemoteBridge", sim_scripttype_childscript, "drawTrajectory", 1, &bufferVisualSize, xyzBufferVisualSize, &refTrajectoryVisual(0), 0, nullptr, 0, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, simx_opmode_blocking);

	// Initialize streaming for joint position readings
	for (int i = 0; i < PSM_FULL_JOINTS; i++) {
		simxGetJointPosition(clientID, jointHandles[i], &psmJointPositions[i], simx_opmode_streaming);
		simxGetJointForce(clientID, jointHandles[i], &psmJointMsrTorques[i], simx_opmode_streaming);
	}
	psmJointActivePositions = (activeJointsIdxs.unaryExpr(psmJointPositions)).matrix();
	tauMsr = (activeJointsIdxs.unaryExpr(psmJointMsrTorques)).matrix();

	// Forward the type of chosen: retrieve the flag to select kinematic or dynamic simulation
	while (ctrlType == CTRL_TYPE::UNDEFINED) {
		simxGetFloatSignal(clientID, "dynamicCtrlFlagSignal", &ctrlType, simx_opmode_blocking);

		// Trigger the next simulation step on CoppeliaSim (if synchronous mode is enabled)
		if (sync) {
			simxSynchronousTrigger(clientID);
		}
	}
	std::cout << "ctrlType = " << ctrlType << std::endl;

	// Set control gains
	if (choose == TASK_TYPE::REGULATION_TASK) {
		if (ctrlType == CTRL_TYPE::KINEMATIC) {

			Kpp.diagonal() << 5.0, 5.0, 5.0;
			Kop.diagonal() << 5.0, 5.0, 5.0;
			Kpd.diagonal() << 0.0, 0.0, 0.0;
			Kpi.diagonal() << 0.1, 0.1, 0.1;

		}
		else if (ctrlType == CTRL_TYPE::DYNAMIC) {

			Kpp.diagonal() << 100.0, 100.0, 100.0;
			Kpd.diagonal() = 2.0 * Kpp.diagonal().cwiseSqrt();

		}
	}
	else if (choose == TASK_TYPE::LINEAR_TRAJECTORY_TASK || choose == TASK_TYPE::SPIRAL_TRAJECTORY_TASK) {

		if (ctrlType == CTRL_TYPE::KINEMATIC) {
			Kpp.diagonal() << 15.0, 15.0, 15.0;
			Kop.diagonal() << 15.0, 15.0, 15.0;
			Kpd.diagonal() << 0.0, 0.0, 0.0;
			Kpi.diagonal() << 0.0005, 0.0005, 0.0005;
		}
		else if (ctrlType == CTRL_TYPE::DYNAMIC) {
			Kpp.diagonal() << 250.0, 250.0, 500.0; // 10ms
			Kpd.diagonal() << 31.0, 31.0, 20.0;
		}

	}
	else if (choose == TASK_TYPE::EXCITING_JOINT_TRAJECTORY_TASK) {

		if (ctrlType == CTRL_TYPE::KINEMATIC) {
			Kqp.diagonal() << 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0;
			Kqd.diagonal() << 30.0, 30.0, 10.0, 1.0, 1.0, 1.0, 1.0;

		}
		else if (ctrlType == CTRL_TYPE::DYNAMIC) {
			errP_q.setZero(3);
			errD_q.setZero(3);
			qdes.setZero(3);
			qd_des.setZero(3);
			qdd_des.setZero(3);
			Kqp.setZero(3, 3);
			Kqd.setZero(3, 3);
			Kqp.diagonal() << 200.0, 200.0, 100.0;
			Kqd.diagonal() << 30.0, 30.0, 10.0;

		}

	}

	// ---  Start control loop
	double excitTrajRate = 20.0;
	double excitTrajTime = 0.0;
	int excitCurIdx = 0.0;
	int sampleRep = floor((1.0/Tsim)/ excitTrajRate);
	double alpha = 0.1;

	for (int i = 0; i < iterations; i++) {

		// Get the current joint readings
		for (int k = 0; k < PSM_FULL_JOINTS; k++) {
			simxGetJointPosition(clientID, jointHandles[k], &psmJointPositions[k], simx_opmode_buffer);
			simxGetJointForce(clientID, jointHandles[k], &psmJointMsrTorques[k], simx_opmode_buffer);
		}

		// Extract active joint position and torques
		psmJointActivePositions = (activeJointsIdxs.unaryExpr(psmJointPositions)).matrix();
		tauMsr = (activeJointsIdxs.unaryExpr(psmJointMsrTorques)).matrix();

		// Filter measured torques
		//tauMsr = alpha * tauMsr + (1 - alpha) * tauMsr_old;
		tauMsr_old = tauMsr;

		// Online simulation
		if (!runPlayback) {

			// Compute measured joint velocities
			psmJointActiveMsrVelocities = (psmJointActivePositions - psmJointActivePrevPositions) / Tsim;

			// Filter joint velocity
			psmJointActiveMsrVelocities = filter_velocity(psmJointActiveMsrVelocities, qdState);
			//psmJointActiveMsrVelocities = alpha * psmJointActiveMsrVelocities + (1 - alpha) * qdot_f_prev;
			qdot_f_prev = psmJointActiveMsrVelocities;
			psmJointActiveVelocities = psmJointActiveMsrVelocities;

			// Get the current dVRK-PSM EE pose from the measured joint configuration
			psmKin.first.clear();
			psmKin = psmKinematics(psmJointPositions, Tw0);
			Tee = psmKin.first[PSM_ACTIVE_JOINTS - 1];
			Ree = Tee.topLeftCorner(3, 3);
			pee = Tee.block<3, 1>(0, 3);
			rpy_ee = rot2rpy(Ree);
			quat_ee = Eigen::Quaternionf(Ree.transpose());

			// Debug: Set pose of the dummy to verify the resulting current EE pose
			simxSetObjectPosition(clientID, dummyHandle, -1, &pee[0], simx_opmode_oneshot);
			simxSetObjectQuaternion(clientID, dummyHandle, -1, &quat_ee.coeffs()[0], simx_opmode_oneshot);

			// Get the current Jacobian matrix and compute the pseudo-inverse
			J = psmKin.second.leftCols(6);
			Jl = J.topLeftCorner(3, 3);
			JlT = Jl.transpose();
			Jldot = (Jl - Jlprev) / Tsim;
			Jlprev = Jl;
			Jlpinv = Jl.inverse();

			// Compute the position error
			if (choose != TASK_TYPE::REGULATION_TASK && choose != TASK_TYPE::EXCITING_JOINT_TRAJECTORY_TASK) {
				pdes << refTrajectory_x(i), refTrajectory_y(i), refTrajectory_z(i);
				pdot_des << refTrajectoryDiff_x(i), refTrajectoryDiff_y(i), refTrajectoryDiff_z(i);
				pdotdot_des << refTrajectoryDiffDiff_x(i), refTrajectoryDiffDiff_y(i), refTrajectoryDiffDiff_z(i);
			}
			else if (choose == TASK_TYPE::EXCITING_JOINT_TRAJECTORY_TASK) {

				if (ctrlType == CTRL_TYPE::KINEMATIC) {
					qdes << refTrajectory_q.col(i);
					qd_des << refTrajectoryDiff_q.col(i);
					qdd_des << refTrajectoryDiffDiff_q.col(i);
					errP_q = qdes - psmJointActivePositions;
					errD_q = qd_des - psmJointActiveMsrVelocities;
				}
				else if (ctrlType == CTRL_TYPE::DYNAMIC) {
					qdes << refTrajectory_q.col(i).head(3);
					qd_des << refTrajectoryDiff_q.col(i).head(3);
					qdd_des << refTrajectoryDiffDiff_q.col(i).head(3);
					errP_q = qdes - psmJointActivePositions.head(3);
					errD_q = qd_des - psmJointActiveMsrVelocities.head(3);
				}
			}
			err_p = pdes - pee;
			errI_p += err_p * Tsim;

			// Compute the orientation error
			err_quat = quatdes * quat_ee;
			rpy_ee = rot2rpy(err_quat.toRotationMatrix());

			// Compute the Cartesian PID action
			u_pos = (Kpp * err_p + Kpd * errD_p + Kpi * errI_p) + pdot_des;
			u_ori = Kop * err_quat.vec() * err_quat.w();
			u_cart << u_pos, u_ori;

			// Compute the final control law depending on the chosen control type 
			if (ctrlType == CTRL_TYPE::KINEMATIC) {

				if (choose != TASK_TYPE::EXCITING_JOINT_TRAJECTORY_TASK) {
					psmJointActiveCmdVelocities.topRows(6) = J.inverse() * u_cart; // [velocities, rad/s]
				}
				else {
					psmJointActiveCmdVelocities = Kqp * errP_q + qd_des;
				}
				psmJointActiveVelocities = psmJointActiveCmdVelocities;


				// Retrieve the matrices of the Lagrangian dynamic model
				g = computeDVRKPSMGravityVector(psmJointActivePositions, withCW, withFriction);
				coriolis = computeDVRKPSMCoriolisVector(psmJointActivePositions, psmJointActiveVelocities, withCW, withFriction);
				M = computeDVRKPSMMassMatrix(psmJointActivePositions, withCW, withFriction);

				// Send commands to the simulated robot in CoppeliaSim
				for (int i = 0; i < PSM_ACTIVE_JOINTS; i++) {
					std::string qdSig = std::string("qLdot_") + std::to_string(i + 1);
					simxSetFloatSignal(clientID, qdSig.c_str(), psmJointActiveCmdVelocities(i), simx_opmode_oneshot);
					simxSetFloatSignal(clientID, (qdSig + std::string("_msr")).c_str(), psmJointActiveMsrVelocities(i), simx_opmode_oneshot);
				}

			}
			else if (ctrlType == CTRL_TYPE::DYNAMIC) {

				// active joints J1, J2, J3 (RRP) -> dynamic control
				// active joints J4, J5, J6+J7 (PPP, wrist) -> no control

				// Retrieve the matrices of the Lagrangian dynamic model (execution time evaluation: compute previous and subsequent time instants)
				// Measure time before retrieving the dynamic model ...
				auto tic_ = std::chrono::high_resolution_clock::now();

				g = computeDVRKPSMGravityVector(psmJointActivePositions, withCW, withFriction);
				coriolis = computeDVRKPSMCoriolisVector(psmJointActivePositions, psmJointActiveVelocities, withCW, withFriction);
				M = computeDVRKPSMMassMatrix(psmJointActivePositions, withCW, withFriction);

				// ... and after having retrieved the dynamic model
				auto toc_ = std::chrono::high_resolution_clock::now();

				// Evaluate the elapsed time
				auto tictoc_ = std::chrono::duration_cast<std::chrono::microseconds>(toc_ - tic_).count() * 1e-6;
				//std::cout << "Dynamic model retrieval exec. time = " << tictoc_ << " [s]" << std::endl;

				// Dynamic control law
				psmActiveJointsVectorf qdot = psmJointActiveVelocities;
				tauCmd.setZero();

				if (choose == TASK_TYPE::REGULATION_TASK) {

					Eigen::Vector3f a = Kpd * (-J.topLeftCorner(3, 3) * qdot.topRows(3)) + Kpp * err_p + Kpi * errI_p;
					tauCmd.head(3) = M.topLeftCorner(3, 3) * Jlpinv * a + coriolis.head(3) + g.head(3) - M.topLeftCorner(3, 3) * Jlpinv * Jldot * qdot.head(3);

				}
				else if (choose < TASK_TYPE::EXCITING_JOINT_TRAJECTORY_TASK) {
					Eigen::Vector3f a = pdotdot_des + Kpd * (pdot_des - J.topLeftCorner(3, 3) * qdot.topRows(3)) + Kpp * err_p + Kpi * errI_p;
					tauCmd.head(3) = M.topLeftCorner(3, 3) * Jlpinv * a + coriolis.head(3) + g.head(3) - M.topLeftCorner(3, 3) * Jlpinv * Jldot * qdot.head(3);
				}
				else if (choose == TASK_TYPE::EXCITING_JOINT_TRAJECTORY_TASK) {
					Eigen::Vector3f a = qdd_des + Kqp * errP_q + Kqd * errD_q;
					tauCmd.head(3) = M.topLeftCorner(3, 3) * a + coriolis.head(3) + g.head(3);
				}


			}

		}
		else {
			psmJointActivePositions = offData[excitCurIdx].head(PSM_ACTIVE_JOINTS);
			psmJointActiveVelocities = offData[excitCurIdx].segment(PSM_ACTIVE_JOINTS, 2 * PSM_ACTIVE_JOINTS - 1);
			psmJointActiveCmdVelocities = psmJointActiveVelocities;
			tauCmd = offData[excitCurIdx].tail(PSM_ACTIVE_JOINTS);
		}

		// Send model-based joint velocities to the simulated robot in CoppeliaSim
		for (int i = 0; i < PSM_ACTIVE_JOINTS; i++) {
			std::string qdSig = std::string("qLdot_") + std::to_string(i + 1);
			simxSetFloatSignal(clientID, qdSig.c_str(), psmJointActiveCmdVelocities(i), simx_opmode_oneshot);
			simxSetFloatSignal(clientID, (qdSig + std::string("_msr")).c_str(), psmJointActiveMsrVelocities(i), simx_opmode_oneshot);
		}

		// Send commands to the simulated robot in CoppeliaSim
		for (int i = 0; i < PSM_ACTIVE_JOINTS; i++) {
			std::string tauSig = std::string("tauL_") + std::to_string(i + 1);
			simxSetFloatSignal(clientID, tauSig.c_str(), tauCmd(i), simx_opmode_oneshot);
		}

		// Send joint measurement signals to CoppeliaSim
		for (int i = 0; i < PSM_ACTIVE_JOINTS; i++) {
			std::string tauSig = std::string("tauLMsr_") + std::to_string(i + 1);
			simxSetFloatSignal(clientID, tauSig.c_str(), tauMsr(i), simx_opmode_oneshot);
		}

		// Compute measured joint accelerations
		psmJointActiveAccelerations = (psmJointActiveVelocities - psmJointActivePrevVelocities) / Tsim;

		// Filter acceleration
		psmJointActiveAccelerations = filter_accelleration(psmJointActiveAccelerations, qddState);
		//psmJointActiveAccelerations = alpha * psmJointActiveAccelerations + (1.0 - alpha) * qddot_f_prev;
		qddot_f_prev = psmJointActiveAccelerations;
		//psmJointActiveAccelerations = filter_accelleration(psmJointActiveAccelerations, qddState);

		// Compute model torque
		tauModel = computePredictedTorque(psmJointActivePositions, psmJointActiveVelocities, psmJointActiveAccelerations, withCW, withFriction);

		// Send current EE position signals
		for (int i = 0; i < 3; i++) {
			simxSetFloatSignal(clientID, posSig[i].c_str(), pee(i), simx_opmode_oneshot);
		}

		// Send current desired joint position signals
		for (int i = 0; i < 3; i++) {
			simxSetFloatSignal(clientID, qdesSig[i].c_str(), qdes(i), simx_opmode_oneshot);
		}


		// Send error signals
		for (int i = 0; i < 3; i++) {
			simxSetFloatSignal(clientID, errPosSig[i].c_str(), err_p(i), simx_opmode_oneshot);
		}

		for (int i = 0; i < 3; i++) {
			simxSetFloatSignal(clientID, errOriSig[i].c_str(), (u_ori)(i), simx_opmode_oneshot);
		}

		// Send joint velocity and accelerations
		for (int i = 0; i < PSM_ACTIVE_JOINTS; i++) {
			std::string qddSig = std::string("qLdd_") + std::to_string(i + 1);
			simxSetFloatSignal(clientID, qddSig.c_str(), psmJointActiveAccelerations(i), simx_opmode_oneshot); // accelerations 
		}

		// Send model torques on CoppeliaSim for plotting
		for (int i = 0; i < 3; i++) {
			simxSetFloatSignal(clientID, jointTorqueSigNames[i].c_str(), -tauModel(i), simx_opmode_oneshot); // tau from the model with kinematic control 
		}

		// Trigger the next simulation step on CoppeliaSim (if synchronous mode is enabled)
		if (sync) {
			simxInt pingTime;
			simxSynchronousTrigger(clientID);
			simxGetPingTime(clientID, &pingTime);
		}

		// Update logs
		// Store data of the current iteration
		qmsrSS << t << ", ";
		for (int i = 0; i < PSM_ACTIVE_JOINTS; i++) {
			qmsrSS << psmJointActivePositions(i) << ((i != PSM_ACTIVE_JOINTS - 1) ? "," : ";\n");
		}

		qdmsrSS << t << ", ";
		for (int i = 0; i < PSM_ACTIVE_JOINTS; i++) {
			qdmsrSS << psmJointActiveMsrVelocities(i) << ((i != PSM_ACTIVE_JOINTS - 1) ? "," : ";\n");
		}

		qddmsrSS << t << ", ";
		for (int i = 0; i < PSM_ACTIVE_JOINTS; i++) {
			qddmsrSS << psmJointActiveAccelerations(i) << ((i != PSM_ACTIVE_JOINTS - 1) ? "," : ";\n");
		}

		qdcmdSS << t << ", ";
		for (int i = 0; i < PSM_ACTIVE_JOINTS; i++) {
			qdcmdSS << psmJointActiveCmdVelocities(i) << ((i != PSM_ACTIVE_JOINTS - 1) ? "," : ";\n");
		}

		tauModelSS << t << ", ";
		for (int i = 0; i < PSM_ACTIVE_JOINTS; i++) {
			tauModelSS << tauModel(i) << ((i != PSM_ACTIVE_JOINTS - 1) ? "," : ";\n");
		}

		tauMsrSS << t << ", ";
		for (int i = 0; i < PSM_ACTIVE_JOINTS; i++) {
			tauMsrSS << tauMsr(i) << ((i != PSM_ACTIVE_JOINTS - 1) ? "," : ";\n");
		}

		tauCmdSS << t << ", ";
		for (int i = 0; i < PSM_ACTIVE_JOINTS; i++) {
			tauCmdSS << tauCmd(i) << ((i != PSM_ACTIVE_JOINTS - 1) ? "," : ";\n");
		}

		gSS << t << ", ";
		for (int i = 0; i < PSM_ACTIVE_JOINTS; i++) {
			gSS << g(i) << ((i != PSM_ACTIVE_JOINTS - 1) ? "," : ";\n");
		}


		posRefSS << t << ", " << pdes(0) << ", " << pdes(1) << ", " << pdes(2) << "; " << std::endl;
		posSS << t << ", " << pee(0) << ", " << pee(1) << ", " << pee(2) << "; " << std::endl;
		oriSS << t << ", " << Ree(0, 0) << ", " << Ree(0, 1) << ", " << Ree(0, 2) << ", "
			<< Ree(1, 0) << ", " << Ree(1, 1) << ", " << Ree(1, 2) << ", "
			<< Ree(2, 0) << ", " << Ree(2, 1) << ", " << Ree(2, 2) << "; " << std::endl;

		oriRefSS << t << ", " << Rdes(0, 0) << ", " << Rdes(0, 1) << ", " << Rdes(0, 2) << ", "
				 << Rdes(1, 0) << ", " << Rdes(1, 1) << ", " << Rdes(1, 2) << ", "
				 << Rdes(2, 0) << ", " << Rdes(2, 1) << ", " << Rdes(2, 2) << "; " << std::endl;

		// Update prev variables
		psmJointActivePrevPositions = psmJointActivePositions;
		psmJointActivePrevVelocities = psmJointActiveVelocities;

		// Update time
		t += Tsim;
		excitTrajTime += 1.0/excitTrajRate;
		excitCurIdx++;

		// Check user requests from keyboard
		if (_kbhit()) {

			if (_getch() == ESCAPE)
				i = iterations; // Exit from the loop

		}


	}

	/// Save logs

	// Get the current program directory
	if (saveLogs) {

		int sessionNum = 0;
		bool folderCreated = false;
		std::string folderName = std::string("LogSession_");
		std::string candidateName, currentLogFolderName;
		std::ofstream outFile;


		// Iterate to create a folder with a not pre-existing name (e.g., from previous logs)
		while (folderCreated == false) {

			candidateName = std::string(curdir) + "\\" + folderName + std::to_string(sessionNum);
			sessionNum++;


			// Try to create the candidate folder
			folderCreated = CreateDirectory(candidateName.c_str(), NULL);

			if (folderCreated) {
				std::cout << "Log Folder #" << sessionNum << " created. " << std::endl;
				currentLogFolderName = candidateName;
			}
		}


		outFile.open((currentLogFolderName + "\\qmsr.txt").c_str(), std::ofstream::out);
		outFile << qmsrSS.str();
		outFile.close();

		outFile.open((currentLogFolderName + "\\qdmsr.txt").c_str(), std::ofstream::out);
		outFile << qdmsrSS.str();
		outFile.close();

		outFile.open((currentLogFolderName + "\\qdcmd.txt").c_str(), std::ofstream::out);
		outFile << qdcmdSS.str();
		outFile.close();

		outFile.open((currentLogFolderName + "\\qddmsr.txt").c_str(), std::ofstream::out);
		outFile << qddmsrSS.str();
		outFile.close();

		outFile.open((currentLogFolderName + "\\pdes.txt").c_str(), std::ofstream::out);
		outFile << posRefSS.str();
		outFile.close();

		outFile.open((currentLogFolderName + "\\Rdes.txt").c_str(), std::ofstream::out);
		outFile << oriRefSS.str();
		outFile.close();

		outFile.open((currentLogFolderName + "\\pee.txt").c_str(), std::ofstream::out);
		outFile << posSS.str();
		outFile.close();

		outFile.open((currentLogFolderName + "\\Ree.txt").c_str(), std::ofstream::out);
		outFile << oriSS.str();
		outFile.close();

		outFile.open((currentLogFolderName + "\\tauModel.txt").c_str(), std::ofstream::out);
		outFile << tauModelSS.str();
		outFile.close();

		outFile.open((currentLogFolderName + "\\tauMsr.txt").c_str(), std::ofstream::out);
		outFile << tauMsrSS.str();
		outFile.close();

		outFile.open((currentLogFolderName + "\\tauCmd.txt").c_str(), std::ofstream::out);
		outFile << tauCmdSS.str();
		outFile.close();

		outFile.open((currentLogFolderName + "\\g.txt").c_str(), std::ofstream::out);
		outFile << gSS.str();
		outFile.close();

		// Ask for notes to save
		std::cout << "Type a note on this log session: (press ENTER to confirm)" << std::endl;
		std::stringstream note;
		std::string input("");

		std::cin.ignore();
		std::getline(std::cin, input);

		std::cout << "note inserted: " << input << std::endl;
		note = std::stringstream(input);
		outFile.open((currentLogFolderName + "\\LogSessionInfo.txt").c_str(), std::ofstream::out);
		outFile << input;
		outFile.close();

	}


	
	// ------------------------------------------------ ^^^ Do stuff ^^^ ----------------------------------------------------------------- //


	// Stop the CoppeliaSim simulation remotely
	simxStopSimulation(clientID, simx_opmode_blocking);

	// Close remote API communications with CoppeliaSim
	simxFinish(clientID);

	return 0;

}


/**
* @brief Get function
* Get the initial fixed transformation Tw0 expressing the pose of the base frame of the arm specified by the input string, wrt the CoppeliaSim world frame
* @param clientID: the remote CoppeliaSim ID
* @param ref0: string with the name of the dummy object representing the base frame of the considered arm
* @return the 4x4 homogeneous transformation with the requested pose
*/
Eigen::Matrix4f getTw0(const int& clientID, const std::string& ref0) {

	Eigen::Matrix4f T;
	Eigen::Vector3f angles_;
	Eigen::Vector3f position_;
	int handle;

	simxGetObjectHandle(clientID, ref0.c_str(), &handle, simx_opmode_blocking);
	simxGetObjectOrientation(clientID, handle, -1, &angles_[0], simx_opmode_blocking);
	simxGetObjectPosition(clientID, handle, -1, &position_[0], simx_opmode_blocking);

	Eigen::Quaternionf R;
	T.setIdentity();

	R = Eigen::AngleAxisf(angles_[0], Eigen::Vector3f::UnitX())
		* Eigen::AngleAxisf(angles_[1], Eigen::Vector3f::UnitY())
		* Eigen::AngleAxisf(angles_[2], Eigen::Vector3f::UnitZ());

	T.block(0, 0, 3, 3) << R.toRotationMatrix();
	T.block(0, 3, 4, 1) << position_[0], position_[1], position_[2], 1;

	return T;
}
