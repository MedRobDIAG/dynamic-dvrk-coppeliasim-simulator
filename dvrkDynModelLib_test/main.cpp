// dvrk dyn model library header file
#include "dvrkDynamics.hpp"
#include "dvrkKinematics.hpp"
#include "Timer.hpp"
#include "utils.hpp"

#define ESCAPE 27

// Remote CoppeliaSim API functions
extern "C" {
#include "extApi.h"
}

// Integer flags specifying the type of control used for the robot motion
// Kinematic: computes and sends joint velocity inputs
// Dynamic: computes and sends joint torque inputs
enum CTRL_TYPE {UNDEFINED = -1, KINEMATIC, DYNAMIC};
enum TASK_TYPE {REGULATION_TASK = 1, LINEAR_TRAJECTORY_TASK, SPIRAL_TRAJECTORY_TASK, JOINT_REGULATION_TASK, JOINT_TRAJECTORY_TASK};

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

	int choose;
	bool saveLogs = false; //  Set this flag to true if you want to log data from the simulation
	float ctrlType = CTRL_TYPE::UNDEFINED;
	bool sync = true; // synchronous flag for CoppeliaSim simulation
	bool runPlayback = false;
	bool withCW = true;
	bool withFriction = false;
	int logNum = 0;
	char curdir[256];
	int j3MirrorHandle = -1;
	simxFloat tau3Mirror = 0.0;
	std::vector < Eigen::VectorXf > offData;
	std::string offlineDatasetPath;
	//std::string jointNames[PSM_FULL_JOINTS] = { "J1_PSM1","J2_PSM1","","Jp11_PSM1","Jp1_PSM1","Jp2_PSM1","Jp21_PSM1","J3_PSM1","","J1_TOOL1","J2_TOOL1","J3_sx_TOOL1","J3_dx_TOOL1","" }; // Joint names in CoppeliaSim simulation
	std::string jointNames[PSM_FULL_JOINTS] = { "J1_PSM1","J2_PSM1","","J23_PSM1","J22_PSM1","J24_PSM1","J25_PSM1","J3_PSM1","","J4_PSM1","J5_PSM1","J6_PSM1","J7_PSM1","" }; // Joint names in CoppeliaSim simulation
	std::string jointTorqueSigNames[3] = { "Torque_J1","Torque_J2","Force_J3" };
	std::string jointMirror3Name = "J3_mirror_PSM1";
	std::string posSig[3] = { "x_ee","y_ee","z_ee" };
	std::string errPosSig[3] = { "error_x","error_y","error_z" };
	std::string errJointSig[3] = { "error_q1","error_q2","error_q3" };
	std::string errOriSig[3] = { "error_alfa","error_beta","error_gamma" };
	psmFullJointsVectori jointHandles;
	psmFullJointsVectorf psmJointPositions, psmJointMsrTorques;
	psmActiveJointsVectorf psmJointActivePositions, psmJointActivePrevPositions, psmJointActiveVelocities, psmJointActivePrevVelocities, psmJointActiveCmdVelocities, psmJointActiveMsrVelocities, psmJointActiveAccelerations;
	psmActiveJointsVectorf tauModel, tauLagModel, tauMsr, tauCmd;
	psmActiveJointsVectorf qdes, qdotdes, qdotdotdes, errq, errqD;
	psmActiveJointsVectorf g;
	Eigen::VectorXf qdState, qddState;
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
	std::stringstream qErrSS("");
	std::stringstream qDesSS("");
	std::stringstream gSS("");

	// Reference trajectory variables
	Eigen::VectorXf refTrajectoryVisual;
	Eigen::VectorXf draw_RefTrajectory_x, draw_RefTrajectory_y, draw_RefTrajectory_z;
	Eigen::VectorXf refTrajectory_x, refTrajectory_y, refTrajectory_z, refTrajectory_q1, refTrajectory_q2, refTrajectory_q3;
	Eigen::VectorXf refTrajectoryDiff_x, refTrajectoryDiff_y, refTrajectoryDiff_z, refTrajectoryDiff_q1, refTrajectoryDiff_q2, refTrajectoryDiff_q3;
	Eigen::VectorXf refTrajectoryDiffDiff_x, refTrajectoryDiffDiff_y, refTrajectoryDiffDiff_z, refTrajectoryDiffDiff_q1, refTrajectoryDiffDiff_q2, refTrajectoryDiffDiff_q3;;
	Eigen::Vector3f startingPoint;
	Eigen::Vector3f refOrientation;
	Eigen::Vector3f pdes;
	Eigen::Vector3f pdot_des;
	Eigen::Vector3f pdotdot_des;
	Eigen::Matrix3f Rdes;
	Eigen::Quaternionf quatdes;
	double t = 0.0;
	double trajDuration = 10.0; // [s]
	double Tsim = 0.001; // [s]
	int totSteps;
	int startHandle; // Handle of dummy in simulated scene specifying the starting point of the trajectory

	// Robot and control variables
	Eigen::Matrix4f Tee;
	Eigen::Matrix3f Ree;
	Eigen::Vector3f pee, rpy_ee;
	Eigen::Vector3f err_p, errI_p, errD_p;
	Eigen::Quaternionf quat_ee, err_quat;
	//Eigen::Matrix<float, 6, PSM_ACTIVE_JOINTS> J;
	Eigen::Matrix<float, 6, PSM_ACTIVE_JOINTS-1> J;
	Eigen::Matrix<float, 3, 3> Jl, Jo, Jlprev,Jldot;
	Eigen::Matrix<float, PSM_ACTIVE_JOINTS, 6> JT, Jpinv;
	Eigen::Matrix<float, 3, 3> JlT, Jlpinv;
	Eigen::Vector6f u_cart;
	Eigen::Vector3f u_pos, u_ori;
	Eigen::Matrix3f Kpp, Kpd, Kpi, Kop, Kqp, Kqd;
	Eigen::Matrix<float, PSM_ACTIVE_JOINTS, PSM_ACTIVE_JOINTS> Kqp7, Kqd7;
	psmActiveJointsVectorf u_ctrl;
	psmActiveJointsVectorf qddmodel;
	psmActiveJointsVectorf coriolis;
	psmMassMatrixf M;
	Eigen::Vector4f quat_eev;
	Kpp.setZero();
	Kpd.setZero();
	Kpi.setZero();
	Kop.setZero();
	Kqp.setZero();
	Kqd.setZero();
	Kqp7.setZero();
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


	// Get the PSM joint handles
	psmJointActivePositions.setZero();
	psmJointActivePrevPositions.setZero();
	psmJointActiveVelocities.setZero();
	psmJointActivePrevVelocities.setZero();
	psmJointActiveAccelerations.setZero();
	psmJointMsrTorques.setZero();
	tauModel.setZero();
	tauLagModel.setZero();
	tauCmd.setZero();
	jointHandles.setZero();
	psmJointPositions.setZero();
	qddmodel.setZero();

	// Get robot base frame CoppeliaSim handle
	simxInt psm1BaseHandle;
	simxInt dummyHandle;
	simxGetObjectHandle(clientID, "DH0_ref", &psm1BaseHandle, simx_opmode_blocking);
	simxGetObjectHandle(clientID, "eeDummy", &dummyHandle, simx_opmode_blocking);
	simxGetObjectHandle(clientID, "J3mirror_PSM1", &j3MirrorHandle, simx_opmode_blocking);
	std::cout << "dummy handles taken. " << std::endl;
	// Get initial joint configuration
	for (int i = 0; i < PSM_FULL_JOINTS; i++) {
		simxGetObjectHandle(clientID, jointNames[i].c_str(), &jointHandles[i], simx_opmode_blocking);
		simxGetJointPosition(clientID, jointHandles[i], &psmJointPositions(i), simx_opmode_blocking);
		simxGetJointForce(clientID, jointHandles[i], &psmJointMsrTorques(i), simx_opmode_blocking);
	}
	std::cout << "jointHandles = " << jointHandles.transpose() << std::endl;
	// Extract active joint position
	psmJointActivePositions = (activeJointsIdxs.unaryExpr(psmJointPositions)).matrix();
	tauMsr = (activeJointsIdxs.unaryExpr(psmJointMsrTorques)).matrix();
	psmJointActivePrevPositions = psmJointActivePositions;

	// Get the reference trajectory starting point coordinates
	simxGetObjectHandle(clientID, "Start", &startHandle, simx_opmode_blocking);
	simxGetObjectPosition(clientID, startHandle, -1, &startingPoint[0], simx_opmode_blocking);

	// Let the user choose which trajectory should be executed by the robot autonomously
	// r: rectilinear
	// s: spiral
	std::cout << "Select task type: " <<
		"\n\t1. regulation task" <<
		"\n\t2. tracking task -- linear trajectory" <<
		"\n\t3. tracking task -- spiral trajectory" << 
		"\n\t4. joint regulation task " << 
		"\n\t5. joint trajectory tracking task " << std::endl;
	std::cin >> choose;

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
	totSteps = int(trajDuration / Tsim);

	int num_microstep = 50; // Set to 1 to make this step uneffective
	double Ts = Tsim / num_microstep;
	int bufferVisualSize = trajDuration / Ts;
	int bufferSize = int(trajDuration / Tsim);
	int xyzBufferVisualSize = 3 * bufferVisualSize;

	refTrajectoryVisual.setZero(xyzBufferVisualSize);
	draw_RefTrajectory_x.setZero(bufferVisualSize);
	draw_RefTrajectory_y.setZero(bufferVisualSize);
	draw_RefTrajectory_z.setZero(bufferVisualSize);

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
	refTrajectoryDiffDiff_q1.setZero(bufferSize);
	refTrajectoryDiffDiff_q2.setZero(bufferSize);
	refTrajectoryDiffDiff_q3.setZero(bufferSize);

	if (choose == TASK_TYPE::REGULATION_TASK) {
		simxInt eeRef;
		simxGetObjectHandle(clientID, "eeRef_regulation", &eeRef, simx_opmode_blocking);
		simxGetObjectPosition(clientID, eeRef, -1, &pdes[0], simx_opmode_blocking);
		std::cout << "pdes = " << pdes.transpose() << std::endl;

	}
	else if (choose == TASK_TYPE::LINEAR_TRAJECTORY_TASK) { // Linear trajectory
	
		std::cout << "\n\nA straigth line trajectory has been selected, continue on Coppelia Sim\n";
		Eigen::Vector3d angular_coef(-0.003, -0.003, 0.003);
		//Eigen::Vector3d angular_coef(-0.003, 0.005, 0.005);

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
		//Eigen::Vector3d Frequency(0.75, 1.5, 0.3);
		Eigen::Vector3d Frequency(0.75, 0.75, 0.1);
		//Eigen::Vector3d Frequency(2, 2, 1);
		//Eigen::Vector3d Frequency(0.2, 0.2, 0.2);
		
		//Eigen::Vector3d Amplitude(0.025, 0.025, 0.02);
		//Eigen::Vector3d Phase(0, 0, 0);
		//Eigen::Vector3d Frequency(3.0, 2.0, 5.0);


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
		for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectoryDiffDiff_x(ix) = -pow(2.0 * M_PI * Frequency(0),2) * Amplitude(0) * cos(2.0 * M_PI * Frequency(0) * i + Phase(0)); ix++; }
		for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectoryDiffDiff_y(iy) = -pow(2.0 * M_PI * Frequency(1),2) * Amplitude(1) * sin(2.0 * M_PI * Frequency(1) * i + Phase(1)); iy++; }
		for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectoryDiffDiff_z(iz) = -pow(2.0 * M_PI * Frequency(2),2) * Amplitude(2) * cos(2.0 * M_PI * Frequency(2) * i + Phase(2)); iz++; }

	}
	else if (choose == TASK_TYPE::JOINT_REGULATION_TASK) { // Joint regulation

		float q1des = 30.0 * M_PI / 180.0;
		float q2des = 20.0 * M_PI / 180.0;
		//qdes << q1des, q2des, 0.0, 0.0, 0.0, 0.0, 0.0;
		qdes << 0.462104, -0.402096, 0.114593;
	}
	else if(choose == TASK_TYPE::JOINT_TRAJECTORY_TASK) { // Joint trajectory tracking task
		
		for (int i = 0; i < PSM_FULL_JOINTS; i++) {
			simxGetJointPosition(clientID, jointHandles[i], &psmJointPositions[i], simx_opmode_blocking);
		}
		psmActiveJointsVectorf qInit = (activeJointsIdxs.unaryExpr(psmJointPositions)).matrix();

		Eigen::Vector3d Amplitude(0.1, 0.1, 0.02);
		Eigen::Vector3d Phase(0, 0, 0);
		//Eigen::Vector3d Frequency(2.0, 3.0, 5.0);
		Eigen::Vector3d Frequency(0.75, 0.75, 0.1);
		double ix = 0, iy = 0, iz = 0;
		for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectory_q1(ix) = Amplitude(0) * sin(2.0 * M_PI * Frequency(0) * i + Phase(0)) + qInit(0); ix++; }
		for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectory_q2(iy) = Amplitude(1) * sin(2.0 * M_PI * Frequency(1) * i + Phase(1)) + qInit(1); iy++; }
		for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectory_q3(iz) = Amplitude(2) * sin(2.0 * M_PI * Frequency(2) * i + Phase(2)) + qInit(2); iz++; }

		// Trajectory derivative buffer for control (sample time Tsim)
		ix = 0; iy = 0; iz = 0;
		for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectoryDiff_q1(ix) = 2.0 * M_PI * Frequency(0) * Amplitude(0) * cos(2.0 * M_PI * Frequency(0) * i + Phase(0)); ix++; }
		for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectoryDiff_q2(iy) = 2.0 * M_PI * Frequency(1) * Amplitude(1) * cos(2.0 * M_PI * Frequency(1) * i + Phase(1)); iy++; }
		for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectoryDiff_q3(iz) = 2.0 * M_PI * Frequency(2) * Amplitude(2) * cos(2.0 * M_PI * Frequency(2) * i + Phase(2)); iz++; }

		// Trajectory second derivative buffer for control (sample time Tsim)
		ix = 0; iy = 0; iz = 0;
		for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectoryDiffDiff_q1(ix) = -pow(2.0 * M_PI * Frequency(0), 2) * Amplitude(0) * sin(2.0 * M_PI * Frequency(0) * i + Phase(0)); ix++; }
		for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectoryDiffDiff_q2(iy) = -pow(2.0 * M_PI * Frequency(1), 2) * Amplitude(1) * sin(2.0 * M_PI * Frequency(1) * i + Phase(1)); iy++; }
		for (double i = 0; i < trajDuration - Tsim; i = i + Tsim) { refTrajectoryDiffDiff_q3(iz) = -pow(2.0 * M_PI * Frequency(2), 2) * Amplitude(2) * sin(2.0 * M_PI * Frequency(2) * i + Phase(2)); iz++; }

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
	std::pair < std::vector<Eigen::Matrix4f>, Eigen::Matrix<float,6,PSM_ACTIVE_JOINTS> > psmKin = psmKinematics(psmJointPositions, Tw0);

	// Draw the reference trajectory on CoppeliaSim through remote call of script function
	simxCallScriptFunction(clientID, "ApiRemoteBridge", sim_scripttype_childscript, "drawTrajectory", 1, &bufferVisualSize, xyzBufferVisualSize, &refTrajectoryVisual(0), 0, nullptr, 0, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, simx_opmode_blocking);

	// Forward the type of chosen
	// Initialize streaming for joint position readings
	for (int i = 0; i < PSM_FULL_JOINTS; i++) {
		simxGetJointPosition(clientID, jointHandles[i], &psmJointPositions[i], simx_opmode_streaming);
		simxGetJointForce(clientID, jointHandles[i], &psmJointMsrTorques[i], simx_opmode_streaming);
	}
	psmJointActivePositions = (activeJointsIdxs.unaryExpr(psmJointPositions)).matrix();
	tauMsr = (activeJointsIdxs.unaryExpr(psmJointMsrTorques)).matrix();

	simxGetJointForce(clientID, j3MirrorHandle, &tau3Mirror, simx_opmode_streaming);


	while (ctrlType == CTRL_TYPE::UNDEFINED) {
		simxGetFloatSignal(clientID, "dynamicCtrlFlagSignal", &ctrlType, simx_opmode_blocking);

		// Trigger the next simulation step on CoppeliaSim (if synchronous mode is enabled)
		if (sync) {
			simxSynchronousTrigger(clientID);
		}
	}

	std::cout << "ctrlType = " << ctrlType << std::endl;//*/

	if (choose == TASK_TYPE::REGULATION_TASK) {
		if (ctrlType == CTRL_TYPE::KINEMATIC) {

			Kpp.diagonal() << 5.0, 5.0, 5.0;
			Kop.diagonal() << 5.0, 5.0, 5.0;
			Kpd.diagonal() << 0.0, 0.0, 0.0;
			Kpi.diagonal() << 0.1, 0.1, 0.1;

		}
		else if (ctrlType == CTRL_TYPE::DYNAMIC) {

			//Kpp.diagonal() << 800.0, 800.0, 200.0;
			//Kpd.diagonal() << 5.0, 5.0, 10.0;
			//Kpi.diagonal() << 0.1, 0.1, 0.1;
			Kqd.diagonal() << 5.0, 5.0, 10.0;

			//Kpp.diagonal() << 400.0, 400.0, 400.0;
			//Kpd.diagonal() << 40.0, 40.0, 40.0;
			Kpp.diagonal() << 100.0, 100.0, 100.0;
			Kpd.diagonal() = 2.0 * Kpp.diagonal().cwiseSqrt();
			//Kpi.diagonal() << 10.0, 10.0, 10.0;

		}
	}
	else if (choose == TASK_TYPE::JOINT_REGULATION_TASK) {

		if (ctrlType == CTRL_TYPE::KINEMATIC) {
			Kqp7.diagonal() << 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0;
		}
		else if (ctrlType == CTRL_TYPE::DYNAMIC) {
			Kqp.diagonal() << 40.0, 40.0, 40.0;
			Kqd.diagonal() << 5.0, 5.0, 5.0;// 2.0 * Kpp.diagonal().cwiseSqrt();
		}
	}
	else if (choose == TASK_TYPE::LINEAR_TRAJECTORY_TASK || choose == TASK_TYPE::SPIRAL_TRAJECTORY_TASK){

		if (ctrlType == CTRL_TYPE::KINEMATIC) {
			Kpp.diagonal() << 15.0, 15.0, 15.0;
			Kop.diagonal() << 15.0, 15.0, 15.0;
			Kpd.diagonal() << 0.0, 0.0, 0.0;
			Kpi.diagonal() << 0.0005, 0.0005, 0.0005;
		}
		else if (ctrlType == CTRL_TYPE::DYNAMIC) {
			//Kpp.diagonal() << 800.0, 800.0, 800.0;
			//Kpd.diagonal() << 56.56, 56.56, 20.0;
			//Kpi.diagonal() << 0.1, 0.1, 0.1;
			Kqd.diagonal() << 5.0, 5.0, 10.0;

			Kpp.diagonal() << 250.0, 250.0, 500.0;
			Kpd.diagonal() << 31.0, 31.0, 20.0;// 2.0 * Kpp.diagonal().cwiseSqrt();
			//Kpi.diagonal() << 10.0, 10.0, 10.0;

		}

	}
	else if (choose == TASK_TYPE::JOINT_TRAJECTORY_TASK) {
		if (ctrlType == CTRL_TYPE::KINEMATIC) {
			Kqp7.diagonal() << 50.0, 50.0, 50.0, 20.0, 20.0, 20.0, 20.0;
			Kqd7.diagonal() << 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2;
		}
		else if (ctrlType == CTRL_TYPE::DYNAMIC) {
			Kqp.diagonal() << 250.0, 250.0, 250.0;
			Kqd.diagonal() << 31.0, 31.0, 31.0;
			Kqp7.diagonal() << 50.0, 50.0, 50.0, 20.0, 20.0, 20.0, 20.0;
			Kqd7.diagonal() << 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2;
		}

	}

	int iterations = (runPlayback) ? offData.size() : totSteps;
	
	if(ctrlType == CTRL_TYPE::KINEMATIC){
		offlineDatasetPath = std::string(curdir) + "\\LogSession_" + std::to_string(logNum) + "\\tauMsr.txt";
	}
	else if (ctrlType == CTRL_TYPE::DYNAMIC) {
		offlineDatasetPath = std::string(curdir) + "\\LogSession_" + std::to_string(logNum) + "\\qdcmd.txt";
	}
	offData = loadCSVLog(offlineDatasetPath);
	std::cout << "Offline data loaded. Size = " << offData.size() << std::endl;

	float errq_I = 0.0;

	/////////////////// -------------- Start control loop
	for (int i = 0; i < iterations; i++) {

		//std::cout << "iterations = " << iterations << std::endl;
		//std::cout << "i = " << i << " (" << (float)i / (float)iterations << " %) " << std::endl;
		//std::cout << std::endl;

		// Get the current joint readings
		for (int k = 0; k < PSM_FULL_JOINTS; k++) {
			simxGetJointPosition(clientID, jointHandles[k], &psmJointPositions[k], simx_opmode_buffer);
			simxGetJointForce(clientID, jointHandles[k], &psmJointMsrTorques[k], simx_opmode_buffer);
		}
		simxGetJointForce(clientID, j3MirrorHandle, &tau3Mirror, simx_opmode_buffer);

		// Extract active joint position and torques
		psmJointActivePositions = (activeJointsIdxs.unaryExpr(psmJointPositions)).matrix();
		tauMsr = (activeJointsIdxs.unaryExpr(psmJointMsrTorques)).matrix();
		//tauMsr[2] += tau3Mirror;

		// Compute measured joint velocities
		psmJointActiveMsrVelocities = (psmJointActivePositions - psmJointActivePrevPositions) / Tsim;
		psmJointActiveVelocities = psmJointActiveMsrVelocities;

		psmJointActiveMsrVelocities = filter_velocity(psmJointActiveMsrVelocities, qdState);

		// Compute measured joint accelerations
		//psmJointActiveAccelerations = (psmJointActiveMsrVelocities - psmJointActivePrevVelocities) / Tsim;

		if (!runPlayback) {

			// Get the current dVRK-PSM EE pose from the measured joint configuration
			psmKin.first.clear();
			psmKin = psmKinematics(psmJointPositions, Tw0);
			Tee = psmKin.first[PSM_ACTIVE_JOINTS - 1];
			Ree = Tee.topLeftCorner(3, 3);
			pee = Tee.block<3, 1>(0, 3);
			rpy_ee = rot2rpy(Ree);
			quat_ee = Eigen::Quaternionf(Ree.transpose());

			// Set position of the dummy
			simxSetObjectPosition(clientID, dummyHandle, -1, &pee[0], simx_opmode_oneshot);
			// Set orientation of the dummy
			simxSetObjectQuaternion(clientID, dummyHandle, -1, &quat_ee.coeffs()[0], simx_opmode_oneshot);

			// Get the current Jacobian matrix and compute the pseudo-inverse
			J = psmKin.second.leftCols(6);
			Jl = J.topLeftCorner(3, 3);
			JlT = Jl.transpose();
			Jldot = (Jl - Jlprev) / Tsim;
			Jlprev = Jl;
			Jlpinv = Jl.inverse();

			// UNCOMMENT THIS IF YOU WANT TO USE 6X7 JACOBIAN MATRIX INSTEAD OF 6X6
			/*JT = J.transpose();
			Jpinv = JT * (J * JT).inverse();

			Jl = J.topLeftCorner(3, 3);
			JlT = Jl.transpose();
			Jlpinv = JlT * (Jl * JlT).inverse();//*/
			
			// Compute the configuration error
			if (choose == TASK_TYPE::JOINT_REGULATION_TASK) {
				errq = (qdes - psmJointActivePositions);
			}
			else if (choose == TASK_TYPE::JOINT_TRAJECTORY_TASK) {
				qdes << refTrajectory_q1(i), refTrajectory_q2(i), refTrajectory_q3(i), 0.0, 0.0, 0.0, 0.0;
				qdotdes << refTrajectoryDiff_q1(i), refTrajectoryDiff_q2(i), refTrajectoryDiff_q3(i), 0.0, 0.0, 0.0, 0.0;
				qdotdotdes << refTrajectoryDiffDiff_q1(i), refTrajectoryDiffDiff_q2(i), refTrajectoryDiffDiff_q3(i), 0.0, 0.0, 0.0, 0.0;
				errq = (qdes - psmJointActivePositions);
				errqD = (qdotdes - psmJointActiveMsrVelocities);
			}

			// Compute the position error
			if (choose != TASK_TYPE::REGULATION_TASK) {
				pdes << refTrajectory_x(i), refTrajectory_y(i), refTrajectory_z(i);
				pdot_des << refTrajectoryDiff_x(i), refTrajectoryDiff_y(i), refTrajectoryDiff_z(i);
				pdotdot_des << refTrajectoryDiffDiff_x(i), refTrajectoryDiffDiff_y(i), refTrajectoryDiffDiff_z(i);
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

			/*std::cout << "err_quat = " << err_quat.vec().transpose() << ", " << err_quat.w() << std::endl;
			std::cout << "quat_ee = " << quat_ee.vec().transpose() << ", " << quat_ee.w() << std::endl;
			std::cout << "qdes = " << qdes.vec().transpose() << ", " << qdes.w() << std::endl;
			std::cout << "u_ori = " << u_ori.transpose() << std::endl;*/
			//std::cout << "J = \n" << J << std::endl;
			//std::cout << std::endl;



			if (ctrlType == CTRL_TYPE::KINEMATIC) {
				
				// Kinematic control law
				if (choose == TASK_TYPE::JOINT_REGULATION_TASK) {			
					psmJointActiveCmdVelocities = Kqp7 * errq;
				}
				else if (choose == TASK_TYPE::JOINT_TRAJECTORY_TASK) {
					psmJointActiveCmdVelocities = Kqp7 * errq + qdotdes;

					std::cout << "qdes = " << qdes.transpose() << std::endl;
				} 
				else{
					psmJointActiveCmdVelocities.topRows(6) = J.inverse() * u_cart; // [velocities, rad/s]
				}

				// In case of kinematic control use commanded joint velocities to compute the predicted torque below (*)
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
				
				// In case of dynamic control use measured joint velocities to compute the predicted torque below (*)
				//psmJointActiveVelocities = psmJointActiveMsrVelocities;
				// Retrieve the matrices of the Lagrangian dynamic model
				g = computeDVRKPSMGravityVector(psmJointActivePositions, withCW, withFriction);
				coriolis = computeDVRKPSMCoriolisVector(psmJointActivePositions, psmJointActiveVelocities, withCW, withFriction);
				M = computeDVRKPSMMassMatrix(psmJointActivePositions, withCW, withFriction);

				// Dynamic control law
				psmActiveJointsVectorf qdot = psmJointActiveVelocities;
				tauCmd.setZero();

				if (choose == TASK_TYPE::REGULATION_TASK) {

					Eigen::Vector3f a = Kpd * (- J.topLeftCorner(3, 3) * qdot.topRows(3)) + Kpp * err_p + Kpi * errI_p;
					tauCmd.head(3) = M.topLeftCorner(3, 3) * Jlpinv * a + coriolis.head(3) + g.head(3) - M.topLeftCorner(3, 3) * Jlpinv * Jldot * qdot.head(3);

				}
				else if (choose == TASK_TYPE::JOINT_REGULATION_TASK) {
					tauCmd.head(3) = Kqp * errq.head(3) - Kqd * qdot.head(3) + g.head(3);
				}
				else {
					Eigen::Vector3f a = pdotdot_des + Kpd * (pdot_des - J.topLeftCorner(3,3) * qdot.topRows(3)) + Kpp * err_p + Kpi * errI_p;
					tauCmd.head(3) = M.topLeftCorner(3, 3) * Jlpinv * a + coriolis.head(3) + g.head(3) - M.topLeftCorner(3, 3) * Jlpinv * Jldot * qdot.head(3);
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
			}

		}
		else {

			tauCmd = offData[i];

		}


		// Compute measured joint accelerations
		if (ctrlType == CTRL_TYPE::KINEMATIC) {
			psmJointActiveAccelerations = (psmJointActiveVelocities - psmJointActivePrevVelocities) / Tsim;
		}
		else if (ctrlType == CTRL_TYPE::DYNAMIC) {
			//psmJointActiveAccelerations.head(3) = M.topLeftCorner(3, 3).inverse() * (tauCmd.head(3) - coriolis.head(3) - g.head(3));
			//psmJointActiveVelocities = psmJointActiveMsrVelocities + psmJointActiveAccelerations * Tsim;
			psmJointActiveAccelerations = (psmJointActiveVelocities - psmJointActivePrevVelocities) / Tsim;

		}

		psmJointActiveAccelerations = filter_accelleration(psmJointActiveAccelerations, qddState);

		if (choose != TASK_TYPE::JOINT_TRAJECTORY_TASK) {
			// Evaluate model torques vs measured torques (*)
			tauModel = computePredictedTorque(psmJointActivePositions, psmJointActiveVelocities, psmJointActiveAccelerations, withCW, withFriction);
			tauLagModel = M * psmJointActiveAccelerations + coriolis + g;
		}
		else {
			tauModel = computePredictedTorque(psmJointActivePositions, psmJointActiveVelocities, psmJointActiveAccelerations, withCW, withFriction);
			tauLagModel = M * psmJointActiveAccelerations + coriolis + g;

			//tauModel = computePredictedTorque(qdes, qdotdes, qdotdotdes, withCW, withFriction);
			//tauLagModel = M * qdotdotdes + coriolis + g;
		}

		//std::cout << "tauModel = " << tauModel.transpose() << std::endl;
		//std::cout << "tauLagModel = " << tauLagModel.transpose() << std::endl;
		//std::cout << std::endl;

		// Send current EE position signals
		for (int i = 0; i < 3; i++) {
			simxSetFloatSignal(clientID, posSig[i].c_str(), pee(i), simx_opmode_oneshot);
		}

		// Send error signals
		if (choose == TASK_TYPE::JOINT_REGULATION_TASK || choose == TASK_TYPE::JOINT_TRAJECTORY_TASK) {
			for (int i = 0; i < 3; i++) {
				simxSetFloatSignal(clientID, errJointSig[i].c_str(), errq(i), simx_opmode_oneshot);
			}
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
			std::string tauLag = std::string("tauModLag_") + std::to_string(i + 1);
			simxSetFloatSignal(clientID, jointTorqueSigNames[i].c_str(), -tauModel(i), simx_opmode_oneshot); // tau from the model with kinematic control 
			simxSetFloatSignal(clientID, tauLag.c_str(), -tauLagModel(i), simx_opmode_oneshot);   // tau from dynamic control
		}

		// Trigger the next simulation step on CoppeliaSim (if synchronous mode is enabled)
		if (sync) {
			simxInt pingTime;
			simxSynchronousTrigger(clientID);
			simxGetPingTime(clientID,&pingTime);
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

		qDesSS << t << ", ";
		for (int i = 0; i < PSM_ACTIVE_JOINTS; i++) {
			qDesSS << qdes(i) << ((i != PSM_ACTIVE_JOINTS - 1) ? "," : ";\n");
		}

		qErrSS << t << ", ";
		for (int i = 0; i < PSM_ACTIVE_JOINTS; i++) {
			qErrSS << errq(i) << ((i != PSM_ACTIVE_JOINTS - 1) ? "," : ";\n");
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

		outFile.open((currentLogFolderName + "\\qDes.txt").c_str(), std::ofstream::out);
		outFile << qDesSS.str();
		outFile.close();

		outFile.open((currentLogFolderName + "\\qErr.txt").c_str(), std::ofstream::out);
		outFile << qErrSS.str();
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
