# dynamic-dvrk-coppeliasim-simulator

## Requirements

- [CoppeliaSim](https://www.coppeliarobotics.com/) simulation software^*^
- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)^*^
- (optional) [Visual Studio](https://visualstudio.microsoft.com/)^*^


## Build

- Set the environment variables "EIGEN" and "VREPx64_coppelia" with the corresponding paths in your system (e.g., EIGEN="C:\\Users\\[username]\\Software\\eigen-3.3.9\\eigen-3.3.9" and VREPx64_coppelia="C:\\Program Files\\CoppeliaRobotics\\CoppeliaSimEdu")

- If you use Visual Studio, open the file _dvrkDynModelLib.sln_ in the folder _dvrkDynModelLib_ and build the solution with the following configuration: x64-Release 

- Otherwise, build the library with the dynamic model of the PSM and the test executable with CMake through the provided CMakeLists

^*^ Please note that, at the current stage, the simulator has been validated on systems with the following configurations: CoppeliaSim 4.2.0, Eigen 3.3.9, Visual Studio 2019 (Platform Toolset v142). 
It is possible that build errors may arise with newer versions of the provided list of requirements. For instance, versions of CoppeliaSim >= 4.2.0 may have a different path for storing the remote API files, required to build the simulator test program. 


## Usage and Testing

- Open one of the CoppeliaSim scene located in _scenes_

- Run the executable _dvrkDynModelLib\_test_. The simulation in CoppeliaSim will start automatically

- Select the trajectory to be executed by the arm on the right. Please note that if you want to execute the offline joint trajectories in (option 4), you must paste the .csv file located in _data_ to the path of your executable. In addition, set the simulation time step in CoppeliaSim to dt=5.0ms to comply with the sample time of the data logged in the .csv file.
 


