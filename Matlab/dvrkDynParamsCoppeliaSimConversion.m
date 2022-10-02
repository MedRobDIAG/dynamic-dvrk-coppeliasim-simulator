% dvrkDynParamsCoppeliaSimConversion.m - main file for the conversion of
% the dynamic parameters of the dVRK PSM in the CoppeliaSim dVRK dynamic 
% simulator, retrieved from the list of parameters identified in "Y. Wang, 
% R. Gondokaryono, A. Munawar, and G. S. Fischer, “A convex 
% optimization-based dynamic model identification package for the da vinci 
% research kit,” IEEE Robotics and Automation Letters, 2019."
%
% This script comes with the CoppeliaSim dVRK dynamic simulator presented
% in:
% "M. Ferro, A. Mirante, F. Ficuciello, M. Vendittelli, "A CoppeliaSim
% dynamic simulator of the daVinci Research Kit", submitted to IEEE
% Robotics and Automation Letters, 2022.
%
% The script requires the MATLAB files for remote communication with
% CoppeliaSim through remote legacy APIs (remApi.m, remoteApiProto.m,
% remoteApi.dll/so/dylib)
%
% Authors: Marco Ferro, Alessandro Mirante, Fanny Ficuciello, Marilena
% Vendittelli
% Contact: marco.ferro@irisa.fr

% Oct 2022; Last revision: 02-Sep-2022
close all;
clear;
clc;

F0frame_name = 'DH0_ref';
fulljoint_names = {'J1_PSM1','J2_PSM1','','J22_PSM1','J23_PSM1','J24_PSM1','J25_PSM1','J3_PSM1','J31_PSM1','J4_PSM1','J5_PSM1','J6_PSM1','J7_PSM1',''};
link_names = {'L1_respondable_PSM1','Lp11_respondable_PSM1','','Lp2_respondable_PSM1','Lp1_respondable_PSM1','L2_respondable_PSM1','Lp21_respondable_PSM1','L3_respondable_PSM1','CW_respondable'};
activejoint_idxs = [1 2 8 10 11 12 13];
activejoint_names = fulljoint_names(activejoint_idxs);
activejoint_handles = zeros(1,numel(activejoint_idxs));
link_handles = zeros(1,numel(link_names));


%%                     Retrived Dynamic Parameters
% Wang barycentric paramters
% Il=[1.318395, 4.47e-5, 0.0, 1.3192812, 0.0, 0.00093;
%     0.0769277, -0.00156765, 0.0069185, 0.0316068, -0.0039539, 0.0603379;
%     0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
%     0.0021128, -0.0103894, -0.0090333, 90.1865526, -0.001028, 90.1868409;
%     0.0009005, 0.0042873, -0.004286, 37.088943, 0.0004319, 37.0889432;
%     0.0003618, 2.07e-5, 3.59e-5, 0.0001666, 0.0001348, 0.0003194;
%     0.0100262, -0.0031696, 0.0004199, 0.0013188, 0.0014119, 0.0109082;
%     0.0037648, -0.0001809, -0.0007013, 0.0038083, -0.0007102, 0.0003704;
%     0.0036083, 0.000204, 0.0009268, 0.0037699, -0.0006548, 0.0004368];

% Wang standard parameters
Il=[1.3183903, 0.0, 0.0, 1.3183907, 0.0, 3.73e-5;
    0.0235898, -0.0015899, 0.0090817, 0.0265699, 0.0040484, 0.0044209;
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    3.61e-5, -2.87e-5, -2.4e-5, 90.0948874, 4.0e-7, 90.0948873;
    3.69e-5, 2e-7, -1e-7, 37.04596, 0.0, 37.04596;
    2.56e-5, 2.07e-5, 3.59e-5, 8.04e-5, -1.2e-5, 6.94e-5;
    2.61e-5, 2.06e-5, -3.58e-5, 8.04e-5, 1.19e-5, 6.96e-5;
    2.54e-5, -3.58e-5, 2.07e-5, 6.9e-5, 1.2e-5, 8.02e-5;
    2.58e-5, 3.6e-5, 2.11e-5, 6.97e-5, -1.23e-5, 8.07e-5];

% link masses
m =[6.0667823
    10.0001371
    0
    2.956747
    0.4796773
    0.1000014
    0.5000069
    0.3630268
    0.3499395]; 

% CoMs
CoM_l =[ 0.0121153,-0.0006075,0.0;
      -0.0195141,-0.0721861,-0.0110856;
       0.0, 0.0, 0.0;
      0.1752131, 0.0199989, 0.0173905;
      0.2978392,-0.030008,0.0299994;
      -7.0e-7,0.0500018,-0.0293552;
       -0.0455733,-0.1400001, 0.02;
       -0.0199883, -0.0199936, -0.0995037;
       0.0260167, -0.018457, -0.0994826]';
   
%% CoppeliaSim connection
disp('Program started');
sim=remApi('remoteApi'); % using the prototype file (remoteApiProto.m)
sim.simxFinish(-1); % just in case, close all opened connections
clientID=sim.simxStart('127.0.0.1',19997,true,true,5000,5);

if (clientID>-1)
    disp('Connected to remote API server');
else
    disp('Could not connect to the remote API server');
end

if (clientID>-1)

    % Read the handles of the robot joints from the scene
    joint_position = zeros(numel(activejoint_names),1);
    for i = 1 : length(activejoint_names)
        [ret, activejoint_handles(i)] = sim.simxGetObjectHandle(clientID,activejoint_names{i},sim.simx_opmode_blocking);
        [ret, joint_position(i)] = sim.simxGetJointPosition(clientID,activejoint_handles(i),sim.simx_opmode_blocking);
    end
    
    % Read the handles of the robot links from the scene
    for i = 1 : length(link_names)
        [ret, link_handles(i)] = sim.simxGetObjectHandle(clientID,link_names{i},sim.simx_opmode_blocking);
    end
    
    Tw0 = checkRefDH0(clientID,sim);   % Receive DH0_ref position: position of the DV's zero frame 
end
   
DH=completeDH(joint_position);     % DH table creation considers the complete One (i.e. Include also closed branch)
chain = chain_position_complete(DH,Tw0); % Trasformation chain

%% Transform the values to be used in CoppeliaSim
w_CoM = zeros(size(CoM_l));
J = zeros(3,3,length(joint_position));
a =['x';'y';'z'];
for i=1:9
    aux = chain(:,:,i)*[CoM_l(:,i);1];
    w_CoM(:,i) = aux(1:3);
    msg = strcat('dVRK-PSM_link ',num2str(i+1));
    disp(msg)
    msg = strcat('mass =  ',num2str(m(i)),' [Kg]');
    disp(msg)
    disp('Inertia matrix [m^2]= ')
    R = chain(1:3,1:3,i);
    T = [R,CoM_l(:,i) ;0 0 0 1];
    % Rotate the inertia tensor / apply Generalized Steiner in case of
    % barycentric coordinates
    J(:,:,i) = R * vec2sm(Il(i,:)) * R';%GeneralizedSteiner(T,vec2sm(Il(i,:)),m(i));
    msg = num2str(J(:,:,i)/m(i));
    T = [chain(1:3,1:3,i),w_CoM(:,i);0 0 0 1];
    disp(msg)
    disp('CoM in World Coord. = ')
    msg = strcat(a,' = ',num2str(w_CoM(:,i)),' [m]');
    disp(msg)
    disp('-------------------------------------------')
end

% Set the inertia matrices and COMs
% sim.simxCallScriptFunction(clientID, 'ApiRemoteBridge', sim.sim_scripttype_childscript, 'setLinkMassAndInertias', link_handles, [boxSizes(i,:),boxCenters(i,:)], '', [], sim.simx_opmode_blocking);

% Now close the connection to CoppeliaSim
sim.simxFinish(clientID);
sim.delete(); % call the destructor!
disp('Program ended');

