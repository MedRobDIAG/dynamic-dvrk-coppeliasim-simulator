function sysCall_init()
    -- do some initialization here
    
    -- Getting simulation time step
    Ts = sim.getSimulationTimeStep()

    -- flags
    dynamicCtrlFlag = 0
    lowlevelCtrlLoopEnabled = 0

    -- Object handles
	J1_PSM = sim.getObjectHandle("J1_PSM1");
	J2_PSM = sim.getObjectHandle("J2_PSM1");
	J3_PSM = sim.getObjectHandle("J3_PSM1");
	J3m_PSM = sim.getObjectHandle("J3mirror_PSM1");
	JCW_PSM = sim.getObjectHandle("J31_PSM1");
	J1_TOOL1 = sim.getObjectHandle("J4_PSM1");
	J2_TOOL1 = sim.getObjectHandle("J5_PSM1");
	J3_dx_TOOL1 = sim.getObjectHandle("J6_PSM1");
    J3_sx_TOOL1 = sim.getObjectHandle("J7_PSM1");
    cwShapeHandle = sim.getObjectHandle("CW_respondable");
    cwShapeMass = sim.getShapeMass(cwShapeHandle)

    -- Graph handles
    graph_q = sim.getObjectHandle("qGraph");
    graph_qd = sim.getObjectHandle("qdotGraph");
    graph_orientation_error = sim.getObjectHandle("Orientation_Error");
    graph_position_error = sim.getObjectHandle("Position_Error");
    graph_joint_reg_error = sim.getObjectHandle("Joint_Error");
    graph_Traj = sim.getObjectHandle("eeTraj");

    graphJ1_tauComp = sim.getObjectHandle("J1_TorqueComparison");
    graphJ2_tauComp = sim.getObjectHandle("J2_TorqueComparison");
    graphJ3_tauComp = sim.getObjectHandle("J3_TorqueComparison");

    -- Plot streams
    --fcwm_meas_str = sim.addGraphStream(graph_Force_J3,'fcw_mir_meas','N',0,{0,1,0})
    q1meas_str = sim.addGraphStream(graph_q,'q1_meas','N',0,{1,0,0})
    q2meas_str = sim.addGraphStream(graph_q,'q2_meas','N',0,{0,1,0})
    q3meas_str = sim.addGraphStream(graph_q,'q3_meas','N',0,{0,0,1})
    q1cmd_str = sim.addGraphStream(graph_q,'q1_des','N',0,{1,0,1})
    q2cmd_str = sim.addGraphStream(graph_q,'q2_des','N',0,{1,1,0})
    q3cmd_str = sim.addGraphStream(graph_q,'q3_des','N',0,{0,1,1})
    qd1cmd_str = sim.addGraphStream(graph_qd,'qd1cmd','rad/s',0,{1,0,0})
    qd2cmd_str = sim.addGraphStream(graph_qd,'qd2cmd','rad/s',0,{0,1,0})
    qd3cmd_str = sim.addGraphStream(graph_qd,'qd3cmd','rad/s',0,{0,0,1})
    qd4cmd_str = sim.addGraphStream(graph_qd,'qd4cmd','rad/s',0,{1,0,1})
    qd5cmd_str = sim.addGraphStream(graph_qd,'qd5cmd','rad/s',0,{1,1,0})
    qd6cmd_str = sim.addGraphStream(graph_qd,'qd6cmd','rad/s',0,{0,1,1})
    err_x_str = sim.addGraphStream(graph_position_error,'err_x','m',0,{1,0,0})
    err_y_str = sim.addGraphStream(graph_position_error,'err_y','m',0,{0,1,0})
    err_z_str = sim.addGraphStream(graph_position_error,'err_z','m',0,{0,0,1})
    err_a_str = sim.addGraphStream(graph_orientation_error,'err_a','rad',0,{1,0,0})
    err_b_str = sim.addGraphStream(graph_orientation_error,'err_b','rad',0,{0,1,0})
    err_g_str = sim.addGraphStream(graph_orientation_error,'err_g','rad',0,{0,0,1})
    err_q1_str = sim.addGraphStream(graph_joint_reg_error,'err_q1','rad',0,{1,0,0})
    err_q2_str = sim.addGraphStream(graph_joint_reg_error,'err_q2','rad',0,{0,1,0})
    err_q3_str = sim.addGraphStream(graph_joint_reg_error,'err_q3','rad',0,{0,0,1})

    xEE = sim.addGraphStream(graph_Traj,'xEE','m',0,{1,0,0})
    yEE = sim.addGraphStream(graph_Traj,'yEE','m',0,{0,1,0})
    zEE = sim.addGraphStream(graph_Traj,'zEE','m',0,{0,0,1})
    sim.addGraphCurve(graph_Traj,'EE-Trajectory',3,{xEE,yEE,zEE},{0,0,0},'[m]',2,{0,0,1})

    tau1mod_comp_str = sim.addGraphStream(graphJ1_tauComp,'tau1_mod','Nm',0,{0,0,1})
    tau2mod_comp_str = sim.addGraphStream(graphJ2_tauComp,'tau2_mod','Nm',0,{0,0,1})
    f3mod_comp_str = sim.addGraphStream(graphJ3_tauComp,'f3_mod','N',0,{0,0,1})

    tau1meas_comp_str = sim.addGraphStream(graphJ1_tauComp,'tau1_meas','Nm',0,{1,0,0})
    tau2meas_comp_str = sim.addGraphStream(graphJ2_tauComp,'tau2_meas','Nm',0,{1,0,0})
    f3meas_comp_str = sim.addGraphStream(graphJ3_tauComp,'f3_meas','N',0,{1,0,0})

    tau1cmd_comp_str = sim.addGraphStream(graphJ1_tauComp,'tau1_cmd','Nm',0,{0,1,0})
    tau2cmd_comp_str = sim.addGraphStream(graphJ2_tauComp,'tau2_cmd','Nm',0,{0,1,0})
    f3cmd_comp_str = sim.addGraphStream(graphJ3_tauComp,'f3_cmd','N',0,{0,1,0})

    --initialization
    sim.setFloatSignal("simulationTimeStep",Ts)
    t = 0.0
    qdot = {0, 0, 0, 0, 0, 0}
    qdot_msr = {0,0,0,0,0,0,0}    
    qdd_msr = {0,0,0,0,0,0,0}    
    error = {0,0,0,0,0,0}
    err_q = {0,0,0,0,0,0,0}
    friction = {0,0,0}
    pee = {0,0,0}
    qmsr = {nil,nil,nil,nil,nil,nil}
    qmsr_prev = {nil,nil,nil,nil,nil,nil}
    qcmd = {nil,nil,nil,nil,nil,nil}
    taucmd = {nil,nil,nil,nil,nil,nil,nil}    
    taumsr = {nil,nil,nil,nil,nil,nil,nil} 
    taumsr_prev = {0,0,0,0,0,0,0} 
    tau_mod = {nil,nil,nil}
    tau_modLag = {nil,nil,nil}
    f_cw = 0.0
    q_cw = 0.0
    qdot_cw = 0.0
    q3 = 0.0
    qdes = {nil,nil,nil}
    tauMsr = {0,0,0}

    qmsr[1] = sim.getJointPosition(J1_PSM)
    qmsr[2] = sim.getJointPosition(J2_PSM)
    qmsr[3] = sim.getJointPosition(J3_PSM)
    qmsr[4] = sim.getJointPosition(J1_TOOL1)
    qmsr[5] = sim.getJointPosition(J2_TOOL1)
    qmsr[6] = sim.getJointPosition(J3_dx_TOOL1)

    for i = 1, 6 do
        qcmd[i] = qmsr[i]
        qmsr_prev[i] = qmsr[i]
    end

    
    -- Enable/Disable low-level joint position control loop on the virtual joints
    sim.setObjectInt32Parameter(J1_PSM,sim.jointintparam_ctrl_enabled,lowlevelCtrlLoopEnabled)
    sim.setObjectInt32Parameter(J2_PSM,sim.jointintparam_ctrl_enabled,lowlevelCtrlLoopEnabled)
    sim.setObjectInt32Parameter(J3_PSM,sim.jointintparam_ctrl_enabled,lowlevelCtrlLoopEnabled)

end

function sysCall_actuation()
    -- put your actuation code here
    sim.setFloatSignal("dynamicCtrlFlagSignal",dynamicCtrlFlag)

    -- Integration of joint velocity commands
    if qdot[1] ~= nil then
        qcmd[1] = qcmd[1] + qdot[1] * Ts
        qcmd[2] = qcmd[2] + qdot[2] * Ts
        qcmd[3] = qcmd[3] + qdot[3] * Ts
        qcmd[4] = qcmd[4] + qdot[4] * Ts
        qcmd[5] = qcmd[5] + qdot[5] * Ts
        qcmd[6] = qcmd[6] + qdot[6] * Ts
    end

    -- Set control inputs based on the selected control mode
    if dynamicCtrlFlag == 0 then -- kinematic control (J1, J2, J3: joint velocity/position commands)

        if qdot[1] ~= nil then
        
            if lowlevelCtrlLoopEnabled == 1 then
                sim.setJointTargetPosition(J1_PSM, qcmd[1])
                sim.setJointTargetPosition(J2_PSM, qcmd[2])
                sim.setJointTargetPosition(J3_PSM, qcmd[3]) 
            else 
                sim.setJointTargetVelocity(J1_PSM, qdot[1])
                sim.setJointTargetVelocity(J2_PSM, qdot[2])
                sim.setJointTargetVelocity(J3_PSM, qdot[3])            
            end
        end
        
        -- Set the lat three joint positions (tool wrist: J4, J5, J6-J7)
        sim.setJointPosition(J1_TOOL1, qcmd[4])
        sim.setJointPosition(J2_TOOL1, qcmd[5])
        sim.setJointPosition(J3_dx_TOOL1, qcmd[6])
        sim.setJointPosition(J3_sx_TOOL1, qcmd[6])

    else -- dynamic control (J1, J2, J3: joint torque commands)

        if taucmd[1] ~= nil then
            sim.setJointMaxForce(J1_PSM,math.abs(taucmd[1]))
            sim.setJointMaxForce(J2_PSM,math.abs(taucmd[2]))
            sim.setJointMaxForce(J3_PSM,math.abs(taucmd[3]))
        
            if lowlevelCtrlLoopEnabled == 1 then
                sim.setJointTargetPosition(J1_PSM,1e3*sign(taucmd[1]))
                sim.setJointTargetPosition(J2_PSM,1e3*sign(taucmd[2]))
                sim.setJointTargetPosition(J3_PSM,1e3*sign(taucmd[3]))
            else 
                sim.setJointTargetVelocity(J1_PSM,1e3*sign(taucmd[1]))
                sim.setJointTargetVelocity(J2_PSM,1e3*sign(taucmd[2]))
                sim.setJointTargetVelocity(J3_PSM,1e3*sign(taucmd[3]))
            end
        end
        
    end

    -- Set the mirror J3 joint force from the weight force of the counterweight
    if f_cw ~= nil then
        sim.setJointMaxForce(J3m_PSM,math.abs(f_cw))
        sim.setJointTargetVelocity(J3m_PSM,1e3*sign(f_cw))
    end

    -- Move the counterweight accordingly with the tool
    if dynamicCtrlFlag then
        if q3 ~= nil and q3 ~= 0.0 then
            sim.setJointTargetPosition(JCW_PSM,q3)
        end
    else
        if qmsr[3] ~= nil then
            sim.setJointTargetPosition(JCW_PSM,q3)
        end
    end

    -- Plots
    -- joint position measurements
    sim.setGraphStreamValue(graph_q,q1meas_str,qmsr[1])
    sim.setGraphStreamValue(graph_q,q2meas_str,qmsr[2])
    sim.setGraphStreamValue(graph_q,q3meas_str,qmsr[3])

    if qdes[1] ~= nil then
        sim.setGraphStreamValue(graph_q,q1cmd_str,qdes[1])
        sim.setGraphStreamValue(graph_q,q2cmd_str,qdes[2])
        sim.setGraphStreamValue(graph_q,q3cmd_str,qdes[3])
    end

    -- joint velocity commands/measurements
    if qdot[1] ~= nil then
        sim.setGraphStreamValue(graph_qd,qd1cmd_str,qdot[1])
        sim.setGraphStreamValue(graph_qd,qd2cmd_str,qdot[2])
        sim.setGraphStreamValue(graph_qd,qd3cmd_str,qdot[3])
        sim.setGraphStreamValue(graph_qd,qd4cmd_str,qdot[4])
        sim.setGraphStreamValue(graph_qd,qd5cmd_str,qdot[5])
        sim.setGraphStreamValue(graph_qd,qd6cmd_str,qdot[6])
    end

    -- joint torque measurements
    if taumsr[1] ~= nil then
        sim.setGraphStreamValue(graphJ1_tauComp,tau1meas_comp_str,taumsr[1])
        sim.setGraphStreamValue(graphJ2_tauComp,tau2meas_comp_str,taumsr[2])
        sim.setGraphStreamValue(graphJ3_tauComp,f3meas_comp_str,taumsr[3])

        if dynamicCtrlFlag == 1 then
        end
    end

    -- joint torque commands
    if taucmd[1] ~= nil then 
        sim.setGraphStreamValue(graphJ1_tauComp,tau1cmd_comp_str,-taucmd[1])
        sim.setGraphStreamValue(graphJ2_tauComp,tau2cmd_comp_str,-taucmd[2])
        sim.setGraphStreamValue(graphJ3_tauComp,f3cmd_comp_str,-taucmd[3])
    end

    -- joint model torques
    if tau_mod[1] ~= nil then
        sim.setGraphStreamValue(graphJ1_tauComp,tau1mod_comp_str,tau_mod[1])
        sim.setGraphStreamValue(graphJ2_tauComp,tau2mod_comp_str,tau_mod[2])
        sim.setGraphStreamValue(graphJ3_tauComp,f3mod_comp_str,tau_mod[3])
    end
    
    -- Cartesian EE position
    if pee[1] ~= nil and pee[1] ~= 0.0 then
        sim.setGraphStreamValue(graph_Traj,xEE,pee[1])
        sim.setGraphStreamValue(graph_Traj,yEE,pee[2])
        sim.setGraphStreamValue(graph_Traj,zEE,pee[3])
    end

    -- Cartesian error
    if error[1] ~= nil then 
        sim.setGraphStreamValue(graph_position_error,err_x_str,error[1])
        sim.setGraphStreamValue(graph_position_error,err_y_str,error[2])
        sim.setGraphStreamValue(graph_position_error,err_z_str,error[3])
        sim.setGraphStreamValue(graph_orientation_error,err_a_str,error[4])
        sim.setGraphStreamValue(graph_orientation_error,err_b_str,error[5])
        sim.setGraphStreamValue(graph_orientation_error,err_g_str,error[6])
    end

    -- Joint regulation error
    if err_q[1] ~= nil then
        sim.setGraphStreamValue(graph_joint_reg_error,err_q1_str,err_q[1])
        sim.setGraphStreamValue(graph_joint_reg_error,err_q2_str,err_q[2])
        sim.setGraphStreamValue(graph_joint_reg_error,err_q3_str,err_q[3])
    end

    
end

function sysCall_sensing()
    -- put your sensing code here

    -- Get simulation time 
    t = sim.getSimulationTime()

    -- Get the EE position
    pee[1] = sim.getFloatSignal('x_ee')
    pee[2] = sim.getFloatSignal('y_ee')
    pee[3] = sim.getFloatSignal('z_ee')

    f_cw = sim.getJointForce(JCW_PSM)
    q_cw = sim.getJointPosition(JCW_PSM)
    qdot_cw = sim.getJointVelocity(JCW_PSM)

    -- Get joint model torques
    tau_mod[1] = sim.getFloatSignal("Torque_J1")
    tau_mod[2] = sim.getFloatSignal("Torque_J2")
    tau_mod[3] = sim.getFloatSignal("Force_J3")

    -- Get joint position measurements 
    qmsr[1] = sim.getJointPosition(J1_PSM)
    qmsr[2] = sim.getJointPosition(J2_PSM)
    qmsr[3] = sim.getJointPosition(J3_PSM)
    q3 = sim.getJointPosition(J3_PSM)
    q3mirror = sim.getJointPosition(J3m_PSM)

    -- Get reference joint position trajectory
    qdes[1] = sim.getFloatSignal("qLdes_1")
    qdes[2] = sim.getFloatSignal("qLdes_2")
    qdes[3] = sim.getFloatSignal("qLdes_3")

    -- Get external joint velocity signals
    qdot[1] = sim.getFloatSignal("qLdot_1")
    qdot[2] = sim.getFloatSignal("qLdot_2")
    qdot[3] = sim.getFloatSignal("qLdot_3")
    qdot[4] = sim.getFloatSignal("qLdot_4")
    qdot[5] = sim.getFloatSignal("qLdot_5")
    qdot[6] = sim.getFloatSignal("qLdot_6")

    -- Get exterenal 6D Cartesian error signals
    error[1] = sim.getFloatSignal("error_x")
    error[2] = sim.getFloatSignal("error_y")
    error[3] = sim.getFloatSignal("error_z")
    error[4] = sim.getFloatSignal("error_alfa")
    error[5] = sim.getFloatSignal("error_beta")
    error[6] = sim.getFloatSignal("error_gamma")

    -- Get external joint regulation error
    err_q[1] = sim.getFloatSignal("error_q1")
    err_q[2] = sim.getFloatSignal("error_q2")
    err_q[3] = sim.getFloatSignal("error_q3")

    -- Get joint force measurements
    taumsr[1] = sim.getFloatSignal("tauLMsr_1")
    taumsr[2] = sim.getFloatSignal("tauLMsr_2")
    taumsr[3] = sim.getFloatSignal("tauLMsr_3")

    -- Get exteranl joint torque signals    
    --if dynamicCtrlFlag == 1 then 
        taucmd[1] = sim.getFloatSignal("tauL_1")
        taucmd[2] = sim.getFloatSignal("tauL_2")
        taucmd[3] = sim.getFloatSignal("tauL_3")
        taucmd[4] = sim.getFloatSignal("tauL_4")
        taucmd[5] = sim.getFloatSignal("tauL_5")
        taucmd[6] = sim.getFloatSignal("tauL_6")
        taucmd[7] = sim.getFloatSignal("tauL_7")
    --end

    tau_modLag[1] = sim.getFloatSignal("tauModLag_1")
    tau_modLag[2] = sim.getFloatSignal("tauModLag_2")
    tau_modLag[3] = sim.getFloatSignal("tauModLag_3")


    -- Get external friction signals  
    friction[1] = sim.getFloatSignal("CF_joint1")
    friction[2] = sim.getFloatSignal("CF_joint2")
    friction[3] = sim.getFloatSignal("CF_joint3")

    for i = 1, 6 do
        qdot_msr[i] = (qmsr[i] - qmsr_prev[i])/Ts
        qmsr_prev[i] = qmsr[i]
    end

end

function sysCall_cleanup()
    -- do some clean-up here
    qdot = {0, 0, 0, 0, 0, 0}
    q    = {0, 0, 0, 0, 0, 0}
    qcmd = {0, 0, 0, 0, 0, 0}

    sim.clearFloatSignal('x_ee')
    sim.clearFloatSignal('y_ee')
    sim.clearFloatSignal('z_ee')
    
    sim.clearFloatSignal("simulationTimeStep")

    sim.clearFloatSignal("Torque_J1")
    sim.clearFloatSignal("Torque_J2")
    sim.clearFloatSignal("Force_J3")
    
    sim.clearFloatSignal("qLdot_1")
    sim.clearFloatSignal("qLdot_2")
    sim.clearFloatSignal("qLdot_3")
    sim.clearFloatSignal("qLdot_4")
    sim.clearFloatSignal("qLdot_5")
    sim.clearFloatSignal("qLdot_6")

    sim.clearFloatSignal("error_x")
    sim.clearFloatSignal("error_y")
    sim.clearFloatSignal("error_z")
    sim.clearFloatSignal("error_alfa")
    sim.clearFloatSignal("error_beta")
    sim.clearFloatSignal("error_gamma")

    sim.clearFloatSignal("tauL_1")
    sim.clearFloatSignal("tauL_2")
    sim.clearFloatSignal("tauL_3")
    sim.clearFloatSignal("tauL_4")
    sim.clearFloatSignal("tauL_5")
    sim.clearFloatSignal("tauL_6")
    sim.clearFloatSignal("tauL_7")
end

-- See the user manual or the available code snippets for additional callback functions and details
function sign(x)
    if(x > 0) then
        return 1
    elseif(x < 0) then
        return -1
    else
        return 0
    end
end