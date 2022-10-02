function sysCall_init()
    -- do some initialization here
    
    -- Getting simulation time step
    Ts = sim.getSimulationTimeStep()

    -- flags
    dynamicCtrlFlag = 0
    lowlevelCtrlLoopEnabled = 0

    -- Object handles
	J1_PSM2 = sim.getObjectHandle("J1_PSM2");
	J2_PSM2 = sim.getObjectHandle("J2_PSM2");
	J3_PSM2 = sim.getObjectHandle("J3_PSM2");
	J3m_PSM2 = sim.getObjectHandle("J3mirror_PSM2");
	JCW_PSM2 = sim.getObjectHandle("J31_PSM2");
	J1_TOOL2 = sim.getObjectHandle("J4_PSM2");
	J2_TOOL2 = sim.getObjectHandle("J5_PSM2");
	J3_dx_TOOL2 = sim.getObjectHandle("J6_PSM2");
    J3_sx_TOOL2 = sim.getObjectHandle("J7_PSM2");
    cwShapeHandle2 = sim.getObjectHandle("CW_respondable2");
    cwShapeMass2 = sim.getShapeMass(cwShapeHandle2)

    --initialization
    sim.setFloatSignal("simulationTimeStep",Ts)
    t = 0.0
    qdot_r = {0, 0, 0, 0, 0, 0}
    qdot_r_msr = {0,0,0,0,0,0,0}    
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
    qdot_r_cw = 0.0
    q3 = 0.0
    qdes = {nil,nil,nil}
    tauMsr = {0,0,0}

    qmsr[1] = sim.getJointPosition(J1_PSM2)
    qmsr[2] = sim.getJointPosition(J2_PSM2)
    qmsr[3] = sim.getJointPosition(J3_PSM2)
    qmsr[4] = sim.getJointPosition(J1_TOOL2)
    qmsr[5] = sim.getJointPosition(J2_TOOL2)
    qmsr[6] = sim.getJointPosition(J3_dx_TOOL2)

    for i = 1, 6 do
        qcmd[i] = qmsr[i]
        qmsr_prev[i] = qmsr[i]
    end

    
    -- Enable/Disable low-level joint position control loop on the virtual joints
    sim.setObjectInt32Parameter(J1_PSM2,sim.jointintparam_ctrl_enabled,lowlevelCtrlLoopEnabled)
    sim.setObjectInt32Parameter(J2_PSM2,sim.jointintparam_ctrl_enabled,lowlevelCtrlLoopEnabled)
    sim.setObjectInt32Parameter(J3_PSM2,sim.jointintparam_ctrl_enabled,lowlevelCtrlLoopEnabled)

end

function sysCall_actuation()
    -- put your actuation code here
    sim.setFloatSignal("dynamicCtrlFlagSignal",dynamicCtrlFlag)

    -- Integration of joint velocity commands
    if qdot_r[1] ~= nil then
        qcmd[1] = qcmd[1] + qdot_r[1] * Ts
        qcmd[2] = qcmd[2] + qdot_r[2] * Ts
        qcmd[3] = qcmd[3] + qdot_r[3] * Ts
        qcmd[4] = qcmd[4] + qdot_r[4] * Ts
        qcmd[5] = qcmd[5] + qdot_r[5] * Ts
        qcmd[6] = qcmd[6] + qdot_r[6] * Ts
    end

    -- Set control inputs based on the selected control mode
    if dynamicCtrlFlag == 0 then -- kinematic control (J1, J2, J3: joint velocity/position commands)

        if qdot_r[1] ~= nil then
        
            if lowlevelCtrlLoopEnabled == 1 then
                sim.setJointTargetPosition(J1_PSM2, qcmd[1])
                sim.setJointTargetPosition(J2_PSM2, qcmd[2])
                sim.setJointTargetPosition(J3_PSM2, qcmd[3]) 
            else 
                sim.setJointTargetVelocity(J1_PSM2, qdot_r[1])
                sim.setJointTargetVelocity(J2_PSM2, qdot_r[2])
                sim.setJointTargetVelocity(J3_PSM2, qdot_r[3])            
            end
        end
        
        -- Set the lat three joint positions (tool wrist: J4, J5, J6-J7)
        sim.setJointPosition(J1_TOOL2, qcmd[4])
        sim.setJointPosition(J2_TOOL2, qcmd[5])
        sim.setJointPosition(J3_dx_TOOL2, qcmd[6])
        sim.setJointPosition(J3_sx_TOOL2, qcmd[6])

    else -- dynamic control (J1, J2, J3: joint torque commands)

        if taucmd[1] ~= nil then
            sim.setJointMaxForce(J1_PSM2,math.abs(taucmd[1]))
            sim.setJointMaxForce(J2_PSM2,math.abs(taucmd[2]))
            sim.setJointMaxForce(J3_PSM2,math.abs(taucmd[3]))
        
            if lowlevelCtrlLoopEnabled == 1 then
                sim.setJointTargetPosition(J1_PSM2,1e3*sign(taucmd[1]))
                sim.setJointTargetPosition(J2_PSM2,1e3*sign(taucmd[2]))
                sim.setJointTargetPosition(J3_PSM2,1e3*sign(taucmd[3]))
            else 
                sim.setJointTargetVelocity(J1_PSM2,1e3*sign(taucmd[1]))
                sim.setJointTargetVelocity(J2_PSM2,1e3*sign(taucmd[2]))
                sim.setJointTargetVelocity(J3_PSM2,1e3*sign(taucmd[3]))
            end
        end
        
    end

    -- Set the mirror J3 joint force from the weight force of the counterweight
    if f_cw ~= nil then
        sim.setJointMaxForce(J3m_PSM2,math.abs(f_cw))
        sim.setJointTargetVelocity(J3m_PSM2,1e3*sign(f_cw))
    end

    -- Move the counterweight accordingly with the tool
    if dynamicCtrlFlag then
        if q3 ~= nil and q3 ~= 0.0 then
            sim.setJointTargetPosition(JCW_PSM2,q3)
        end
    else
        if qmsr[3] ~= nil then
            sim.setJointTargetPosition(JCW_PSM2,q3)
        end
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

    f_cw = sim.getJointForce(JCW_PSM2)
    q_cw = sim.getJointPosition(JCW_PSM2)
    qdot_r_cw = sim.getJointVelocity(JCW_PSM2)

    -- Get joint model torques
    tau_mod[1] = sim.getFloatSignal("Torque_J1")
    tau_mod[2] = sim.getFloatSignal("Torque_J2")
    tau_mod[3] = sim.getFloatSignal("Force_J3")

    -- Get joint position measurements 
    qmsr[1] = sim.getJointPosition(J1_PSM2)
    qmsr[2] = sim.getJointPosition(J2_PSM2)
    qmsr[3] = sim.getJointPosition(J3_PSM2)
    q3 = sim.getJointPosition(J3_PSM2)
    q3mirror = sim.getJointPosition(J3m_PSM2)

    -- Get reference joint position trajectory
    qdes[1] = sim.getFloatSignal("qRdes_1")
    qdes[2] = sim.getFloatSignal("qRdes_2")
    qdes[3] = sim.getFloatSignal("qRdes_3")

    -- Get external joint velocity signals
    qdot_r[1] = sim.getFloatSignal("qRdot_1")
    qdot_r[2] = sim.getFloatSignal("qRdot_2")
    qdot_r[3] = sim.getFloatSignal("qRdot_3")
    qdot_r[4] = sim.getFloatSignal("qRdot_4")
    qdot_r[5] = sim.getFloatSignal("qRdot_5")
    qdot_r[6] = sim.getFloatSignal("qRdot_6")

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
    taumsr[1] = sim.getFloatSignal("tauRMsr_1")
    taumsr[2] = sim.getFloatSignal("tauRMsr_2")
    taumsr[3] = sim.getFloatSignal("tauRMsr_3")

    -- Get exteranl joint torque signals    
    taucmd[1] = sim.getFloatSignal("tauR_1")
    taucmd[2] = sim.getFloatSignal("tauR_2")
    taucmd[3] = sim.getFloatSignal("tauR_3")
    taucmd[4] = sim.getFloatSignal("tauR_4")
    taucmd[5] = sim.getFloatSignal("tauR_5")
    taucmd[6] = sim.getFloatSignal("tauR_6")
    taucmd[7] = sim.getFloatSignal("tauR_7")

    for i = 1, 6 do
        qdot_r_msr[i] = (qmsr[i] - qmsr_prev[i])/Ts
        qmsr_prev[i] = qmsr[i]
    end

end

function sysCall_cleanup()
    -- do some clean-up here
    qdot_r = {0, 0, 0, 0, 0, 0}
    q    = {0, 0, 0, 0, 0, 0}
    qcmd = {0, 0, 0, 0, 0, 0}
   
    sim.clearFloatSignal("simulationTimeStep")
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