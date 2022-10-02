function DH0 = checkRefDH0(clientID,sim)

    [~,handle]=sim.simxGetObjectHandle(clientID,'DH0_ref',sim.simx_opmode_blocking);
 	[~,temp_]=sim.simxGetObjectOrientation(clientID,handle,-1,sim.simx_opmode_blocking);
 
    temp = [double(temp_(1)),double(temp_(2)),double(temp_(3))]; 
    Rx = [   1   0   0 ;  0   cos(temp(1))  -sin(temp(1)) ;    0   sin(temp(1))   cos(temp(1)) ];
    Ry = [  cos(temp(2))  0  sin(temp(2)) ;   0  1  0  ; -sin(temp(2))  0   cos(temp(2))   ];
    Rz = [  cos(temp(3))  -sin(temp(3)) 0 ; sin(temp(3)) cos(temp(3)) 0;  0   0  1 ];
    R_xyz = Rx*Ry*Rz;
    
    [~,temp_]=sim.simxGetObjectPosition(clientID,handle,-1,sim.simx_opmode_blocking);
    position = [double(temp_(1)),double(temp_(2)),double(temp_(3))];
    
    DH0 = [  R_xyz   [position(1) position(2) position(3)]' ; 0  0  0  1];
     
end