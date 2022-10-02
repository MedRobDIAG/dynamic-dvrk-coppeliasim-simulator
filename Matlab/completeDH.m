function DH = completeDH(q)

        %q = [q(1) q(2) q(8) q(10) q(11) q(12) q(13)]; 

        %Parameters of DaVinci robot PSM
        l_RCC=431.8*10^-3;
        l_2L3=40.09*10^-3;
        l_2H1=144.54*10^-3;
        l_2H2=38.08*10^-3;
        l_c1=l_2H1+l_2H2;
        l_c2=-l_RCC+l_2H1;
        l_2L1=96*10^-3;
        l_2L2=516*10^-3;
        l_3=40.09*10^-3;
        l_p2y=0.91*10^-2;
        l_tool=416.2*10^-3;
        l_clip=0.91*10^-2;
       

        %----
     


        %  joint type | d_i | theta_i | a_i  | alpha_i                          
        %       1     | d1  |    q1   |  0   |   pi/2     
        DH = ...
        [
            [1,          0,  q(1)+pi/2,     0,   pi/2]
            [1,          0,  q(2)-pi/2,     0,  -pi/2]
            [1,          0,       pi/2, l_2L3,      0]
            [1,          0,  pi/2-q(2), l_2H1,      0]
            [1,          0,  pi/2-q(2), l_c1,       0]
            [1,          0,       q(2), l_2L2,      0]
            [1,          0,    q(2)+pi, l_2L1,      0]
            [2,  q(3)+l_c2,          0,   l_3,  -pi/2]
            [2,       q(3),          0, l_2L3,  -pi/2]
            [1,     l_tool,       q(4),     0,      0]
            [1,          0,  q(5)+pi/2,     0,   pi/2]
            [1,          0,  q(6)+pi/2,  l_p2y, -pi/2]
            [1,          0,  q(7)+pi/2,  l_p2y, -pi/2]
            [1,         -l_clip,         pi,      0,  -pi/2]    
        ];

        %----



end