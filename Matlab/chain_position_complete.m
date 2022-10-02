function chain = chain_position_complete(DH,Vrep_T_DH0)

        chain = FramesRotationFromDH(DH,'modified'); %from DH-table to T (homogeneus transformation)
      
        chain(:,:,1) = Vrep_T_DH0*chain(:,:,1);
        chain(:,:,2)=  chain(:,:,1)*chain(:,:,2);
        chain(:,:,3)=  chain(:,:,2)*chain(:,:,3);
        chain(:,:,4)=  chain(:,:,3)*chain(:,:,4);
        chain(:,:,5)=  chain(:,:,3)*chain(:,:,5);
        chain(:,:,6)=  chain(:,:,4)*chain(:,:,6);
        chain(:,:,7)=  chain(:,:,4)*chain(:,:,7);
        chain(:,:,8)=  chain(:,:,6)*chain(:,:,8);
        chain(:,:,9)=  chain(:,:,2)*chain(:,:,9);
        chain(:,:,10)=  chain(:,:,8)*chain(:,:,10);
        chain(:,:,11)=  chain(:,:,10)*chain(:,:,11);
        chain(:,:,12)=  chain(:,:,11)*chain(:,:,12);
        chain(:,:,13)=  chain(:,:,11)*chain(:,:,13);
        chain(:,:,14)=  chain(:,:,13)*chain(:,:,14);
        
end