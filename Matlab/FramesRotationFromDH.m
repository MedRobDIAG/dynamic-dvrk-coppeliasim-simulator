
function [ framesChain ] = FramesRotationFromDH(DH,Convention)

if nargin < 2
   Convention = 'classic';  % default 'classic D-H convention'
else
    if ~strcmp(Convention,'modified')  % 'modified D-H convention'
        Convention = 'classic';
    end
end

% framesChain = sym(zeros(4,4,size(DH,1)));

framesChain = zeros(4,4,size(DH,1));


for i=1:size(DH,1)
    
    if strcmp(Convention,'classic')
        framesChain(:,:,i) =  [cos(DH(i,3)), -cos(DH(i,5))*sin(DH(i,3)),  sin(DH(i,5))*sin(DH(i,3)),DH(i,4)*cos(DH(i,3));
                              sin(DH(i,3)),  cos(DH(i,5))*cos(DH(i,3)), -sin(DH(i,5))*cos(DH(i,3)),DH(i,4)*sin(DH(i,3));
                               0     ,            sin(DH(i,5))   ,          cos(DH(i,5))     ,     DH(i,2)        ;
                                     0     ,                 0         ,              0            ,         1         ];
    else
        framesChain(:,:,i) =  [  cos(DH(i,3))    ,     -sin(DH(i,3))     ,      0       ,      DH(i,4);
                                 cos(DH(i,5))*sin(DH(i,3)),  cos(DH(i,5))*cos(DH(i,3)), -sin(DH(i,5)),-sin(DH(i,5))*DH(i,2);
                                 sin(DH(i,5))*sin(DH(i,3)),  sin(DH(i,5))*cos(DH(i,3)), cos(DH(i,5)) , cos(DH(i,5))*DH(i,2);
                                    0               ,                 0         ,      0       ,         1         ];
    end
end


end