function [Rot]=RotationMat(alpha,a)
    Rot=[];
    for iter=1:length(alpha)
        tmp=cos(alpha(iter))*eye(3,3)+(1-cos(alpha(iter)))*a*a'+sin(alpha(iter))*[0,-a(3),a(2);a(3),0,-a(1);-a(2),a(1),0];
        Rot=[Rot;tmp];
    end
end