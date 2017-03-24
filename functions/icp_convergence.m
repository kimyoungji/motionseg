function is_converge=icp_convergence(T)
    m=T(1:3,1:3);
    t=T(1:3,4);
    r = SpinCalc('DCMtoEA123',m,0.0001,0);
    if abs(r(1))>180
        r(1)=360-abs(r(1));
    end
    if abs(r(2))>180
        r(2)=360-abs(r(2));
    end
    if abs(r(3))>180
        r(3)=360-abs(r(3));
    end
    r=r*pi/180;
    is_converge=norm([t;r']);
%     if norm([t;r'])<0.005
%         is_converge=1;
%     else
%         is_converge=0;
%     end
end