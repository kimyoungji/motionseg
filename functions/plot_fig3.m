%test

load('./synthetic/plane_rot1/pt_5.mat');
% N=size(a,2);

    plot3(a(1,:)',a(2,:)',a(3,:)','k.');
    xhandle=xlabel('x(m)');
    set(xhandle,'Fontsize',30);
    set(xhandle,'Fontangle','italic');
    set(xhandle,'Fontname','Times new roman');
    yhandle=ylabel('y(m)');
    set(yhandle,'Fontsize',30);
    set(yhandle,'Fontangle','italic');
    set(yhandle,'Fontname','Times new roman');
    zhandle=zlabel('z(m)');
    set(zhandle,'Fontsize',30);
    set(zhandle,'Fontangle','italic');
    set(zhandle,'Fontname','Times new roman');
    set(gca,'Fontname','Times new roman');
    set(gca,'Fontsize',15);
%     grid on;
    axis equal;
%     axis([90,110,40,100,-20,20]);
