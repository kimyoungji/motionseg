function plot_mesh(mesh,X)
    [N,~]=size(mesh);
    for iter=1:N
        a_1=mesh(iter,1);
        a_2=mesh(iter,2);
        a_3=mesh(iter,3);
        plot3([X(a_1,1);X(a_2,1)],[X(a_1,2);X(a_2,2)],[X(a_1,3);X(a_2,3)],'-or','MarkerFaceColor','k','MarkerSize',2);
        hold on;
        plot3([X(a_1,1);X(a_3,1)],[X(a_1,2);X(a_3,2)],[X(a_1,3);X(a_3,3)],'-or','MarkerFaceColor','k','MarkerSize',2);
        hold on;
        plot3([X(a_3,1);X(a_2,1)],[X(a_3,2);X(a_2,2)],[X(a_3,3);X(a_2,3)],'-or','MarkerFaceColor','k','MarkerSize',2);
        hold on;
    end
end