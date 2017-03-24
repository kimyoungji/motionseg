function plot_corr(X,clusters)
    [M,N]=size(X);
    m=M/4-1;
    color=colormap('prism');
    for iter=1:m
        for iter2=1:N
            plot3([X(4*(iter-1)+1,iter2),X(4*(iter)+1,iter2)],[X(4*(iter-1)+2,iter2),X(4*(iter)+2,iter2)],[X(4*(iter-1)+3,iter2),X(4*(iter)+3,iter2)],'-o','Color',color(clusters(iter2),:),'MarkerSize',2);
            hold on;
        end
    end
end