function [clusters,energy]= Gcut_update(R_inv,smoothcost,x_mesh,clusters,M,depth)
[N,group_num]=size(R_inv);
h= GCO_Create(N,group_num);
datacost=((10^3))*R_inv';%set datacost
% smoothcost=(10^3)*smoothcost;
% datacost=datacost+ones(size(datacost));

%set neighbors
neighbors=zeros(N,N);
n_cost=(10^3)*1;
% n_cost=(10^3);

% x=contrast(1:end-1,:);
beta1=0;
beta2=0;
beta_num1=0;%3*length(x_mesh);
beta_num2=0;
diff1=zeros(size(x_mesh));%ones(length(x_mesh),3)*100;
diff2=zeros(size(x_mesh));
for iter=1:length(x_mesh)
    diff1(iter,1)=(10^(-10))+norm(M(x_mesh(iter,1),:)-M(x_mesh(iter,2),:));
%     diff1(iter,1)=abs(M(x_mesh(iter,1),clusters(x_mesh(iter,1)))-M(x_mesh(iter,2),clusters(x_mesh(iter,1))))+abs(M(x_mesh(iter,1),clusters(x_mesh(iter,2)))-M(x_mesh(iter,2),clusters(x_mesh(iter,2))));
    diff2(iter,1)=(10^(-10))+abs(depth(x_mesh(iter,1))-depth(x_mesh(iter,2)));
    beta1=beta1+diff1(iter,1);
    beta_num1=beta_num1+1;
%     if diff2(iter,1)<0.1
        beta2=beta2+diff2(iter,1);
        beta_num2=beta_num2+1;
%     end

%     diff1(iter,2)=abs(M(x_mesh(iter,1),clusters(x_mesh(iter,1)))-M(x_mesh(iter,3),clusters(x_mesh(iter,1))))+abs(M(x_mesh(iter,1),clusters(x_mesh(iter,3)))-M(x_mesh(iter,3),clusters(x_mesh(iter,3))));%
    diff1(iter,2)=(10^(-10))+norm(M(x_mesh(iter,1),:)-M(x_mesh(iter,3),:));
    diff2(iter,2)=(10^(-10))+abs(depth(x_mesh(iter,1))-depth(x_mesh(iter,3)));
    beta1=beta1+diff1(iter,2);
    beta_num1=beta_num1+1;
%     if diff2(iter,2)<0.1
        beta2=beta2+diff2(iter,2);
        beta_num2=beta_num2+1;
%     end

    diff1(iter,3)=(10^(-10))+norm(M(x_mesh(iter,2),:)-M(x_mesh(iter,3),:));
%     diff1(iter,3)=abs(M(x_mesh(iter,2),clusters(x_mesh(iter,2)))-M(x_mesh(iter,3),clusters(x_mesh(iter,2))))+abs(M(x_mesh(iter,2),clusters(x_mesh(iter,3)))-M(x_mesh(iter,3),clusters(x_mesh(iter,3))));
    diff2(iter,3)=(10^(-10))+abs(depth(x_mesh(iter,2))-depth(x_mesh(iter,3)));
    beta1=beta1+diff1(iter,3);
    beta_num1=beta_num1+1;
%     if diff2(iter,3)<0.1
        beta2=beta2+diff2(iter,3);
        beta_num2=beta_num2+1;
%     end

end
beta1=(2*beta_num1)/(beta1);%10^-2;%1;%
beta2=(2*beta_num2)/(beta2);
for iter=1:length(x_mesh)
%         if clusters(x_mesh(iter,1))~=clusters(x_mesh(iter,2))
        neighbors(x_mesh(iter,1),x_mesh(iter,2))=n_cost*(10*exp(-(beta1*diff1(iter,1)))+1*exp(-(beta2*diff2(iter,1))));
%         end
%         if clusters(x_mesh(iter,2))~=clusters(x_mesh(iter,3))
        neighbors(x_mesh(iter,2),x_mesh(iter,3))=n_cost*(10*exp(-(beta1*diff1(iter,3)))+1*exp(-(beta2*diff2(iter,3))));
%         end
%         if clusters(x_mesh(iter,1))~=clusters(x_mesh(iter,3))
        neighbors(x_mesh(iter,1),x_mesh(iter,3))=n_cost*(10*exp(-(beta1*diff1(iter,2)))+1*exp(-(beta2*diff2(iter,2))));
%         end
end

labelcost=zeros(group_num,1);
for iter=1:group_num
    labelcost(iter)=sum(clusters==iter)/N*(10^3);%(10^5)*sum(init_clusters==iter)/N;%
end


GCO_SetDataCost(h,datacost);
GCO_SetSmoothCost(h,smoothcost);
GCO_SetNeighbors(h,neighbors);
% GCO_SetLabelCost(h,labelcost);
GCO_Expansion(h);
clusters=GCO_GetLabeling(h);

energy=zeros(size(R_inv));
for iter1=1:N
    for iter2=1:group_num
        energy(iter1,iter2)=GCO_GetEnergy(h,iter1,iter2);
    end
end

GCO_Delete(h);


end