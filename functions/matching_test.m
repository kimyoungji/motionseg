%% matching test
addpath('./functions');
addpath('./functions/vbgm');
addpath('./functions/kdtree');
addpath('./functions/PairwiseMatching');
addpath('./functions/estimateRigidTransform');
addpath('./gco-v3.0/matlab');
addpath('./helper_functions');
addpath('./WOBJ_toolbox_Version2b');
%% initialize
N=100;% number of all points
d=3;% dimension of data [x;y;z]

ref_ang=0.2*pi/180;
dif_ang=-5*pi/180;
T=[cos(ref_ang),-sin(ref_ang),0,0;
    sin(ref_ang),cos(ref_ang),0,0;
    0,0,1,0;
    0,0,0,1];

diff=1;
ref=[10,10,21,10];
X1=make_dpoints(20,[ref(1),ref(1)+diff;ref(2),ref(2)+diff],[1,3,2]);
X2=T*X1;
X1=X1(:,201:400);
t=1:100;
X1=X1(:,t.*2);
X2=X2(:,101:300);
% X2=X2(:,101:400);%+randn(4,300)/1000;
%% matching
% select overlapping region
tree=kdtree_build(X2(1:3,:)');%create kd tree
% x_selected=X1;
[corr,~]=kdtree_nearest_neighbor(tree,X1(1:3,:)');
ransacCoef.minPtNum=3;
ransacCoef.iterNum=100;
ransacCoef.thInlrRatio=0.3;
ransacCoef.thDist=0.5;
[~, inlierIdx, corr] = ransac1( X1(1:3,:),X2(1:3,:),corr,tree,ransacCoef,@estimateRigidTransform,@estimateDistance);
x_selected=X1(:,inlierIdx);

%affinity matrix for assignment problem
N=size(x_selected,2);
in_idx=10;
nodes=zeros(1,N*in_idx);
labels=zeros(1,N*in_idx);
label_table=[];
A=zeros(1,N*in_idx);
for iter=1:N
    [candidates,dist]=kdtree_k_nearest_neighbors(tree,x_selected(1:3,iter)',in_idx);
    nodes(in_idx*(iter-1)+1:in_idx*iter)=ones(1,in_idx).*iter;
    
    for iter2=1:in_idx
        [idx,~]=find(label_table==candidates(iter2));
        if isempty(idx)
            idx=length(label_table)+1;
            label_table=[label_table;candidates(iter2)];
        end
        labels(in_idx*(iter-1)+iter2)=idx;
        A(in_idx*(iter-1)+iter2)=100*exp(-100*dist(iter2));
%         if iter==candidates(iter2)
%             A(in_idx*(iter-1)+iter2)=100;
%         end
    end
end
kdtree_delete(tree);

%construct pairwise terms
M=zeros(N*in_idx);
M(sub2ind(size(M), 1:in_idx*N,1:in_idx*N))=A;
dist=0;
for iter=1:N
        x11=repmat(x_selected(1:3,iter),(N-iter)*in_idx,in_idx);
        x12=repmat(x_selected(1:3,iter+1:end),in_idx,1);
        x10=x11-repmat(x12(:),1,in_idx);
        x1 = reshape(sum(reshape(x10,3, (N-iter)*in_idx*in_idx).^2),(N-iter)*in_idx,in_idx);
        x1=sqrt(x1);
        x21=repmat(X2(1:3,label_table(labels(in_idx*(iter-1)+1:iter*in_idx))),(N-iter)*in_idx,1 );
        x22=X2(1:3,label_table(labels(in_idx*(iter)+1:end)));
        x20=x21-repmat(x22(:),1,in_idx);
        x2 = reshape(sum(reshape(x20,3, (N-iter)*in_idx*in_idx).^2),(N-iter)*in_idx,in_idx);
        x2=sqrt(x2);
        dist=(x1-x2).*(x1-x2);
        sigma=0.01;
        tmp_dist=4.5-dist/(2*sigma^2);
        tmp_dist(dist>9*sigma^2)=0;
        tmp_dist(tmp_dist<0)=0;
        dist=tmp_dist;
        
        M(iter*in_idx+1:end,in_idx*(iter-1)+1:in_idx*(iter))=dist;
        M(in_idx*(iter-1)+1:in_idx*(iter),iter*in_idx+1:end)=dist';

end
% tic
% sol  = spectral_matching_yj(M, labels, nodes);
% toc
[sol, ~]  = spectral_matching_ipfp(M, labels, nodes);
matching=reshape(sol,in_idx,N);

for iter=1:N
    label_tmp=labels(in_idx*(iter-1)+1:in_idx*iter);
    if sum(matching(:,iter))==0
       neighbors(iter)=0; 
    else
       neighbors(iter)=label_table(label_tmp(matching(:,iter)==1));
    end
end

% neighbors=labels(matching==1);
%% plot

for iter=1:N
    idx=neighbors(iter);
    if idx==0
        plot3(x_selected(1,iter),x_selected(2,iter),x_selected(3,iter),'-or','MarkerFaceColor','k','MarkerSize',5);
        hold on;
    else
        plot3([x_selected(1,iter);X2(1,idx)],[x_selected(2,iter);X2(2,idx)],[x_selected(3,iter);X2(3,idx)],'-or','MarkerFaceColor','k','MarkerSize',5);
        hold on;
    end
end

hold on;
for iter=1:size(X2,2)
    plot3(X2(1,iter),X2(2,iter),X2(3,iter),'-ob','MarkerFaceColor','k','MarkerSize',2);
    hold on;
end

hold on;
for iter=1:size(X1,2)
    plot3(X1(1,iter),X1(2,iter),X1(3,iter),'-ok','MarkerFaceColor','k','MarkerSize',2);
    hold on;
end
grid on;
axis equal;
