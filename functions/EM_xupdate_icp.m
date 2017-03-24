function [T_mat,X_pre,X_nn,d_for_outliers,occ_idx]=EM_xupdate_icp(clusters,group_num,x,y)
%     group_num=numel(em_pack.rank_q);
    k=3;
    [~,N]=size(x);
    tree=kdtree_build(y');%create kd tree
    X_pre = zeros(3*group_num,N);
    X_nn = zeros(3*group_num,N);
    X_1 = zeros(3,N);
    X_2 = zeros(3,N);
    d_for_outliers=zeros(group_num,N);
    T_mat={};
    for g_idx=1:group_num
        clust_idx=find(clusters==g_idx);
        sampled_N=length(clust_idx);
        if sampled_N>10
            [T,~,~,~] = icp(x, y, tree);          
        else
%             em_pack.p_clusters(clust_idx)=0;
            T=eye(4);
        end
        T_mat{g_idx}=T;
        tmp=T*[x;ones(1,N)];
        X_pre(3*(g_idx-1)+1:3*g_idx,:)=tmp(1:3,:);
        [corr,dist]=kdtree_nearest_neighbor(tree,tmp(1:3,:)');
        d_for_outliers(g_idx,:)=dist;
        X_nn(3*(g_idx-1)+1:3*g_idx,:) = y(:,corr);
        X_1(:,clust_idx)=y(:,corr(clust_idx));
        X_2(:,clust_idx)=tmp(1:3,clust_idx);
    end
    kdtree_delete(tree);
    %% occlusion detection
    occ_idx=zeros(N,1);
    tree=kdtree_build(X_1');
    for iter=1:N
        [~,dist]=kdtree_k_nearest_neighbors(tree,X_1(:,iter)',2);
        if dist(2)<10^-3
            occ_idx(iter)=1;
        end
    end
    kdtree_delete(tree);
%     tree=kdtree_build(X_1');
%     [~,dist]=kdtree_nearest_neighbor(tree,X_2');
%     occ_idx=(dist>0.1);
%     kdtree_delete(tree);

        
end