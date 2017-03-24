function [W,rr,inliers_1,C,dists,em_pack]=new_subspace(out_idx,dataMatrix,X,neighbors,em_pack)
    [M,~]=size(em_pack.X_mat);
    nn=M/4-1;
    k=4;
    tree=kdtree_build(dataMatrix');
    [T,~,inlier,~] = icp(X(1:3,out_idx), neighbors,dataMatrix, tree, 0.5);
    X_tmp=T*X(1:4,:);
    [corres_1,dist_1]=kdtree_nearest_neighbor(tree,X_tmp(1:3,:)');
    kdtree_delete(tree);
    
    inliers=out_idx;
    inliers_1=out_idx(inlier);
    new_X=[em_pack.X_mat(1:end-k,:);dataMatrix(:,corres_1);ones(1,numel(corres_1))];
    
    % find W
    X_norm=normc(new_X(:,inliers));
%     X_norm=(new_X(:,inliers));
    N_in=numel(inliers);
    A_a=zeros(M,M);
    for iter2=1:N_in
        A_a=A_a+X_norm(:,iter2)*X_norm(:,iter2)';
    end
    [~,D,W]=svd(A_a);
    s_values=diag(D);
    rr=find_rank(s_values(1:end-nn));
    if rr>4
        rr=4;
    end
            
    S=zeros(M-rr,M-rr);
    for iter2=1:N_in
        S=S+(W(:,rr+1:end)'*normc(new_X(:,inliers(iter2))))*(W(:,rr+1:end)'*normc(new_X(:,inliers(iter2))))';
    end
    C=S/N_in;
%     C=C+eye(M-rr)*(10^-5);
    C_inv=pinv(C);
    
     Mdists=zeros(numel(inliers),1);
     for iter2=1:numel(inliers)
        xx=W(:,rr+1:end)'*normc(em_pack.X_mat(:,inliers(iter2))); 
        Mdists(iter2)=(xx'*C_inv*xx);
     end
    
   
%     [idx1,idx2]=deleteoutliers(Mdists,0.05);
%     if mean(Mdists(idx2))<mean(Mdists(idx1))
%         idx1=idx2;
%     end
%     inliers=out_idx(idx1);
    dists=mean(dist_1(inliers_1));

end