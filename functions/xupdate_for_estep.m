function X_candi=xupdate_for_estep(X,dataMatrix,T_mat)
    N=size(X,2);
    group_num=size(T_mat,2);
    X_candi=repmat(X(end-3:end-1,:),group_num,1);
    tree=kdtree_build(dataMatrix');%create kd tree
    for iter=1:group_num
        X_tmp=T_mat{iter}*[X(end-7:end-5,:);ones(1,N)];
        [corr,dists]=kdtree_nearest_neighbor(tree,X_tmp(1:3,:)');
%         X_candi(3*(iter-1)+1:3*iter,dists<0.01)=dataMatrix(:,corr(dists<0.01));
        X_candi(3*(iter-1)+1:3*iter,:)=dataMatrix(:,corr);
    end
     kdtree_delete(tree);
end