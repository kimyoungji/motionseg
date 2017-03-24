function normals=PointNormal(x)
    N=size(x,1);
    normals=zeros(N,3);
    tree=kdtree_build(x);%create kd tree
    for iter=1:N
     idx=kdtree_k_nearest_neighbors( tree, x(iter,:), 10);
     normals(iter,:)=normnd(x(idx,:));
    end
    kdtree_delete(tree);
end