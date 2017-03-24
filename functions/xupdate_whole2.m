function [X,neighbors,em_pack]=xupdate_whole2(X,neighbors,clusters,group_num,em_pack)
    [M,N]=size(em_pack.X_mat);
    m=M/4-1;
    for g_idx=1:group_num
        em_pack.T_pre{g_idx}=em_pack.T_cur{g_idx}*em_pack.T_pre{g_idx};
        em_pack.T_cur{g_idx}=eye(4);
    end
    %% initial ICP
%     for g_idx=1:group_num
%         [c_idx,~]=find(clusters==g_idx);
%         for iter=1:m-1
%             tree=kdtree_build(em_pack.D_mat(3*(iter-1)+1:3*iter,:)');
%             if numel(c_idx)>3
%                 [T,~,~,~] = icp(em_pack.X_mat(4*(iter-1)+1:4*(iter-1)+3,c_idx), em_pack.Neighbors(iter,c_idx),em_pack.D_mat(3*(iter-1)+1:3*iter,:), tree, 0.5);
%             else
%                 T=eye(4);
%             end
%             X_tmp=T*[em_pack.X_mat(4*(iter-1)+1:4*(iter-1)+3,:);ones(1,N)];
%             [corr,~]=kdtree_nearest_neighbor(tree,X_tmp(1:3,:)');
%             em_pack.X_mat(4*(iter)+1:4*(iter)+3,c_idx)=em_pack.D_mat(3*(iter-1)+1:3*iter,corr(c_idx));
%             if numel(c_idx)>3
%                 em_pack.Neighbors(iter,c_idx)=corr(c_idx);
%             end
%             kdtree_delete(tree);
%         end
%         tree=kdtree_build(em_pack.D_mat(end-2:end,:)');
%         if numel(c_idx)>3
%             [T,~,~,~] = icp(X(1:3,c_idx), neighbors(end,c_idx),em_pack.D_mat(end-2:end,:), tree, 0.5);
%         else
%             T=eye(4);
%         end
%         X_tmp=T*[X(1:3,:);ones(1,N)];
%         [corr,~]=kdtree_nearest_neighbor(tree,X_tmp(1:3,:)');
%         em_pack.X_mat(end-3:end-1,c_idx)=em_pack.D_mat(end-2:end,corr(c_idx));
%         X(end-3:end-1,c_idx)=em_pack.D_mat(end-2:end,corr(c_idx));
%         if numel(c_idx)>3
%             neighbors(end,c_idx)=corr(c_idx);
%             em_pack.Neighbors(end,c_idx)=corr(c_idx);
%         end
%         kdtree_delete(tree);
%     end
    if sum(em_pack.is_new)==0
%         for g_idx=1:group_num
%             [c_idx,~]=find(clusters==g_idx);
%             for iter=1:m-1
%                 tree=kdtree_build(em_pack.D_mat(3*(iter-1)+1:3*iter,:)');
%                 if numel(c_idx)>3
%                     [T,~,~,~] = icp(em_pack.X_mat(4*(iter-1)+1:4*(iter-1)+3,c_idx), em_pack.Neighbors(iter,c_idx),em_pack.D_mat(3*(iter-1)+1:3*iter,:), tree, 0.5);
%                 else
%                     T=eye(4);
%                 end
%                 X_tmp=T*[em_pack.X_mat(4*(iter-1)+1:4*(iter-1)+3,:);ones(1,N)];
%                 [corr,~]=kdtree_nearest_neighbor(tree,X_tmp(1:3,:)');
%                 em_pack.X_mat(4*(iter)+1:4*(iter)+3,c_idx)=em_pack.D_mat(3*(iter-1)+1:3*iter,corr(c_idx));
%                 if numel(c_idx)>3
%                     em_pack.Neighbors(iter,c_idx)=corr(c_idx);
%                 end
%                 kdtree_delete(tree);
%             end
%             tree=kdtree_build(em_pack.D_mat(end-2:end,:)');
%             if numel(c_idx)>3
%                 [T,~,~,~] = icp(X(1:3,c_idx), neighbors(end,c_idx),em_pack.D_mat(end-2:end,:), tree, 0.5);
%             else
%                 T=eye(4);
%             end
%             X_tmp=T*[X(1:3,:);ones(1,N)];
%             [corr,~]=kdtree_nearest_neighbor(tree,X_tmp(1:3,:)');
%             em_pack.X_mat(end-3:end-1,c_idx)=em_pack.D_mat(end-2:end,corr(c_idx));
%             X(end-3:end-1,c_idx)=em_pack.D_mat(end-2:end,corr(c_idx));
%             if numel(c_idx)>3
%                 neighbors(end,c_idx)=corr(c_idx);
%                 em_pack.Neighbors(end,c_idx)=corr(c_idx);
%             end
%             kdtree_delete(tree);
%         end

        em_pack.X_mat=em_pack.X_mat(1:end-4,:);
        em_pack.D_mat=em_pack.D_mat(1:end-3,:);
        em_pack.Neighbors=em_pack.Neighbors(1:end-1,:);
    else
%         for g_idx=1:group_num
%             [c_idx,~]=find(clusters==g_idx);
%             tree=kdtree_build(em_pack.D_mat(end-2:end,:)');
%             [tmp_corr,~]=kdtree_nearest_neighbor(tree,em_pack.X_mat(end-7:end-5,c_idx)');
%             if numel(c_idx)>3
%                 [T,~,~,~] = icp(X(1:3,c_idx), tmp_corr,em_pack.D_mat(end-2:end,:), tree, 0.5);
%             else
%                 T=eye(4);
%             end
%             X_tmp=T*[X(1:3,:);ones(1,N)];
%             [corr,~]=kdtree_nearest_neighbor(tree,X_tmp(1:3,:)');
%             em_pack.X_mat(end-3:end-1,c_idx)=em_pack.D_mat(end-2:end,corr(c_idx));
%             X(end-3:end-1,c_idx)=em_pack.D_mat(end-2:end,corr(c_idx));
%             if numel(c_idx)>3
%                 em_pack.Neighbors(end,c_idx)=corr(c_idx);
%                 neighbors(end,c_idx)=corr(c_idx);
%             end
%             kdtree_delete(tree);
%         end
    end

end