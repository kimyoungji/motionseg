function [X,neighbors,em_pack,param_pack,sum_seq]=xupdate_whole(X,Data_mat,neighbors,clusters,group_num,em_pack,sum_seq)
    [M,N]=size(em_pack.X_mat);
    for g_idx=1:group_num
        em_pack.T_pre{g_idx}=em_pack.T_cur{g_idx}*em_pack.T_pre{g_idx};
        em_pack.T_cur{g_idx}=eye(4);
    end

%% update motions according to motion signatures
    m=M/4;
    param_pack=cell(1,group_num);
    if m>16
        if sum_seq(2)==0
        M_mat={};
        for g_idx=2:group_num
            [c_idx,~]=find(clusters==g_idx);
            M_mat{g_idx}=zeros(m-1,12);
            for iter=2:m-1
                %icp
                tree=kdtree_build(em_pack.D_mat(3*(iter-1)+1:3*iter,:)');
                if numel(c_idx)>3
                    [T,~,~,~] = icp(em_pack.X_mat(5:7,c_idx), em_pack.Neighbors(iter,c_idx),em_pack.D_mat(3*(iter-1)+1:3*iter,:), tree, 0.5);
                else
                    T=eye(4);
                end
                X_tmp=T*[em_pack.X_mat(5:7,:);ones(1,N)];
                [corres_1,~]=kdtree_nearest_neighbor(tree,X_tmp(1:3,:)');
%                 em_pack.X_mat(4*(iter)+1:4*(iter)+3,c_idx)=em_pack.D_mat(3*(iter-1)+1:3*iter,corres_1(c_idx));
                if numel(c_idx)>3
                    em_pack.Neighbors(iter,c_idx)=corres_1(c_idx);
                end
                kdtree_delete(tree);
                
                R_part=T(1:3,1:3)-eye(3);
                t_part=T(1:3,4);
                M_mat{g_idx}(iter,:)=[R_part(:)',t_part'];
            end
            [simple_M,param]=MotionSignatures(M_mat{g_idx});
            
            if param.signature(1)<9 && param.signature(2)<3
                if param.signature(1)==0 && param.signature(2)==0
                else
                    param_pack{g_idx}=param;
%                     for iter=1:m-1
%                         %icp
%                         tree=kdtree_build(em_pack.D_mat(3*(iter-1)+1:3*iter,:)');
%                         X_pre=simple_M(4*(iter-1)+1:4*iter,:)*[em_pack.X_mat(1:3,:);ones(1,N)];
%                         if numel(c_idx)>3
%                             [T,~,~,~] = icp(X_pre(1:3,c_idx), em_pack.Neighbors(iter,c_idx),em_pack.D_mat(3*(iter-1)+1:3*iter,:), tree, 0.5);
%                         else
%                             T=eye(4);
%                         end
%                         X_tmp=T*[X_pre(1:3,:);ones(1,N)];
%                         [corres_1,~]=kdtree_nearest_neighbor(tree,X_tmp(1:3,:)');
%                         em_pack.X_mat(4*(iter)+1:4*(iter)+3,c_idx)=em_pack.D_mat(3*(iter-1)+1:3*iter,corres_1(c_idx));
%                         if numel(c_idx)>3
%                             em_pack.Neighbors(iter,c_idx)=corres_1(c_idx);
%                         end
%                         kdtree_delete(tree);
%                     end
                end
            end
        end
        end
        sum_seq(2)=sum_seq(2)+1;
        if sum_seq(2)==10
            sum_seq(2)=0;
        end
    end
%% update X_mat
        if sum(em_pack.is_new)==0 && sum_seq(1)<2
            em_pack.X_mat=em_pack.X_mat(1:end-4,:);
            em_pack.D_mat=em_pack.D_mat(1:end-3,:);
            em_pack.Neighbors=em_pack.Neighbors(1:end-1,:);
            sum_seq(1)=sum_seq(1)+1;
        else
            sum_seq(1)=0;
        end
%         if sum(em_pack.is_new)==0
%             em_pack.X_mat=em_pack.X_mat(1:end-4,:);
%             em_pack.D_mat=em_pack.D_mat(1:end-3,:);
%             em_pack.Neighbors=em_pack.Neighbors(1:end-1,:);
%         end
end