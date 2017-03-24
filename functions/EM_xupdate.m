function [X,neighbors,em_pack,X_candi,x_mesh,occ_idx]=EM_xupdate(em_pack,X,neighbors,dataMatrix,x_mesh)

    group_num=numel(em_pack.rank_q);
    k=4;
    [~,N]=size(X);
    n=size(dataMatrix,1)/3;
    p_idx=n;
    
    ransacCoef.minPtNum=4;
    ransacCoef.iterNum=10;
    ransacCoef.thInlrRatio=0.5;
    ransacCoef.thDist=1;
                
    tree=kdtree_build(dataMatrix(3*(p_idx-1)+1:3*p_idx,:)');%create kd tree
    T_mat=cell(1,group_num);
    occ_idx=zeros(N,1);
    X_candi=repmat(X(end-3:end-1,:),group_num,1);
    dists=zeros(N,1);
    for g_idx=1:group_num
        
        clust_idx=find(em_pack.clusters==g_idx);
        
        if numel(clust_idx)>3
                if em_pack.is_new(g_idx)==1
                    [corres_1,~]=kdtree_nearest_neighbor(tree,X(end-2*k+1:end-k-1,clust_idx)');
%                     [corres_1,~]=kdtree_nearest_neighbor(tree,em_pack.X_mat(end-7:end-5,clust_idx)');
                    m_iter=20;
                else
                    corres_1=neighbors(end,clust_idx);
                    m_iter=10;
                end
                length(clust_idx)
                % update neighbors
                %x_tmp=X(1:3,clust_idx);
                x_t=em_pack.T_pre{g_idx}*X(1:4,clust_idx);
                x_tmp=x_t(1:3,:);
                
                T_mat{g_idx}=eye(4);
                for iter=1:m_iter
                    [T, ii, oo, corres_1] = ransac2(x_tmp,dataMatrix(3*(p_idx-1)+1:3*p_idx,corres_1),tree,dataMatrix(3*(p_idx-1)+1:3*p_idx,:),ransacCoef,@estimateRigidTransform,@estimateDistance);
%                     if em_pack.is_new(g_idx)==1
%                         T_mat{g_idx}=T*T_mat{g_idx};
%                         % update X
%                         x_tmp=T*[x_tmp;ones(1,length(clust_idx))];
%                         x_tmp=x_tmp(1:3,:);
%                         if icp_convergence(T)<0.005
%                             break;
%                         end
%                     else
%                         if numel(ii)>numel(oo)
%                             T_mat{g_idx}=T*T_mat{g_idx};
%                             % update X
%                             x_tmp=T*[x_tmp;ones(1,length(clust_idx))];
%                             x_tmp=x_tmp(1:3,:);
%                             if icp_convergence(T)<0.005
%                                 break;
%                             end
%                         end 
%                     end
%                     if numel(ii)>numel(oo)
                    T_mat{g_idx}=T*T_mat{g_idx};
                    % update X
                    x_tmp=T*[x_tmp;ones(1,length(clust_idx))];
                    x_tmp=x_tmp(1:3,:);
                    if icp_convergence(T)<0.0001
                        break;
                    end
%                     end
                end
                temp_X=T_mat{g_idx}*em_pack.T_pre{g_idx}*[X(1:3,:);ones(1,N)];
                
                % save X_nn
                [corres_1,dist_1]=kdtree_nearest_neighbor(tree,temp_X(1:3,:)');
                neighbors(p_idx,clust_idx)=corres_1(clust_idx);
                em_pack.Neighbors(end,clust_idx)=corres_1(clust_idx);
                dists(clust_idx)=dist_1(clust_idx);
                em_pack.dists{g_idx}(p_idx)=mean(dist_1(clust_idx));
                X_candi(3*(g_idx-1)+1:3*g_idx,:)=dataMatrix(3*(p_idx-1)+1:3*p_idx,corres_1);
        else
            T_mat{g_idx}=eye(4);
        end

    end

    %determine
    X(k*p_idx+1:k*p_idx+3,:)=dataMatrix(3*(p_idx-1)+1:3*p_idx,neighbors(p_idx,:));
    em_pack.X_mat(end-3:end-1,:)=dataMatrix(3*(p_idx-1)+1:3*p_idx,neighbors(p_idx,:));
    em_pack.T_cur=T_mat;
    
%     %find occluded regions or outliers
%     for iter=1:N
%         n_idx=neighbors(end,iter);
%         if dists(iter)>0.3 %%sum(neighbors(end,:)==n_idx)>2 %|| 
%             occ_idx(iter)=1;
%         end
%     end

%     %find occluded regions or outliers
%     prev_max=size(X,2);
%     c_idx=find(dists<=0.3);
%     X=X(:,c_idx);
%     neighbors=neighbors(:,c_idx);
%     X_candi=X_candi(:,c_idx);
%     em_pack.R=em_pack.R(c_idx,:);
%     em_pack.Mdist=em_pack.Mdist(c_idx,:);
%     em_pack.clusters=em_pack.clusters(c_idx,:);
%     
%     %recalculate x_mesh
%     c_idx=find(dists>0.3);
%     x_mesh=rep_x_mesh(x_mesh,c_idx,prev_max);
    
    
    kdtree_delete(tree);
end