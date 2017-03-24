function [X,neighbors,em_pack,occ_idx,T_mat]=EM_xupdate(em_pack,X,neighbors,dataMatrix)

    group_num=numel(em_pack.rank_q);
    k=3;
    [~,N]=size(X);
    n=size(dataMatrix,1)/k;
    in_idx=3;
    p_idx=n;
    
    ransacCoef.minPtNum=3;
    ransacCoef.iterNum=10;
    ransacCoef.thInlrRatio=0.5;
    ransacCoef.thDist=0.5;
                
    tree=kdtree_build(dataMatrix(k*(p_idx-1)+1:k*p_idx,:)');%create kd tree
    T_mat=cell(1,group_num);
    occ_idx=zeros(N,1);
    for g_idx=1:group_num
        
        clust_idx=find(em_pack.clusters==g_idx);
        sampled_N=length(clust_idx);
        
        if length(clust_idx)>3 

                U=em_pack.W_mat{g_idx}(:,em_pack.rank_q(g_idx)+1:end);
                C_inv=em_pack.covinv{g_idx};
                matching=zeros(1,sampled_N);
                for iter=1:sampled_N
                    [candidates,~]=kdtree_k_nearest_neighbors(tree,X(end-5:end-3,clust_idx(iter))',in_idx);
                    dist=zeros(in_idx,1);
                    for iter2=1:in_idx
                        x=U'*[X(1:end-3,clust_idx(iter));dataMatrix(end-2:end,candidates(iter2))];
                        dist(iter2)=(x'*C_inv*x)./sqrt(length(C_inv));%
                    end
                    [~,m_idx]=min(dist);
                    matching(iter)=candidates(m_idx);
                end

                %correspondence propagation
                no_idx=find(em_pack.clusters==g_idx);
                length(no_idx)
                % update neighbors
                x_tmp=X(k*(p_idx-1)+1:k*p_idx,clust_idx);
                corr=matching;

                T_mat{g_idx}=eye(4);
                for iter=1:10
                    [T, ~, ~, corr] = ransac2( x_tmp,dataMatrix(k*(p_idx-1)+1:k*p_idx,corr),tree,dataMatrix(k*(p_idx-1)+1:k*p_idx,:),ransacCoef,@estimateRigidTransform,@estimateDistance);
                    T_mat{g_idx}=T*T_mat{g_idx};
                    % update X
                    x_tmp=T*[x_tmp;ones(1,length(clust_idx))];
                    x_tmp=x_tmp(1:3,:);
                end
                temp_X=T_mat{g_idx}*[X(k*(p_idx-1)+1:k*p_idx,:);ones(1,N)];
                
                % save X_nn
                [corr,dist]=kdtree_nearest_neighbor(tree,temp_X(1:3,:)');
                neighbors(p_idx,no_idx)=corr(no_idx);
                em_pack.dists(g_idx)=mean(dist(no_idx));

        end

    end

    %determine
    X(k*p_idx+1:k*p_idx+3,:)=dataMatrix(k*(p_idx-1)+1:k*p_idx,neighbors(p_idx,:));
    
    %find occluded region
    for iter=1:N
        n_idx=neighbors(end,iter);
        if sum(neighbors(end,:)==n_idx)>1
            occ_idx(iter)=1;
        end
    end

    kdtree_delete(tree);
end