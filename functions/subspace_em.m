function [res_clusters,res_neighbors,res_dist]=subspace_em(X_ori,group_num,x1_mesh)
    % x: input data k x N
    % y: data matrix (k*n) x N
    % init_clusters: initial clusters
    % group_num: number of groups
    %[k,N]=size(x);

    k=4;
    N=size(X_ori,2);
    n=300;
    eltime=zeros(n,1);
    
    %% initialize
    
    em_pack.rank_q=ones(1,group_num);
    em_pack.dists={};%zeros(1,group_num);
    em_pack.W_mat={};
    em_pack.T_pre={};
    em_pack.T_cur={};
    em_pack.X_mat=[];
    em_pack.D_mat=[];
    em_pack.Neighbors=[];
    em_pack.cov={};
    em_pack.P_c=ones(1,group_num);
    em_pack.R=zeros(N,group_num);
    em_pack.Mdist=zeros(N,group_num);
    em_pack.e_criteria=zeros(group_num,1);
    em_pack.is_new=zeros(1,group_num);
    em_pack.N=N;
    param_pack={};
 
    
    %% EM
    neighbors=[];%init_neighbors;
    Data_mat=[];
    X=[X_ori(1:3,:);ones(1,em_pack.N)];
    em_pack.X_mat=[X_ori(1:3,:);ones(1,em_pack.N)];
    m_seq_num=5;
    sum_seq=zeros(2,1);
    em_pack.dists{1}=0;
    em_pack.T_pre{1}=eye(4);
    em_pack.T_cur{1}=eye(4);
    param_pack{1}=[];
    for iter=1:n
        tic;
        %arrange X_mat, Data_mat and neighbors
        x_for_neighbors=X(end-3:end-1,:);
        if iter>m_seq_num%%mod(iter,m_seq_num)==0
            X=[X(1:k,:);X(k*2+1:end,:)];%X_mat(5:end,:);%
            %X_mat=X_mat(5:end,:);
            Data_mat=Data_mat(4:end,:);
            neighbors=neighbors(2:end,:);
            
            for gg=1:group_num
                em_pack.dists{gg}(1:end-1)=em_pack.dists{gg}(2:end);
                em_pack.dists{gg}(end)=0;
            end
        end
        
        if (size(em_pack.X_mat,1)/k)==17
            em_pack.X_mat=[em_pack.X_mat(1:k,:);em_pack.X_mat(k*2+1:end,:)];
            em_pack.D_mat=em_pack.D_mat(4:end,:);
            em_pack.Neighbors=em_pack.Neighbors(2:end,:);
        end

        %load data
        % clouds from img
        ff=sprintf('./bonn/watercan/%d.png',200+1*iter);
%         ff=sprintf('./rigid_dataset/car_d1/%d.png',60+1*iter);
        IMG=imread(ff);
        x=img_to_xyz(IMG,1,1);%./100;


        if ~isempty(Data_mat)
            if size(x,1)>size(Data_mat,2)
                Data_mat=[Data_mat,zeros(size(Data_mat,1),size(x,1)-size(Data_mat,2))];
                em_pack.D_mat=[em_pack.D_mat,zeros(size(em_pack.D_mat,1),size(x,1)-size(em_pack.D_mat,2))];
            elseif size(x,1)<size(Data_mat,2)
                x=[x;zeros(size(Data_mat,2)-size(x,1),k-1)];
            end
        end

        %update nearest neighbors
        Data_mat=[Data_mat;x'];
        em_pack.D_mat=[em_pack.D_mat;x'];
        tree=kdtree_build(x);%create kd tree
        tmp_neighbors=[];
        [tmp_neighbors,dist]=kdtree_nearest_neighbor(tree,x_for_neighbors');
        kdtree_delete(tree);
        neighbors=[neighbors;tmp_neighbors'];
        em_pack.Neighbors=[em_pack.Neighbors;tmp_neighbors'];
        X=[X;Data_mat(end-2:end,tmp_neighbors);ones(1,em_pack.N)];
        em_pack.X_mat=[em_pack.X_mat;Data_mat(end-2:end,tmp_neighbors);ones(1,em_pack.N)];
        
        if iter==1
            em_pack.clusters=ones(em_pack.N,1);
            em_pack.R=zeros(em_pack.N,group_num);%/(group_num+1);
            tree=kdtree_build(X(end-2*k+1:end-k-1,:)');
            for iter3=1:em_pack.N
                for iter2=1:group_num
                    if(em_pack.clusters(iter3)==iter2)
                        em_pack.R(iter3,iter2)=1;
                    end
                end
            end
            kdtree_delete(tree);

        end
        
        em_pack.is_new=zeros(1,group_num);
        [X,neighbors,em_pack,last_clusters,x1_mesh]=em_iteration(X,Data_mat,neighbors,em_pack,x1_mesh);
        %group_num=max(tmp_clusters);%size(em_pack.R,2);
        group_num=max(last_clusters);%max(em_pack.clusters);
        str=sprintf('exp_%d.ply',iter);
        save_ply(str,last_clusters,X(end-k+1:end,:));

%         [X,neighbors,em_pack]=xupdate_whole2(X,neighbors,last_clusters,group_num,em_pack);
        [X,neighbors,em_pack,param_pack,sum_seq]=xupdate_whole(X,Data_mat,neighbors,last_clusters,group_num,em_pack,sum_seq);
        % find ax_points
        
        for g_idx=1:group_num
            param=param_pack{g_idx};
            if ~isempty(param)
                ax_points=zeros(3,400);
                if param.signature(1)>0
                    ax_points(1:3,1:100)=repmat(param.t_a,1,100)+param.a*[-10:0.2:9.8];
                    ax_points(1:3,101:200)=repmat(param.t_b,1,100)+param.b*[-10:0.2:9.8];
                end
                if param.signature(2)>0
                    for t_idx=1:param.signature(2)
                        ax_points(1:3,101+t_idx*100:200+t_idx*100)=param.T(:,t_idx)*[-10:0.2:9.8];
                    end
                end
                str=sprintf('ax_%d_%d.ply',iter,g_idx);
                save_ax_ply(str,ax_points);
            end
        end
        eltime(iter)=toc;
    end
    save('eltime.mat','eltime');
    res_clusters=last_clusters;
    res_neighbors=neighbors(end,:);
    res_dist=zeros(1,em_pack.N);
    for iter=1:em_pack.N
        d=(Data_mat(end-2:end,res_neighbors(iter))-Data_mat(end-2:end,iter));
        res_dist(iter)=(d'*d);
    end
    plot_corres(X,res_clusters);

end
