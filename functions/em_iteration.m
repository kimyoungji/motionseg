function [X,neighbors,em_pack,last_clusters,x_mesh]=em_iteration(X,dataMatrix,neighbors,em_pack,x_mesh)
%% initialize
    end_e=3;
    ee=0;
    is_converge=0;
    X_candi=repmat(X(end-3:end-1,:),numel(em_pack.rank_q),1);
    
    while ~is_converge
        group_num=size(em_pack.R,2);

       %% clustering with graphcut
        if group_num>1
            
            gc_dist=(ones(size(em_pack.R))-(em_pack.R));%-log(em_pack.R)*10^(5);%%exp(-em_pack.R*100);%
            smoothcost=(ones(group_num)-eye(group_num));

            M_tmp=normc(em_pack.X_mat)';%em_pack.Mdist;
            [em_pack.clusters,energy]=Gcut_update(gc_dist,smoothcost,x_mesh,em_pack.clusters,M_tmp,X(3,:));
            
            m_energy=max(energy(:));
            em_pack.R=(m_energy*ones(size(energy))-energy);%*(10^5);
            em_pack.R=em_pack.R./repmat(sum(em_pack.R,2),1,group_num);
            
            m_clusters=em_pack.clusters;
            group_num=size(energy,2);
            for iter=1:group_num
                c_idx=find(em_pack.clusters==iter);
%                 numel(c_idx)
                if numel(c_idx)>10
                    [idx2,idx1]=deleteoutliers(energy(c_idx,iter)./(10^3),1);
                    if numel(idx1)>10 && numel(idx2)>10
                        if mean(energy(c_idx(idx1),iter))<mean(energy(c_idx(idx2),iter))
                            idx1=idx2;
                        end
                        m_clusters(c_idx(idx1))=0;
                    end
                end
            end
            
        else
            em_pack.clusters=ones(size(em_pack.clusters));
            em_pack.R=ones(em_pack.N,1);
            m_clusters=em_pack.clusters;
        end

        

       %% M-Step
        [check_bin,bin_contents,em_pack,group_num]=EM_mstep(em_pack,m_clusters);

        % find null cluster
        em_pack.e_criteria
        check_bin
        em_pack=reject_nc(em_pack,check_bin,bin_contents);
        if group_num==1
            em_pack.R=ones(em_pack.N,1);
        end
       
        %str=sprintf('test%d_%d_1.ply',s_num,ee+1);
        %save_ply(str,em_pack.clusters,X(end-3:end-1,:));

        
       %% X update
       %if ee==1
        [X,neighbors,em_pack,X_candi,x_mesh,occ_idx]=EM_xupdate(em_pack,X,neighbors,dataMatrix,x_mesh);% X update
        em_pack.N=size(X,2);
       %end
       %% E-Step
        em_pack=EM_estep(em_pack,X_candi);%

       %% check convergence
        em_pack.e_criteria
        ee=ee+1;
        if ee>=end_e
            is_converge=1;
        end
        last_clusters=em_pack.clusters;
       %% outlier handling
        %find o_clusters
%         o_idx=[];%
        o_idx=find(occ_idx==1);
        tmp_Mdist=zeros(em_pack.N,1);
        for iter=1:group_num
            c_check=find(em_pack.clusters==iter);
            tmp_Mdist(c_check)=em_pack.Mdist(c_check,iter);
        end
        [idx2,idx1]=deleteoutliers(tmp_Mdist,10);
        if mean(tmp_Mdist(idx1))<mean(tmp_Mdist(idx2))
            idx1=idx2;
        end
        out_idx=idx1;
        o_idx=[o_idx;out_idx];
        em_pack.clusters(o_idx)=0;
        
        if is_converge<1
        %find new clusters
        if numel(o_idx)>5
            [W,rr,inliers,C,dists,em_pack]=new_subspace(o_idx,dataMatrix(end-2:end,:),X,neighbors(end,o_idx),em_pack);
            d=1000;
            min_disp=iter;
            big_c=0;
            c_sum=0;
            for iter3=1:group_num
                W1=em_pack.W_mat{iter3};
                rr1=em_pack.rank_q(iter3);
                tmp_d=subspace_disparity(W1(:,1:rr1),W(:,1:rr));
                if tmp_d<d
                    d=tmp_d;
                    min_disp=iter3;
                end
                if sum(last_clusters(inliers)==iter3)>c_sum
                    c_sum=sum(last_clusters(inliers==iter3));
                    big_c=iter3;
                end
                
            end
            
            str=sprintf('%d: %f, %f, %f',iter, norm(C), d, dists )%
            if  group_num<10 && d>0.5 %&& norm(C)<10^(-2) && dists>10^(-3)
                   %create a new cluster
                    group_num=group_num+1;
                    em_pack.R(:,group_num)=zeros(em_pack.N,1);%ones(em_pack.N,1)*0.01;%
                    em_pack.R(inliers,:)=0;
                    em_pack.R(inliers,group_num)=1;%ones(em_pack.N,1);
%                     em_pack.R=em_pack.R./repmat(sum(em_pack.R,2),1,group_num);

                    em_pack.Mdist(:,group_num)=zeros(em_pack.N,1);
                    em_pack.W_mat{group_num}=W;
                    em_pack.T_pre{group_num}=em_pack.T_pre{big_c};
                    em_pack.T_cur{group_num}=eye(4);
                    em_pack.rank_q(group_num)=rr;
                    em_pack.cov{group_num}=C;%+eye(size(C))*(10^-10);
                    em_pack.covinv{group_num}=pinv(em_pack.cov{group_num});
                    em_pack.e_criteria(group_num)=norm(em_pack.cov{group_num});
                    em_pack.P_c(group_num)=numel(inliers)/(em_pack.N);%0.1;%1/group_num;%+numel(out_idx)/em_pack.N;%A_d;%/em_pack.N;

                    em_pack.dists{group_num}= em_pack.dists{1}*0;
                    em_pack.clusters(inliers)=group_num;
                    em_pack.is_new(group_num)=1;
            else
                em_pack=EM_estep(em_pack);%,X_candi
            end
        else
            em_pack=EM_estep(em_pack);%,X_candi
        end
        end

        
       %% save data for clusters
        %str=sprintf('test%d_%d_2.ply',s_num,ee);
        %save_ply(str,em_pack.clusters,X(end-3:end-1,:));

    end


    ee
    em_pack.e_criteria
end
