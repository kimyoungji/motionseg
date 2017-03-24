function [check_bin,bin_contetns,em_pack,g_num]=EM_mstep(em_pack,m_clusters)
    [~,group_num]=size(em_pack.R);
%     [M,N]=size(X);
    [M,N]=size(em_pack.X_mat);
    nn=M/4-1;
    
%     X_norm=normc(X);
    X_norm=normc(em_pack.X_mat);
    check_bin=zeros(1,group_num);

    
    for iter=1:group_num
        [out_idx,~]=find(m_clusters==iter);
        R_tmp=zeros(N,1);%em_pack.R(:,iter);%*10^(-3);%
        R_tmp(out_idx)=1;%em_pack.R(out_idx,iter);
        R_tmp2=zeros(N,1);%em_pack.R(:,iter);%*10^(-3);%
        R_tmp2(out_idx)=em_pack.R(out_idx,iter);
        
        if numel(out_idx)>10 %numel(out_idx)>30 %
 
            check_bin(iter)=iter;
            % construct A
            A_a=zeros(M,M);A_d=0;
            for iter2=1:N
                A_a=A_a+R_tmp2(iter2)*X_norm(:,iter2)*X_norm(:,iter2)';
                A_d=A_d+R_tmp2(iter2);
            end
            
            % find W
            [~,D,W]=svd(A_a);
            s_values=diag(D);

            em_pack.rank_q(iter)=find_rank(s_values(1:end-nn));
            if em_pack.rank_q(iter)>4
                em_pack.rank_q(iter)=4;
            end
            if em_pack.rank_q(iter)<=1 || em_pack.rank_q(iter)>=M-1
                em_pack.rank_q(iter)=4;
            end

            
            %assign parameters
            em_pack.W_mat{iter}=W;
            S=zeros(M-em_pack.rank_q(iter),M-em_pack.rank_q(iter));
            for iter2=1:N
                S=S+R_tmp2(iter2)*(em_pack.W_mat{iter}(:,em_pack.rank_q(iter)+1:end)'*X_norm(:,iter2))*(em_pack.W_mat{iter}(:,em_pack.rank_q(iter)+1:end)'*X_norm(:,iter2))';
            end
            C=S/A_d;
            
%             C=C+eye(M-em_pack.rank_q(iter))*(10^-10);
            em_pack.cov{iter}=C;
            em_pack.e_criteria(iter)=norm(C);
            em_pack.P_c(iter)=sum(R_tmp2);
            em_pack.covinv{iter}=pinv(em_pack.cov{iter});
    

            if em_pack.e_criteria(iter)==0 %|| em_pack.e_criteria(iter)<10^-10 %|| em_pack.e_criteria(iter)>1
                check_bin(iter)=0;
            end

        else
            check_bin(iter)=0;
        end
    end
    em_pack.P_c=em_pack.P_c./sum(em_pack.P_c);


    %merge clusters
    for iter=1:group_num
        m_dist=100;
        m_ids=iter;
        for iter2=1:group_num
            if iter~=iter2 && check_bin(iter)>0 && check_bin(iter2)>0 
                d=subspace_disparity(em_pack.W_mat{iter}(:,1:em_pack.rank_q(iter)),em_pack.W_mat{iter2}(:,1:em_pack.rank_q(iter2)));
                if m_dist>d
                    m_dist=d;
                    m_ids=iter2;
                end
            end
        end
        if m_dist<0.1 %&& m_ids~=iter
            if em_pack.e_criteria(iter)>em_pack.e_criteria(m_ids)
                check_bin(iter)=m_ids;
            else
                check_bin(m_ids)=iter;
            end
        end
    end
    
    g_num=0;
    bin_contetns=zeros(group_num,1);
    for iter=1:group_num
        c_idx=find(check_bin==iter);
        if numel(c_idx)>0
            g_num=g_num+1;
            bin_contetns(c_idx)=g_num;
        end
    end
    
end