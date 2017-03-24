function [em_pack]=EM_estep(em_pack,X_candi)
    [~,N]=size(em_pack.X_mat);
%     nn=M/4-1;
    group_num=numel(em_pack.rank_q);
    k=4;
%     rank_c=zeros(group_num,1);
%     X_norm=X-repmat(X(1:3,:),M/3,1);

    
    for iter=1:N
        em_pack.R(iter,:)=0;
        for iter2=1:group_num
            
            C_inv=em_pack.covinv{iter2};
            U=em_pack.W_mat{iter2}(:,em_pack.rank_q(iter2)+1:end);
            if nargin == 1
%                 xx=U'*normc(X(:,iter));
                xx=U'*normc(em_pack.X_mat(:,iter));
            else
%                 xx=U'*normc([X(1:end-k,iter);X_candi(3*(iter2-1)+1:3*iter2,iter);1]);
                xx=U'*normc([em_pack.X_mat(1:end-k,iter);X_candi(3*(iter2-1)+1:3*iter2,iter);1]);
            end
            em_pack.Mdist(iter,iter2)=(xx'*C_inv*xx);%+(norm(x))/length(x);%./sqrt(length(C_inv));%./((length(C_inv)/3)^0.5);%3*(xx'*C_inv*xx)./length(C_inv);%*em_pack.P_c(iter2)
%             em_pack.Mdist(iter,iter2)=(xx'*xx)*norm(C_inv);
            MM=em_pack.Mdist(iter,iter2);
            if MM<0
                em_pack.R(iter,iter2)=0;
            else
                C=em_pack.cov{iter2};
                em_pack.R(iter,iter2)=((2*pi)^(-((size(C,1))/2))/sqrt(norm(C)))*exp(-0.5*MM)*em_pack.P_c(iter2);%
%                 em_pack.R(iter,iter2)=((2*pi)^(-((size(C,1))/2))/sqrt(norm(C)))*exp(-0.5*MM)*em_pack.P_c(iter2);%
            end
  
            if em_pack.R(iter,iter2)==Inf
                em_pack.R(iter,iter2)=10^(100);
            end
        end

        if sum(em_pack.R(iter,:))==0
%             [~,m_clust]=min(em_pack.Mdist(iter,:));
%             em_pack.R(iter,m_clust)=0.1;%ones(1,group_num);
            em_pack.R(iter,:)=ones(1,group_num);
        end
        em_pack.R(iter,:)=em_pack.R(iter,:)/sum(em_pack.R(iter,:));
    end
    
        
end