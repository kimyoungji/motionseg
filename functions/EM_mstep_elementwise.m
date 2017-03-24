function [package]=EM_mstep_elementwise(X, package,p_idx,out_idx)
    [M,N]=size(X);
    [U,D,V]=svd(X(:,out_idx));
    package.rank_q(p_idx)=find_rank(diag(D));
    if package.rank_q(p_idx)>4
        package.rank_q(p_idx)=4;
    end

    package.W_mat{p_idx}=U;
    package.p_clusters(out_idx)=p_idx;

    %mean update
%     R_tmp=package.R(:,p_idx);
    R_tmp=zeros(N,1);
    R_tmp(out_idx)=1;
    U=U(:,package.rank_q(p_idx)+1:end);
    package.m_value{p_idx}=(R_tmp'*(U'*X)')'/sum(R_tmp);
    %covariance update
    S=zeros(M-package.rank_q(p_idx),M-package.rank_q(p_idx));
    for iter2=1:N
        S=S+R_tmp(iter2)*(U'*X(:,iter2)-package.m_value{p_idx})*(U'*X(:,iter2)-package.m_value{p_idx})';
    end
    package.cov{p_idx}=S/(sum(R_tmp));
    package.cov{p_idx}=package.cov{p_idx}+eye(size(package.cov{p_idx}))*10^(-300);
    C=package.cov{p_idx};
    package.e_criteria(p_idx)=norm(C)/length(C);
    package.R(:,p_idx)=R_tmp;
    package.P_c(p_idx)=sum(R_tmp)/N;
    package.Mdist(:,p_idx)=R_tmp;
end