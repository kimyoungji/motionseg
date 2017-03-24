function res_pack=reject_nc(in_pack,check_bin,bin_contents)
    gg=max(bin_contents);
    [N,~]=size(in_pack.R);
    
    res_pack.rank_q=ones(1,gg);
    res_pack.dists={};
    res_pack.W_mat={};
    res_pack.T_pre={};
    res_pack.T_cur={};
    res_pack.X_mat=[];
    res_pack.D_mat=[];
    res_pack.Neighbors=[];
    res_pack.cov={};
    res_pack.covinv={};
    res_pack.P_c=zeros(1,gg);
    res_pack.e_criteria=ones(gg,1);
    res_pack.R=zeros(N,gg);
    res_pack.Mdist=zeros(N,gg);
    res_pack.clusters=zeros(size(in_pack.clusters));
    res_pack.is_new=zeros(1,gg);
    res_pack.N=0;
    
    for iter=1:gg
        tg_idx=find(bin_contents==iter);
        src_idx=check_bin(tg_idx(1));
        res_pack.rank_q(iter)=in_pack.rank_q(src_idx);
        res_pack.dists{iter}=in_pack.dists{src_idx};
        res_pack.W_mat{iter}=in_pack.W_mat{src_idx};
        res_pack.T_pre{iter}=in_pack.T_pre{src_idx};
        res_pack.T_cur{iter}=in_pack.T_cur{src_idx};
        res_pack.cov{iter}=in_pack.cov{src_idx};
        res_pack.covinv{iter}=in_pack.covinv{src_idx};
        res_pack.P_c(iter)=in_pack.P_c(src_idx);
        res_pack.e_criteria(iter)=in_pack.e_criteria(src_idx);
        res_pack.R(:,iter)=in_pack.R(:,src_idx);
        res_pack.Mdist(:,iter)=in_pack.Mdist(:,src_idx);
        res_pack.is_new(iter)=in_pack.is_new(src_idx);
        res_pack.N=in_pack.N;
        
        for iter2=1:numel(tg_idx)    
            res_pack.clusters(in_pack.clusters==tg_idx(iter2))=iter;
        end
    end
    
    zero_idx=find(res_pack.clusters==0);
    if gg==1
        res_pack.clusters(zero_idx)=1;
    else
        [~,m_idx]=max(res_pack.R(zero_idx,:)');
        res_pack.clusters(zero_idx)=m_idx;
    end
    res_pack.P_c=(res_pack.P_c)./sum(res_pack.P_c);
    check_R=sum(res_pack.R');
    res_pack.R(check_R==0,:)=1/gg;
    res_pack.R=(res_pack.R)./repmat(sum(res_pack.R')',1,gg);
    res_pack.X_mat=in_pack.X_mat;
    res_pack.D_mat=in_pack.D_mat;
    res_pack.Neighbors=in_pack.Neighbors;
    
    
end