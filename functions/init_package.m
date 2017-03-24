function out_pack=init_package(em_pack)

        out_pack.rank_q=0;%outlier_pack.rank_q;
        out_pack.W_mat=[];%outlier_pack.W_mat;
        out_pack.m_value=[];%outlier_pack.m_value;
        out_pack.cov=[];%outlier_pack.cov;
        out_pack.P_c=0;%0.1;
        out_pack.R=zeros(length(em_pack.R),1);
        out_pack.Mdist=zeros(length(em_pack.R),1);
        out_pack.clusters=zeros(length(em_pack.R),1);
        out_pack.p_clusters=zeros(length(em_pack.R),1);
        out_pack.e_criteria=0;

end