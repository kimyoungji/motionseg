function [M,labels,nodes,label_table]=construct_affinity_matrix(tree,clust_idx,p_idx,X,dataMatrix,C_inv,A_d,U,m)

    in_idx=10;
    sampled_N=numel(clust_idx);
    k=3;
    
    M=zeros(sampled_N*in_idx);
    nodes=zeros(1,sampled_N*in_idx);
    labels=zeros(1,sampled_N*in_idx);
    label_table=[];
    A=zeros(1,sampled_N*in_idx);
%     C_inv=em_pack.covinv{g_idx};
%     C=em_pack.cov{g_idx};
%     U=em_pack.W_mat{g_idx}(1:3*p_idx+3,em_pack.rank_q(g_idx)+1:end);

    for iter=1:sampled_N
        if nargin==5
            [neighborIds,dist]=kdtree_k_nearest_neighbors(tree,X(k*p_idx+1:k*p_idx+3,clust_idx(iter))',in_idx);
        else
            [candidates,~]=kdtree_k_nearest_neighbors(tree,X(k*p_idx+1:k*p_idx+3,clust_idx(iter))',100);
            dist=zeros(100,1);
            for iter2=1:100
                x=U'*[X(1:3*p_idx,clust_idx(iter));dataMatrix(k*(p_idx-1)+1:k*p_idx,candidates(iter2))]-m;
                dist(iter2)=(x'*C_inv*x)./sqrt(length(C_inv));%./((length(C)/3)^(0.5));%3*(x'*C_inv*x)/length(C);%*A_d
            end
            [dist,n_idx]=sort(dist);
            dist=dist(1:in_idx);
            dist=(max(dist)*ones(in_idx,1)-dist);
            if sum(dist)==0
                dist=ones(in_idx,1);
            end
            neighborIds=candidates(n_idx);
        end

        nodes(in_idx*(iter-1)+1:in_idx*iter)=ones(1,in_idx).*iter;
        for iter2=1:in_idx
            [idx,~]=find(label_table==neighborIds(iter2));
            if isempty(idx)
                idx=length(label_table)+1;
                label_table=[label_table;neighborIds(iter2)];
            end
            labels(in_idx*(iter-1)+iter2)=idx;
            A(in_idx*(iter-1)+iter2)=(dist(iter2));
        end
    end

%     for iter=1:sampled_N
%         x11=repmat( X(1:3,clust_idx(iter)),(sampled_N-iter)*in_idx,in_idx);
%         x12=repmat(X(1:3,clust_idx(iter+1:end)),in_idx,1);
%         x10=x11-repmat(x12(:),1,in_idx);
%         x1 = reshape(sum(reshape(x10,3, (sampled_N-iter)*in_idx*in_idx).^2),(sampled_N-iter)*in_idx,in_idx);
%         x1=sqrt(x1);
%         x21=repmat( dataMatrix(k*(p_idx-1)+1:k*p_idx,label_table(labels(in_idx*(iter-1)+1:iter*in_idx))),(sampled_N-iter)*in_idx,1 );
%         x22=dataMatrix(k*(p_idx-1)+1:k*p_idx,label_table(labels(in_idx*(iter)+1:end)));
%         x20=x21-repmat(x22(:),1,in_idx);
%         x2 = reshape(sum(reshape(x20,3, (sampled_N-iter)*in_idx*in_idx).^2),(sampled_N-iter)*in_idx,in_idx);
%         x2=sqrt(x2);
%         dist=(x1-x2).*(x1-x2);
%         sigma=0.5;
%         tmp_dist=1*((9*sigma^2)-dist);
%         tmp_dist(tmp_dist<0)=0;
%         dist=tmp_dist*100;
%         M(iter*in_idx+1:end,in_idx*(iter-1)+1:in_idx*(iter))=dist;
%         M(in_idx*(iter-1)+1:in_idx*(iter),iter*in_idx+1:end)=dist';
%     end
    M(sub2ind(size(M), 1:in_idx*sampled_N,1:in_idx*sampled_N))=A;
end