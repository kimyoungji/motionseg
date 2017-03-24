function s_pack=split_clusters(package,X)
    %% split
    [~,N]=size(X);
    group_num=numel(package.rank_q);
    ori_dist=package.R';%([package.R';(X)]);%
    s_pack=init_package(package);
    new_idx=0;
    for iter=1:group_num
        pre_idx=find(package.clusters==iter);
        if ((numel(pre_idx)>2)) %&& package.e_criteria(iter)>10^(-4)) || group_num==1
            
            if group_num==1
                [label,~] = vbgm((X(:,pre_idx)), 2);
            else
                [label,~] = vbgm(ori_dist(:,pre_idx), 2);
            end

            if max(label)==1 %|| package.e_criteria(iter)<10^(-1) %|| sum(label==2)<10
                new_idx=new_idx+1; 
                tmp_idx=pre_idx(label==1);
                iter_idx=tmp_idx(package.p_clusters(tmp_idx)~=0);
                if numel(iter_idx)<3
                    new_idx=new_idx-1;
                else
                    s_pack.clusters(tmp_idx)=new_idx;
                    s_pack.p_clusters(iter_idx)=new_idx;
                    s_pack.R(:,new_idx)=package.R(:,iter);
                    s_pack.Mdist(:,new_idx)=package.Mdist(:,iter);
                end
            else
                for iter2=1:max(label)
                    new_idx=new_idx+1;
                    tmp_idx=pre_idx(label==iter2);
                    iter_idx=tmp_idx(package.p_clusters(tmp_idx)~=0);
                    if numel(iter_idx)<3
                        new_idx=new_idx-1;
                    else
                        s_pack.clusters(tmp_idx)=new_idx;
                        s_pack.p_clusters(iter_idx)=new_idx;
                        s_pack.R(:,new_idx)=zeros(N,1);
                        s_pack.R(iter_idx,new_idx)=1;
                        s_pack.Mdist(:,new_idx)=s_pack.R(:,new_idx);
                    end
                end
            end
        else
            new_idx=new_idx+1; 
            tmp_idx=pre_idx;
            iter_idx=tmp_idx(package.p_clusters(tmp_idx)~=0);
            if numel(iter_idx)<4
                new_idx=new_idx-1;
            else
                s_pack.clusters(tmp_idx)=new_idx;
                s_pack.p_clusters(iter_idx)=new_idx;
                s_pack.R(:,new_idx)=package.R(:,iter);
                s_pack.Mdist(:,new_idx)=package.Mdist(:,iter);
            end
        end
    end
    for iter=1:N
        s_sum=sum(s_pack.R(iter,:));
        if s_sum~=0
            s_pack.R(iter,:)=s_pack.R(iter,:)/s_sum;
        else
            s_pack.R(iter,:)=ones(1,new_idx)/new_idx;
%             s_pack.clusters(iter)=0;
%             s_pack.p_clusters(iter)=0;
        end
    end
%     a=(sum(s_pack.R'))';
%     s_pack.R=s_pack.R./repmat(a,1,new_idx);

end