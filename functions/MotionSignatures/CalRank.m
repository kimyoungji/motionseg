function [M_rank]=CalRank(M)
    M_rank=0;
    [~,D,~]=svd(M);
    s_value=diag(D);
    
    ratio_mat=[];
    if ~isSimilar(s_value(1),0,0.05)
        for iter=1:length(s_value)-1
           cur_ratio=s_value(iter)/s_value(iter+1);
           ratio_mat=[ratio_mat;cur_ratio]; 
        end
        [~,idx]=max(ratio_mat);
        M_rank=idx;
        if s_value(end)>0.1
            M_rank=length(s_value);
        end
    end

end

function [res]=isSimilar(a,b,epsilon)
    %epsilon=0.001;
    if a<b+epsilon && a>b-epsilon
        res=1;
    else
        res=0;
    end

end