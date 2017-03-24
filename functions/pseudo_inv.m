function [C_inv]=pseudo_inv(C)
    [U,D,V]=svd(C);
    r=find_rank(D);
    s_values=diag(D);
    for iter=1:length(s_values)
        if iter>r
            s_values(iter)=0;
        else
            s_values(iter)=1/s_values(iter);
        end
    end
    C_inv=V*diag(s_values)*U';
end