function smoothcost=setsmoothcost(W,rank_q)
    n=size(W,2);
    smoothcost=zeros(n);
    for iter=1:n-1
        for iter2=iter+1:n
           d=subspace_disparity(W{iter}(:,1:rank_q(iter)),W{iter2}(:,1:rank_q(iter2)));
           smoothcost(iter,iter2)=1-d;%(10^(-5))+(1-d);
           smoothcost(iter2,iter)=1-d;%(10^(-5))+(1-d);
        end
    end
end