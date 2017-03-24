function is_occluded=OcclusionCheck(X,densities)
    N=size(X,2);
    is_occluded=zeros(N,1);

    for iter=1:N
        x1=repmat(X(:,iter),1,10);
        x2=X(:,densities(iter,2:11));
        xx=(x1-x2).*(x1-x2);
        dist=sqrt(sum(xx));
        cc=sum(dist);
         
        if cc<densities(iter,1)*0.9 || cc>densities(iter,1)*1.1
            is_occluded(iter)=1;
        end
    end
    

end