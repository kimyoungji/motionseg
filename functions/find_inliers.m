function r=find_inliers(X)
%     x1=X(1:end-1);
%     x2=X(2:end);
%     dif=x2-x1;
%     [~,r]=max(dif);
    r=numel(X);
    for iter=1:numel(X)-1
        if X(iter)>10%X(iter+1)/X(iter)>1.5
            r=iter;
            break;
        end
    end

end