function X=normalize_r(x)
    X=zeros(size(x));
    for iter=1:size(x,1)
        if sum(x(iter,:))==0
            x(iter,:)=ones(size(x(iter,:))); 
        end
        X(iter,:)=x(iter,:)./sum(x(iter,:));
        
    end
end