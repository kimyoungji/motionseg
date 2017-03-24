function [neighborIds neighborDistances] = k_nearest_neighbors(dataMatrix, queryMatrix, k)

numDataPoints = size(dataMatrix,2);
numQueryPoints = size(queryMatrix,2);

neighborIds = zeros(k,numQueryPoints);
neighborDistances = zeros(k,numQueryPoints);

D = size(dataMatrix, 1); %dimensionality of points

for i=1:numQueryPoints
    d=zeros(1,numDataPoints);
%     for t=1:numDataPoints
%         x=dataMatrix(:,t)-queryMatrix(:,i);
%         d(t)=abs(x'*pseudo_inv(C)*x);
%     end
    for t=1:D % this is to avoid slow repmat()
        d=d+(dataMatrix(t,:)-queryMatrix(t,i)).^2;
    end
    for j=1:k
        [s,t] = min(d);
        neighborIds(j,i)=t;
        neighborDistances(j,i)=sqrt(s);
        d(t) = NaN; % remove found number from d
    end
    
end