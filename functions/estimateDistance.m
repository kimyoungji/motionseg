function dist = estimateDistance(T, x, y)

Y=T*[x;ones(1,size(x,2))];
diff=y-Y(1:3,:);
dist=diff.*diff;
dist=sum(dist);