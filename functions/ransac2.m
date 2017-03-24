function [f inlierIdx outlierIdx corres_1] = ransac2( x,y,tree,data_mat,ransacCoef,funcFindF,funcDist)
%[f inlierIdx] = ransac1( x,y,ransacCoef,funcFindF,funcDist )
%	Use RANdom SAmple Consensus to find a fit from X to Y.
%	X is M*n matrix including n points with dim M, Y is N*n;
%	The fit, f, and the indices of inliers, are returned.
%
%	RANSACCOEF is a struct with following fields:
%	minPtNum,iterNum,thDist,thInlrRatio
%	MINPTNUM is the minimum number of points with whom can we 
%	find a fit. For line fitting, it's 2. For homography, it's 4.
%	ITERNUM is the number of iteration, THDIST is the inlier 
%	distance threshold and ROUND(THINLRRATIO*n) is the inlier number threshold.
%
%	FUNCFINDF is a func handle, f1 = funcFindF(x1,y1)
%	x1 is M*n1 and y1 is N*n1, n1 >= ransacCoef.minPtNum
%	f1 can be of any type.
%	FUNCDIST is a func handle, d = funcDist(f,x1,y1)
%	It uses f returned by FUNCFINDF, and return the distance
%	between f and the points, d is 1*n1.
%	For line fitting, it should calculate the dist between the line and the
%	points [x1;y1]; for homography, it should project x1 to y2 then
%	calculate the dist between y1 and y2.
%	Yan Ke @ THUEE, 20110123, xjed09@gmail.com


minPtNum = ransacCoef.minPtNum;
iterNum = ransacCoef.iterNum;
thInlrRatio = ransacCoef.thInlrRatio;
thDist = ransacCoef.thDist;
ptNum = size(x,2);
thInlr = round(thInlrRatio*ptNum);

inlrNum = zeros(1,iterNum);
fLib = cell(1,iterNum);

% pre_dist=sum((x-y).*(x-y));
% idx1=find(pre_dist>0.001);
% idx2=find(pre_dist(idx1)<1);
% contents=idx1(idx2);
% if numel(contents)<10
%     contents=1:ptNum;
% end
% contents=1:ptNum;
for p = 1:iterNum
	% 1. fit using  random points
	sampleIdx = randIndex(ptNum,minPtNum);
	f1 = funcFindF(y(:,sampleIdx),x(:,sampleIdx));
	% 1.5 find new y
    X1=f1*[x;ones(1,size(x,2))];
    [corres_1,~]=kdtree_nearest_neighbor(tree,X1(1:3,:)');
    y_tmp=data_mat(:,corres_1);
    
	% 2. count the inliers, if more than thInlr, refit; else iterate
 	dist = funcDist(f1,x,y_tmp);
    
	[inlier1,outlier1] = deleteoutliers(dist,thDist);%find(dist < thDist);%
    if mean(dist(inlier1))>mean(dist(outlier1)) && sum(outlier1>10)
        inlier1=outlier1;
    end
%     if length(inlier1)>3%thInlr
%         fLib{p} = funcFindF(y_tmp(:,inlier1),x(:,inlier1));
%     end


% 	[~,inlier1] = find(dist < thDist);%deleteoutliers(dist,thDist);%find(dist < thDist);%
	inlrNum(p) = sum(dist(inlier1));%length(inlier1);%
	if length(inlier1) < thInlr
        inlrNum(p)=10^100;
        continue; 
    end
    if length(inlier1)>3
        fLib{p} = funcFindF(y_tmp(:,inlier1),x(:,inlier1));
    end
    
end


% 3. choose the coef with the most inliers
[~,idx] = min(inlrNum);%changed max->min
f = fLib{idx};
if isempty(f)
    f=eye(4);
end
X1=f*[x;ones(1,size(x,2))];
[corres_1,~]=kdtree_nearest_neighbor(tree,X1(1:3,:)');
y=data_mat(:,corres_1);
dist = funcDist(f,x,y);
% [~,inlierIdx]=find(dist < 0.001);
% [~,outlierIdx]=find(dist >= 0.001);
[inlierIdx,outlierIdx] = deleteoutliers(dist,thDist);%find(dist < thDist);
if mean(dist(inlierIdx))>mean(dist(outlierIdx))
    tmp=inlierIdx;
    inlierIdx=outlierIdx;
    outlierIdx=tmp;
end
    
% % corres_1espondence
% X1=f*[x;ones(1,size(x,2))];
% [corres_1,~]=kdtree_nearest_neighbor(tree,X1(1:3,:)');
end