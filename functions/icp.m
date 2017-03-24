function [T,corres_1,inlier,outlier] = icp(x, corres_1, data_mat,tree, thratio)
    max_iter=10;
    ransacCoef.minPtNum=4;
    ransacCoef.iterNum=10;
    ransacCoef.thInlrRatio=thratio;
    ransacCoef.thDist=1;
    T=eye(4);
%     tree=kdtree_build(data_mat');
%     [corres_1,~]=kdtree_nearest_neighbor(tree,x');
    for iter=1:max_iter
        y=data_mat(:,corres_1);
        [T_tmp, inlier, outlier, corres_1] = ransac2( x,y,tree,data_mat,ransacCoef,@estimateRigidTransform,@estimateDistance);
        tmp=T_tmp*[x;ones(1,size(x,2))];
        x=tmp(1:3,:);
        T=T_tmp*T;
%         if norm(T_tmp-eye(4))<10^-5
%             break;
%         end
        if icp_convergence(T_tmp)<0.0001
            break;
        end
    end
    
%     tmp=T*[x;ones(1,size(x,2))];
%     [corres_1,dist]=kdtree_nearest_neighbor(tree,tmp(1:3,:)');
%     [inlier,outlier]=deleteoutliers((dist),thdist);
%     if mean(dist(outlier))<mean(dist(inlier))
%         tmp=outlier;
%         outlier=inlier;
%         inlier=tmp;
%     end

end