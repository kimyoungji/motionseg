function [x,idx]=sample_with_dist(dist,X,num)
    [~,sorted_idx]=sort(dist);
    idx=sorted_idx(1:num);
    x=X(idx,:);
end