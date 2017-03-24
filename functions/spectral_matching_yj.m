function [sol]  = spectral_matching_yj(M, labels, nodes)


tic;

n = length(labels);

v = ones(length(nodes),1);

v = v/norm(v);

nNodes = max(nodes);

nLabels = max(labels);

iterClimb_eigen = 30;

%% compute the first eigenvector (iterative power method)

for i = 1:iterClimb_eigen
  
  v = M*v;

  v = v/norm(v);
  
end

%% double-stochastic normalization

v0 = v;
v1 = v;

for k = 1:20

    for j = 1:nNodes

        f = find(nodes == j);

        v1(f) = v0(f)/(sum(v0(f))+eps);

    end

    for j = 1:nLabels

        f = find(labels == j);

        v0(f) = v1(f)/(sum(v1(f))+eps);

    end

end

v = (v1+v0)/2;
%% --------------------------------------
%initialization
x_star = v;
x = zeros(n,1);
L = [1:n]';
while 1
    %step 3
    max_num=0;
    max_idx=0;
    for iter=1:length(L)
        if x_star(L(iter))>max_num
            max_num=x_star(L(iter));
            max_idx=L(iter);
        end
    end
    if max_num>0
       x(max_idx)=1;
    else
        break;
    end
    %step 4
    ii=nodes(max_idx);
    jj=labels(max_idx);
    
    new_L=[];
    for iter=1:length(L)
      if nodes(L(iter))~=ii && labels(L(iter))~=jj
           new_L=[new_L;L(iter)];
      end
    end
    L=new_L;
    %step 5
    if isempty(L)
        break;
    end
end
sol=x;
end