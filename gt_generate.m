%generate ground truth

% IMG=imread('./bonn/chair/170.png');
% x=img_to_xyz(IMG,10,10);
% N=size(x,1);
% clusters=zeros(N,1);
% save_ply('./result/right_chair/gt.ply',clusters,x');

[x,clusters]=load_ply('./result/right_chair/gt.ply',1);
N=size(x,1);
clusters=ones(N,1);
for iter=0:1:1
    ff=sprintf('./result/right_chair/gt%d.ply',iter);
    [tmp_x,tmp_clusters]=load_ply(ff,2);
    n=size(tmp_x,1);
    for iter1=1:n
        for iter2=1:N
            
            if norm(tmp_x(iter1,:)-x(iter2,:))<0.0001
                if iter==0
                    clusters(iter2)=0;
                else
                    clusters(iter2)=iter+1;
                end
            end
        end
    end
end
 save_ply('./result/right_chair/rchair_gt.ply',clusters,x');