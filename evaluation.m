% evaluation
acc=[];
[~,gt_pre]=load_ply('./result/right_chair/rchair_gt.ply',1);
for iter=1:310
    % generate gt and res_clusters
    if iter==146
        [~,gt_pre]=load_ply('./result/right_chair/rchair_gt2.ply',1);
    end
    ff=sprintf('./result/right_chair/exp_%d.ply',iter);
    [~,res_cluster_pre]=load_ply(ff,1);
    % erase don't care clusters
    [ids,~]=find(gt_pre>0);
    gt=gt_pre(ids);
    res_clusters=res_cluster_pre(ids);
    g1=max(gt);
    g2=max(res_clusters);
    % generate res
    res={};
    for iter2=1:g2
        [ids,~]=find(res_clusters==iter2);
        res{iter2}=ids;
    end
    
    acc_tmp=seg_acc(gt,res,g1,g2);% measure seg_acc
    acc=[acc;acc_tmp];
end
mean(acc)
sqrt(var(acc))