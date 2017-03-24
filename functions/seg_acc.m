function [acc]=seg_acc(gt,res,g1,g2)
    aff_mat=zeros(g1,g2);
    corr_mat=zeros(g1,2);
    for iter=1:g1
        for iter2=1:g2
            aff_mat(iter,iter2)=sum(gt(res{iter2})==iter);
        end
        [val,ids]=max(aff_mat(iter,:));
        ol=find(corr_mat(1:end,1)==ids);
        if isempty(ol)
            corr_mat(iter,1)=ids;
            corr_mat(iter,2)=val;
        else
            if val>corr_mat(ol,2)
                corr_mat(iter,1)=ids;
                corr_mat(iter,2)=val;
                corr_mat(ol,1)=0;
                corr_mat(ol,2)=0;
            end
        end
    end
    
    %compute average accuracy
    acc_mat=[];
    for iter=1:g1
        tp=corr_mat(iter,2);
        fp=sum(aff_mat(iter,:))-tp;
        if corr_mat(iter,1)>0
        	fn=sum(aff_mat(:,corr_mat(iter,1)))-tp;
        else
            fn=0;
        end
        acc_mat=[acc_mat;tp/(tp+fp+fn)];
    end
    acc=mean(acc_mat);
%     acc=acc_mat;
end