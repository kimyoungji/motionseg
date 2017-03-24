function [x_mesh]= rep_x_mesh(x_mesh,c_idx,prev_max)
    
    member_sum=ones(1,size(x_mesh,1));
    d_idx=ismember(x_mesh(:,1),c_idx);
    %d_idx=d_idx(d_idx>0);
    member_sum(:,d_idx)=0;
    d_idx=ismember(x_mesh(:,2),c_idx);
    %d_idx=d_idx(d_idx>0);
    member_sum(:,d_idx)=0;
    d_idx=ismember(x_mesh(:,3),c_idx);
    %d_idx=d_idx(d_idx>0);
    member_sum(:,d_idx)=0;
    x_mesh=x_mesh(member_sum==1,:);

    id_table=zeros(1,prev_max);
    cur=0;
    for iter=1:prev_max
        if sum((c_idx==iter))==0
            cur=cur+1;
            id_table(iter)=cur;
        end
    end
    x_mesh(:,1)=id_table(x_mesh(:,1));
    x_mesh(:,2)=id_table(x_mesh(:,2));
    x_mesh(:,3)=id_table(x_mesh(:,3));
end