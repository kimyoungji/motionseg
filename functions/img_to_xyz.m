function pcloud=img_to_xyz(IMG,coef1,coef2)
    
    data_pre=depthToCloud(IMG);

    t1=1:size(IMG,1)/coef1;
    t2=1:size(IMG,2)/coef2;    
    data(:,:,1)=data_pre(t1*coef1,t2*coef2,1);
    data(:,:,2)=data_pre(t1*coef1,t2*coef2,2);
    data(:,:,3)=data_pre(t1*coef1,t2*coef2,3);
    
%     data=bilateral3(data,3,3,5,1,1);
%     data=depthToCloud(IMG);
    xx=data(:,:,1);
    x=[];
    x=[x,xx(:)];
    xx=data(:,:,2);
    x=[x,xx(:)];
    xx=data(:,:,3);
    x=[x,xx(:)];
    pcloud=x(~isnan(x(:,1)),:);
end