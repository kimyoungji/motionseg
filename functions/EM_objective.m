function J=EM_objective(package)
    global obj_R;
    global obj_X;
    r=size(package,2);
%     global obj_N;
    [M,N]=size(obj_X);
%     x=obj_X(r+1:end,:);
    W_t=package(1:M,:)';
    m=package(M+1,:)';
    sigma=package(M+2:end,1:r);
    J=0;
    for iter=1:N
        J=J+obj_R(iter)*(log((2*pi)^(-r/2)*norm(sigma)^(-1/2))+((W_t*obj_X(:,iter)-m)'*pinv(sigma)*(W_t*obj_X(:,iter)-m)));
    end

end