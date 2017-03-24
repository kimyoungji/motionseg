function [simple_M]=MtoSimplifiedM(n,param)
    
    if param.signature(1)==2 || param.signature(1)==8
        [R_a]=RotationMat(param.alpha,param.a);
    end
    if param.signature(1)==8
        [R_b]=RotationMat(param.beta,param.b);
    end
    for iter=1:n
        R_1=eye(4);
        R_2=eye(4);
        if param.signature(1)==2
            tmp_R=R_a(3*(iter-1)+1:3*iter,:);
            R_1=[[tmp_R,param.t_a-tmp_R*param.t_a];[0,0,0,1]];
        elseif param.signature(1)==8
            tmp_Ra=R_a(3*(iter-1)+1:3*iter,:);
            tmp_Rb=R_b(3*(iter-1)+1:3*iter,:);
            R_1=[[tmp_Rb,param.t_b-tmp_Rb*param.t_b];[0,0,0,1]]*[[tmp_Ra,param.t_a-tmp_Ra*param.t_a];[0,0,0,1]];
        end
        if param.signature(2)>0
            td=param.signature(2);
            R_2(:,4)=[param.T*param.t_f(td*(iter-1)+1:td*iter);1];
        end
         simple_M(4*(iter-1)+1:4*iter,:)=R_2*R_1;
    end

end