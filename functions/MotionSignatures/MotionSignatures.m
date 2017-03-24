function [simple_M,parameters]=MotionSignatures(M)

% % set Motion matrix
%     M=[];
%     init_M=[];
%     for iter=1:30
%         theta=5*abs(randn(1))*pi/180;
%         theta2=5*abs(randn(1))*pi/180;
%         theta3=randn(1);
%         theta4=abs(randn(1));
%         a_init=normc([0;1;0]);
%         R_a=RotationMat(theta2,a_init);
%         b_init=normc([1;0;0]);
%         R_b=RotationMat(theta,b_init);
%         R_part=R_b*R_a-eye(3,3);
%         t_a_init=[10;0;0];
%         t_b_init=[0;1;0];
%         tt=[theta4;0;0];
% %         t_part=tt+R_b*t_a_init-R_b*R_a*t_a_init+t_b_init-R_b*t_b_init;
% %         if iter<16
% %             R_part=R_b;%-eye(3,3);
% %             t_part=t_b_init-R_b*t_b_init;
% %         else
% %             R_part=R_a;%-eye(3,3);
% %             t_part=t_a_init-R_a*t_a_init;
% %         end
% %         
% %         init_M=[init_M;[R_part,t_part;0,0,0,1]];
% %         R_part=R_part-eye(3);
%         t_part=tt;
%         R_part=zeros(3,3);
%         M=[M;R_part(:)',t_part'];    
%     end
%     M=M+0.0001*rand(30,12);

%% select type of motion (signature)
    signature=zeros(1,2);
    [~,D,~]=svd(M(:,1:9));
    s_values=diag(D);
    %tt_s=sum(s_values>0.1);
%     sig1=find_rank2(s_values(1:end));%find_rank(s_values);%
    sig1=rank(s_values(1:end),1);
    %sig2=find_rank(s_values);
    %signature(1)=sig1;
    %sig1=sig1+tt_s;
    
    if sig1<1
        signature(1)=0;
    elseif sig1<3
        signature(1)=2;
    elseif sig1<9
        signature(1)=8;
    end
    

    [~,D,~]=svd(M);
    s_values1=diag(D);
    if sig1>=9
        signature(2)=rank(s_values1,1)-sig1;
    else
        signature(2)=rank(s_values1(1:sig1+3),1)-sig1;
    end
    if signature(2)<0
        signature(2)=0;
    elseif signature(2)>3
        signature(2)=3;
    end
    
%     M_19=exact_alm_rpca(M(:,1:9),signature(1));
%% initialization
    a=[0;0;0];
    b=[0;0;0];
    alpha=[];
    beta=[];
    t_a=[0;0;0];
    t_b=[0;0;0];
    T=zeros(3,2);
    t_f=zeros(15,1);
    seq_num=size(M,1);
    simple_M=repmat(eye(4),seq_num,1);
%% rotation axes, angles and T
    
    if signature(1)==2
        %axis
        R_part=[M(:,1),M(:,4),M(:,7)];
        R_part=[R_part;M(:,2),M(:,5),M(:,8)];
        R_part=[R_part;M(:,3),M(:,6),M(:,9)];
%         A2=exact_alm_rpca(R_part,2);
        [~,~,V]=svd(R_part);
        a=V(:,3);
        a=normc(a);
        %angles
        c=normc(cross(a,normc(a+[1;1;1])));
        c_cross=[0,-c(3),c(2);c(3),0,-c(1);-c(2),c(1),0];
        R_a=reshape(M(:,1:9)',3,3*length(M))+repmat(eye(3,3),1,length(M));
        sin_mat=reshape(a'*c_cross*R_a,3,length(M))'*c;
        cos_mat=reshape(c'*R_a,3,length(M))'*c;
        cos_mat(cos_mat>1)=1;
        cos_mat(cos_mat<-1)=-1;
        sin_mat(sin_mat>1)=1;
        sin_mat(sin_mat<-1)=-1;
        alpha=sign(asin(sin_mat)).*abs(acos(cos_mat));
        
    elseif signature(1)==8
        %axes  
%         [~,~,V]=svd(M(:,1:9));
%         tmp_d=0;
%         N=reshape(V(:,9),3,3);
%         for iter=3:9
%             tmp_N=reshape(V(:,iter),3,3);
%             [~,d,~]=svd(tmp_N);
%             if d(1)>tmp_d
%                 tmp_d=d(1);
%                 N=tmp_N;
%                 ss=iter;
%             end 
%         end

        ss=9;
        M2=exact_alm_rpca(M(:,1:9),ss-1);
        [~,~,V]=svd(M2);
%         [~,~,V]=svd(M(:,1:9));
        N=reshape(V(:,9),3,3);
        
        N2=exact_alm_rpca(N,1);
        [u,d,v]=svd(N2);
        
        b=u(:,1);
        a=v(:,1);
        a=normc(a);%./norm(a);
        b=normc(b);%./norm(b);
        
        if abs(a'*b)>0.9
%             b=u(:,2);
            signature=[0,0];
        end
        
        %angles
        a_cross=[0,-a(3),a(2);a(3),0,-a(1);-a(2),a(1),0];
        b_cross=[0,-b(3),b(2);b(3),0,-b(1);-b(2),b(1),0];
        R_f_1=M(:,1:3)';
        R_f_2=M(:,4:6)';
        R_f_3=M(:,7:9)';
        R_f=[R_f_1(:)';R_f_2(:)';R_f_3(:)']+repmat(eye(3,3),1,length(M));
        R_f_transpose=reshape(M(:,1:9)',3,3*length(M))+repmat(eye(3,3),1,length(M));
        cos_alpha=(reshape(b'*R_f,3,length(M))'*b-repmat((b'*a)^2,length(M),1))./(1-(b'*a)^2);
        cos_beta=(reshape(a'*R_f,3,length(M))'*a-repmat((a'*b)^2,length(M),1))./(1-(a'*b)^2);
        sin_beta=(reshape(b'*a_cross*R_f,3,length(M))'*a)./(1-(a'*b)^2);
        sin_alpha=-(reshape(a'*b_cross*R_f_transpose,3,length(M))'*b)./(1-(b'*a)^2);
        cos_alpha(cos_alpha>1)=1;
        cos_alpha(cos_alpha<-1)=-1;
        cos_beta(cos_beta>1)=1;
        cos_beta(cos_beta<-1)=-1;
        sin_alpha(sin_alpha>1)=1;
        sin_alpha(sin_alpha<-1)=-1;
        sin_beta(sin_beta>1)=1;
        sin_beta(sin_beta<-1)=-1;
        alpha=-sign(asin(sin_alpha)).*abs(acos(cos_alpha));
        beta=-sign(asin(sin_beta)).*abs(acos(cos_beta)); 
%         else
%             signature(1)=9;
%         end
    end
    
    [U,~,~]=svd(M(:,1:9));
%     [~,~,V]=svd([U(:,1:signature(1)),M(:,10:12)]);

    if signature(2)==1
        M2=exact_alm_rpca([U(:,1:sig1),M(:,10:12)],1);
        [~,~,V]=svd(M2);
        [~,~,V]=svd([U(:,1:sig1),M(:,10:12)]);
        orth_T=normc(V(end-2:end,end-1:end));
        T=normc(cross(orth_T(:,1),orth_T(:,2)));
%         T=normc(V(end-2:end,end));
    elseif signature(2)==2
        M2=exact_alm_rpca([U(:,1:sig1),M(:,10:12)],2);
        [~,~,V]=svd(M2);
        orth_T=normc(V(end-2:end,3));
        T1=normc(cross(orth_T,normc(orth_T+[1;1;1])));%OrthoVector(orth_T);
        T2=normc(cross(orth_T,T1));
        T=[normc(T1),normc(T2)];
    elseif signature(2)==3
        T=eye(3);
    end

%% t_a, t_b and tilda(t_f)    
    q_tmp=M(:,10:12)';
    Q=q_tmp(:);
    T_mat=[];
    if signature(2)>0
        T_tmp=cell(1,length(M));
        [T_tmp{:}]=deal(T);
        T_mat=full(blkdiag(T_tmp{:}));
    else
        T_mat=[];
    end
    
    if signature(1)==0 || signature(1)==9
        if signature(2)>0
        x=pinv(T_mat'*T_mat)*T_mat'*Q;
        t_f=x;
        else
           t_f=[]; 
        end
    elseif signature(1)==2
        a_ortho1=normc(cross(a,normc(a+[1;1;1])));
        a_ortho2=normc(cross(a,a_ortho1));
        A_tmp=[a_ortho1,a_ortho2];
        [R_a]=RotationMat(alpha,a);

        A=(repmat(eye(3,3),length(M),1)-R_a)*A_tmp;
        P=[A,T_mat];
        x=pinv(P'*P)*P'*Q;
        t_a=A_tmp*x(1:2);
        t_f=x(3:end);
    elseif signature(1)==8
        a_ortho1=normc(cross(a,normc(a+[1;1;1])));
        a_ortho2=normc(cross(a,a_ortho1));
        A_tmp=[a_ortho1,a_ortho2];
        [R_a]=RotationMat(alpha,a);
        
        b_ortho1=normc(cross(b,normc(b+[1;1;1])));
        b_ortho2=normc(cross(b,b_ortho1));
        B_tmp=[b_ortho1,b_ortho2];
        [R_b]=RotationMat(beta,b);
        
        B=(repmat(eye(3,3),length(M),1)-R_b)*B_tmp;
        A=BlkMultiple3by3(R_b,(repmat(eye(3,3),length(M),1)-R_a))*A_tmp;
        
        P=[A,B,T_mat];
        x=pinv(P'*P)*P'*Q;
        t_a=A_tmp*x(1:2);
        t_b=B_tmp*x(3:4);
        t_f=x(5:end);
    end

%% set parameters
    parameters.signature=signature;
    parameters.a=a;
    parameters.b=b;
    parameters.alpha=alpha;
    parameters.beta=beta;
    parameters.t_a=t_a;
    parameters.t_b=t_b;

    parameters.T=T;
    parameters.t_f=t_f;

%% set simple_M
    [simple_M]=MtoSimplifiedM(15,parameters);
%     for iter=1:size(M,1)
%         R_1=eye(4);
%         R_2=eye(4);
%         if signature(1)==2
%             tmp_R=R_a(3*(iter-1)+1:3*iter,:);
%             R_1=[[tmp_R,t_a-tmp_R*t_a];[0,0,0,1]];
%         elseif signature(1)==8
%             tmp_Ra=R_a(3*(iter-1)+1:3*iter,:);
%             tmp_Rb=R_b(3*(iter-1)+1:3*iter,:);
%             R_1=[[tmp_Rb,t_b-tmp_Rb*t_b];[0,0,0,1]]*[[tmp_Ra,t_a-tmp_Ra*t_a];[0,0,0,1]];
%         end
%         if signature(2)>0
%             td=signature(2);
%             R_2(:,4)=[T*t_f(td*(iter-1)+1:td*iter);1];
%         end
%          simple_M(4*(iter-1)+1:4*iter,:)=R_2*R_1;
%     end

    
    
   

end