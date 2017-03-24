function dist=subspace_disparity(U,V)
    M=size(U,2);
    N=size(V,2);
    
    if M<N
        V_tmp=V;
        V=U;
        U=V_tmp;
        N_tmp=N;
        N=M;
        M=N_tmp;
    end

        d=0;
        for iter=1:M
            min_v=10^100;
            v_dist=0;
            for iter2=1:N
                %d(u_i,v);
                v_dist=(norm(U(:,iter)-V(:,iter2)));
                if v_dist<min_v
                    min_v=v_dist;
                end
            end
            d=d+(min_v^2);
        end
        dist=sqrt(d)/M;

%         d=0;
%         for iter=1:M
%             max_v=0;
%             for iter2=1:N
%                 v=((U(:,iter)'*V(:,iter2))^2);
%                 if v>max_v
%                     max_v=v;
%                 end
%             end
%             d=d+max_v;
%         end
%         dist=(M-(d))/M;

end