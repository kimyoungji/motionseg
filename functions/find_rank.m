function [rr]=find_rank(s_values)
%     s_values=s_values./s_values(1);
    N=length(s_values);
    
    if N==0
        rr=0;
    elseif N==1
        rr=1;
    elseif sum(s_values)==0
        rr=0;
    else
        last_n=N;
        start_n=0;
        if s_values(1)>0.1
            start_n=1;
        end
        for iter=1:N
            if s_values(iter)>1
                start_n=iter;
            end
            if s_values(iter)<10^(-3)
                last_n=iter+1;
                break;
            end
        end
        if last_n>N
            last_n=N;
        end
        
        if last_n>0 && start_n>0 && start_n<=last_n
            X=abs(s_values(start_n:last_n-1)./s_values(start_n+1:last_n));
            [m,rr]=max(X);
            rr=rr+start_n-1;
            if m<3
                rr=N;
            end
        else
            rr=start_n;
        end
%         if s_values(1)<0.1
%             rr=0;
%         end

        
    end
    
end