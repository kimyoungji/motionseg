function [v2]=OrthoVector(v1)

        for iter=1:length(v1)
            if v1(iter)~=0
                v2=v1;
                iter_plus=iter+1;
                if iter_plus==length(v1)+1
                    iter_plus=1;
                end
                v2(iter)=-v1(iter_plus);
                v2(iter_plus)=v1(iter);
                break;
            end
        end



end