function [C]=BlkMultiple3by3(A,B)
    size_A=size(A);
    C=[];
    for iter=1:size_A/3
        C=[C;A(3*iter-2:3*iter,:)*B(3*iter-2:3*iter,:)];
    end
end