function [data] = make_dpoints(m,pose,xyz)
    % m : width of grid
    % pose : [  min, max]
    % xyz          
    %       
    
    [X,Y] = meshgrid(linspace(pose(1,1),pose(1,2),m), linspace(pose(2,1),pose(2,2),m));
    
    X = reshape(X,1,[]);
    Y = reshape(Y,1,[]);
    Z = zeros(size(X));%sin(X).*cos(Y);%ones(1,m^2);%
    data=[X;Y;Z];%+randn(3,m^2)/1000;
    data=[data(xyz(1),:);data(xyz(2),:);data(xyz(3),:)];
    data=[data;ones(1,m^2)];
    
    
    
end