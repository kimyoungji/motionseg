function [x,clusters]=load_ply(str,ver)
    fid = fopen(str, 'r');
    for iter=1:3
        tline = fgets(fid);
    end
    tline= textscan(fid, '%s %s %d\n');%sscanf(fgetl(fid), '%s %s %d\n');%sscanf(fid, '%c %c %d\n');
    N=tline{3};
if ver==1
    for iter=1:29
        tline = fgets(fid);
    end
else
    for iter=1:6
        tline = fgets(fid);
    end
end
%x=x2;
color(1:9,:)=[1,0,0;
    0,1,0;
    0,0,1;
    1,1,0;
    0,1,1;
    1,0,1;
    1,0.2,0.2;
    0.2,1,0.2;
    0.2,0.2,1]*255;
clusters=zeros(N,1);
x=[];
for iter=1:N
    if ver==1
        t_line=sscanf(fgetl(fid), '%f %f %f %d %d %d\n');
    else
        t_line=sscanf(fgetl(fid), '%f %f %f\n');
    end
    x=[x;t_line(1:3)'];
    if ver==1
        c_color=t_line(4:6);
        if c_color(1)==color(1,1) && c_color(2)==color(1,2) && c_color(3)==color(1,3)
            clusters(iter)=1;
        elseif c_color(1)==color(2,1) && c_color(2)==color(2,2) && c_color(3)==color(2,3)
            clusters(iter)=2;
        elseif c_color(1)==color(3,1) && c_color(2)==color(3,2) && c_color(3)==color(3,3)
            clusters(iter)=3;
        elseif c_color(1)==color(4,1) && c_color(2)==color(4,2) && c_color(3)==color(4,3)
            clusters(iter)=4;
        elseif c_color(1)==color(5,1) && c_color(2)==color(5,2) && c_color(3)==color(5,3)    
            clusters(iter)=5;
        elseif c_color(1)==color(6,1) && c_color(2)==color(6,2) && c_color(3)==color(6,3) 
            clusters(iter)=6;
        elseif c_color(1)==color(7,1) && c_color(2)==color(7,2) && c_color(3)==color(7,3)
            clusters(iter)=7;
        elseif c_color(1)==color(8,1) && c_color(2)==color(8,2) && c_color(3)==color(8,3)
            clusters(iter)=8;
        elseif c_color(1)==color(9,1) && c_color(2)==color(9,2) && c_color(3)==color(9,3)
            clusters(iter)=9;
        else
            clusters(iter)=0;
        end
    end
end
fclose(fid);

end