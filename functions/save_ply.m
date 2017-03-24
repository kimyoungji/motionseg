function save_ply(str,clusters,x)
N=size(x,2);

fid = fopen(str, 'w');
fprintf(fid, 'ply\n');
fprintf(fid, 'format ascii 1.0\n');
fprintf(fid, 'comment PCL generated\n');
fprintf(fid, 'element vertex %d\n',N);
fprintf(fid, 'property float x\n');
fprintf(fid, 'property float y\n');
fprintf(fid, 'property float z\n');
fprintf(fid, 'property uchar red\n');
fprintf(fid, 'property uchar green\n');
fprintf(fid, 'property uchar blue\n');
fprintf(fid, 'element camera 1\n');
fprintf(fid, 'property float view_px\n');
fprintf(fid, 'property float view_py\n');
fprintf(fid, 'property float view_pz\n');
fprintf(fid, 'property float x_axisx\n');
fprintf(fid, 'property float x_axisy\n');
fprintf(fid, 'property float x_axisz\n');
fprintf(fid, 'property float y_axisx\n');
fprintf(fid, 'property float y_axisy\n');
fprintf(fid, 'property float y_axisz\n');
fprintf(fid, 'property float z_axisx\n');
fprintf(fid, 'property float z_axisy\n');
fprintf(fid, 'property float z_axisz\n');
fprintf(fid, 'property float focal\n');
fprintf(fid, 'property float scalex\n');
fprintf(fid, 'property float scaleyply\n');
fprintf(fid, 'property float centerx\n');
fprintf(fid, 'property float centery\n');
fprintf(fid, 'property int viewportx\n');
fprintf(fid, 'property int viewporty\n');
fprintf(fid, 'property float k1\n');
fprintf(fid, 'property float k2\n');
fprintf(fid, 'end_header\n');
%x=x2;
color=round(rand(100,3).*255);
color(1:9,:)=[1,0,0;
    0,1,0;
    0,0,1;
    1,1,0;
    0,1,1;
    1,0,1;
    1,0.2,0.2;
    0.2,1,0.2;
    0.2,0.2,1]*255;
for iter=1:N
    if clusters(iter)==0
        fprintf(fid, '%f %f %f %d %d %d\n', x(1,iter),x(2,iter),x(3,iter),255,255,255);
    elseif clusters(iter)==-1
        fprintf(fid, '%f %f %f %d %d %d\n', x(1,iter),x(2,iter),x(3,iter),0,0,0);
    else
        fprintf(fid, '%f %f %f %d %d %d\n', x(1,iter),x(2,iter),x(3,iter),color(clusters(iter),1),color(clusters(iter),2),color(clusters(iter),3));%magenta    
    end
end
fprintf(fid, '0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 ');
fprintf(fid, '%d 1 0 0',N);
fclose(fid);

end