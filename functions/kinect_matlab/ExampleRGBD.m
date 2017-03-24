addpath('Mex')
%SAMPLE_XML_PATH='Config/SamplesIRConfig.xml';
SAMPLE_XML_PATH='Config/RGBDConfig.xml';
% Start the Kinect Process
KinectHandles=mxNiCreateContext(SAMPLE_XML_PATH);

% To use the Kinect hardware use :
%KinectHandles=mxNiCreateContext(SAMPLE_XML_PATH);

figure;
I=mxNiPhoto(KinectHandles); I=permute(I,[3 2 1]);
D=mxNiDepth(KinectHandles); D=permute(D,[2 1]);
subplot(1,2,1),h1=imshow(I); 
subplot(1,2,2),h2=imshow(D,[0 9000]); colormap('gray');
    
for i=1:10000
    I=mxNiPhoto(KinectHandles); I=permute(I,[3 2 1]);
    D=mxNiDepth(KinectHandles); D=permute(D,[2 1]);
    mxNiUpdateContext(KinectHandles);
    set(h1,'CDATA',I);
    set(h2,'CDATA',D);
    drawnow; 
end

% Stop the Kinect Process
mxNiDeleteContext(KinectHandles);