% clear all;
% Bonn example
addpath('./functions');
addpath('./functions/vbgm');
addpath('./functions/kdtree');
addpath('./functions/estimateRigidTransform');
addpath('./functions/PairwiseMatching');
addpath('./functions/MotionSignatures');
addpath('./functions/exact_alm_rpca/exact_alm_rpca');
addpath('./functions/exact_alm_rpca/exact_alm_rpca/PROPACK');
addpath('./functions/apg/apg');

addpath('./functions/gco-v3.0/matlab');
addpath('./functions/helper_functions');
addpath('./functions/kinect_matlab');

% clouds from img
% IMG=imread('./rigid_dataset/car_d1/60.png');
IMG=imread('./bonn/watercan/200.png');
% IMG=imread('./loader/depth4/90.png');
x1=img_to_xyz(IMG,7,7);
KinectHandles=0;

% % clouds from mat
% load('./synthetic/plane_rot2/pt_1.mat');
% x1=a';

% % clouds from kinect
% addpath('./kinect_matlab/Mex')
% SAMPLE_XML_PATH='./kinect_matlab/Config/RGBDConfig.xml';
% KinectHandles=mxNiCreateContext(SAMPLE_XML_PATH);
% D=mxNiDepth(KinectHandles); D=permute(D,[2 1]);
% mxNiUpdateContext(KinectHandles);
% IMG = D;
% x1=img_to_xyz(IMG,10,10);

% initialize 
group_num=1;
N=size(x1,1);

x1_mesh=delaunay(x1(:,1)',x1(:,2)');
x1_mesh=sort(x1_mesh')';
% plot_mesh(x1_mesh,x1);

%Run EM
[res_clusters,res_neighbors,res_dist]=subspace_em(x1',group_num,x1_mesh);
%save_ply('exp.ply',res_clusters,x1');

%mxNiDeleteContext(KinectHandles);
clear all;

