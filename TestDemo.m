% Texture Segmentation using Active contours
% Reference : "A Novel Active Contour Model for Texture Segmentation, A. Tatu, S. Bansal"
close all;
clear;
clc;
I=load('test_image.mat');
I=I.im;
opt.Isize=[100 100];
opt.lambda=.8;
opt.mu=.00000001;
opt.dim_patch=11;
opt.max_itr=10;
tic;
[I,phi]=texture_seg_split_bregman( I,opt );
toc;
imshow(I);hold on;[C,h] = contour(phi,'r','LineWidth',1);
savefig('Result.fig');