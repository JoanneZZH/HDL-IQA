%% for test one distorted image
%% read the reference image and the distorted image
clear all;
addpath('.\utils');
I_ref_ori = imread('blur-color-1.png');
I_dist_ori = imread('blur-color-4.png');

%% extract features and calculate the resulting score
score = HDL_IQA(I_ref_ori,I_dist_ori)