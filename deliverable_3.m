%% 2-3
clear all;close all; clc;
addpath(genpath('C:\Users\ngiak\OneDrive\dip_hw2_2020\images'));
RGB = imread('im2.jpg');
RGB = imresize(RGB, 0.1);
I = RGB;
angle1 = 54; angle2 = 213; 
%% MY IMAGE ROTATION
display('Running myImgRotation()');
tic
rotI1 = myImgRotation(I ,angle1);
figure; imshow(rotI1/255);
toc
tic
rotI2 = myImgRotation(I,angle2);
figure; imshow(rotI2/255);
toc
%% MATLAB IMAGE ROTATION
display('Running imrotate()');
tic
rotI3 = imrotate(I ,angle1);
figure; imshow(rotI3);
toc
tic
rotI4 = imrotate(I ,angle2);
figure; imshow(rotI4);
toc