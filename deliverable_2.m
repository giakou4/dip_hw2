%% 2-2 
clear all; close all; clc;
addpath(genpath('C:\Users\ngiak\OneDrive\dip_hw2_2020\images'));
RGB = imread('im2.jpg');
RGB = imresize(RGB, 0.1);
BW = rgb2gray(RGB);
%% MY DETECT HARRIS FEATURES
display('Running myDetectHarrisFeatures()');
tic
corners = myDetectHarrisFeatures(BW);
figure
imshow(BW);
hold on;
plot(corners(:,2),corners(:,1),'r*');
toc
%% MATLAB DETECT HARRIS FEATURES
display('Running detectHarrisFeatures()');
tic
corners2 = detectHarrisFeatures(BW);
corners2 = corners2.Location;
figure
imshow(BW)
hold on;
plot(corners2(:,1),corners2(:,2),'y*');
toc