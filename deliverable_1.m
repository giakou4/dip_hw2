%% 2-1 
clear all; close all; clc; 
addpath(genpath('C:\Users\ngiak\OneDrive\dip_hw2_2020\images'));
RGB = imread('im1.jpg');
RGB = imresize(RGB, 0.1);
BW  = rgb2gray(RGB);
%figure; imshow(BW); title('Initial Image');
BWcanny = edge(BW,'canny');
Drho = 1; Dtheta = 1; n=5;
%% MY HOUGH TRANSFORM
display('Running myHoughTransform()');
tic
[H1, L1, res] = myHoughTransform(BWcanny, Drho, Dtheta, n);
myFunctions.plotHough(H1,L1,-90:Dtheta:89,-floor(norm([size(BW,1) size(BW,2)])):Drho:floor(norm([size(BW,1) size(BW,2)])));
myFunctions.plotLines(BW,L1); title('Hough Lines')
toc
%% MATLAB HOUGH TRANSFORM
% this code was copied from https://www.mathworks.com/help/images/ref/hough.html
display('Running hough(), houghpeaks(), houghlines()');
tic
[H2,theta,rho] = hough(BWcanny,'RhoResolution',Drho,'Theta',-90:Dtheta:89);
figure;
imshow(imadjust(rescale(H2)),[],'XData',theta,'YData',rho,'InitialMagnification','fit');
xlabel('\theta (degrees)'); ylabel('\rho'); title('Hough Transform'); title('Hough Transfrom MATLAB');
axis on, axis normal, hold on;
P2 = houghpeaks(H2,n);
plot(theta(P2(:,2)),rho(P2(:,1)),'s','color','green','MarkerSize',10);
lines2 = houghlines(BW,theta,rho,P2);
L2 = NaN*ones(length(lines2),2);
L2(:,1) = [lines2.rho];
L2(:,2) = [lines2.theta];
myFunctions.plotLines(BW,L2); title('Hough Lines MATLAB');
toc