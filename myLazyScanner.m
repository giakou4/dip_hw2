% This script reads an image from the current directory and then tries to
% identify sub-images which saves as new images in the current directory
% with the suffix _1, _2, etc.
% Bad coding follow :(
clear all; close all; clc
addpath(genpath('C:\Users\ngiak\OneDrive\dip_hw2_2020\images'));
%% LOAD IMAGE HERE
name = 'im1'; % assume prefix .jpg otherwise see line 10
%% PRE WORK (RGB2GRAY, RESIZE FOR FASTER CALLS, DEFINE SUB-IMAGES SET)
RGB = imread([name '.jpg']);
RGB = imresize(RGB, 0.1);
BW = rgb2gray(RGB);
global IMAGES_FOUND; % stores sub images in 3rd dim
global NO_IMAGES_FOUND; % count number of sub images found
NO_IMAGES_FOUND = 0;
IMAGES_FOUND = zeros(size(RGB,1),size(RGB,2),1);
Drho = 1; Dtheta = 1; n=5; % 5 = {all borders} + 1 = 4 + 1 :have at least 1 line to seperate >2 images
%% RECURSIVE FINDER FTW!
display(['CALLING RECURSIVE FINDER 3']); disp([' ']);
[IMAGES_FOUND,NO_IMAGES_FOUND] = myFunctions.recursiveFinder3(BW,0,Drho,Dtheta,n,IMAGES_FOUND,NO_IMAGES_FOUND);
if NO_IMAGES_FOUND == 1 
    display(['CALLING RECURSIVE FINDER 3 AGAIN']); disp([' ']);
    NO_IMAGES_FOUND = 0;
    IMAGES_FOUND = zeros(size(RGB,1),size(RGB,2),1);
    [IMAGES_FOUND,NO_IMAGES_FOUND] = myFunctions.recursiveFinder3(BW,1,Drho,Dtheta,n,IMAGES_FOUND,NO_IMAGES_FOUND);
end
if NO_IMAGES_FOUND == 1 
    display(['CALLING RECURSIVE FINDER 4']); disp([' ']);
    NO_IMAGES_FOUND = 0;
    IMAGES_FOUND = zeros(size(RGB,1),size(RGB,2),1);
    [IMAGES_FOUND,NO_IMAGES_FOUND] = myFunctions.recursiveFinder4(BW,0,Drho,Dtheta,n,IMAGES_FOUND,NO_IMAGES_FOUND);
end
if NO_IMAGES_FOUND == 1 
    display(['CALLING RECURSIVE FINDER 4 AGAIN']); disp([' ']);
    NO_IMAGES_FOUND = 0;
    IMAGES_FOUND = zeros(size(RGB,1),size(RGB,2),1);
    [IMAGES_FOUND,NO_IMAGES_FOUND] = myFunctions.recursiveFinder4(BW,1,Drho,Dtheta,n,IMAGES_FOUND,NO_IMAGES_FOUND);
end
%% SAVE IMAGES TO ALL_IMAGES{i}
for i=1:NO_IMAGES_FOUND
    curr_im = IMAGES_FOUND(:,:,i);
    res_im = myFunctions.removeBlackLines(curr_im);
    ALL_IMAGES{i} = res_im;
end
clear curr_im; clear res_im;

%% SEARCH FOR DUPLICATES
disp(['Checking for correlation among all images where 0: initial'])
idx = [];
if NO_IMAGES_FOUND > 1
    % compare sub-images with other sub-images
    for i=1:NO_IMAGES_FOUND
        [N1,N2]=size(BW);
         temp = imresize(ALL_IMAGES{i},[N1 N2]);
         R = corr2(BW,temp);
         disp(['(',num2str(0),',',num2str(i),') -> ',num2str(R)])
         if R>0.49 idx=[idx i]; end
    end
    % compate initial image with sub-images    
    for i=1:NO_IMAGES_FOUND
        for j=i+1:NO_IMAGES_FOUND
            [N1,N2]=size(ALL_IMAGES{i});
            temp = imresize(ALL_IMAGES{j},[N1 N2]);
            R = corr2(ALL_IMAGES{i},temp);
            disp(['(',num2str(i),',',num2str(j),') -> ',num2str(R)])
            if R>0.5 idx=[idx j]; end
        end
    end
end
clear temp;
for i=1:length(idx)
    ALL_IMAGES{idx(i)}=[];
end
NO_IMAGES_FOUND2 = 0;
%% PLOT REMAINING IMAGES
for i=1:NO_IMAGES_FOUND
    if norm(size(ALL_IMAGES{i}))>0
        NO_IMAGES_FOUND2 = NO_IMAGES_FOUND2 + 1;
        %ALL_IMAGES{i} = myFunctions.rotImageFinal(ALL_IMAGES{i});
        figure;imshow(ALL_IMAGES{i}/255); title(['Sub-Image No ',num2str(NO_IMAGES_FOUND2)]);
        new_name = [name '_' num2str(NO_IMAGES_FOUND2) '.jpg'];
        imwrite(ALL_IMAGES{i}/255,new_name);
    end
end
NO_IMAGES_FOUND = NO_IMAGES_FOUND2; %clear NO_IMAGES_FOUND2;
if NO_IMAGES_FOUND == 0
    display(' ');
    display('NO SUB IMAGES FOUND');
end