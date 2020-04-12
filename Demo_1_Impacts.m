%% The code provided here belongs to the following paper:
% O. Sarrafzadeh, A. Mehri, H. Rabbani, N. Ghane, A. Talebi, "Circlet based
% framework for red blood cells segmentation and counting", in Proc. IEEE
% Workshop on Signal Processing Systems,
%%% Hangzhou, China, Oct. 14-16, 2015.  
%% The Code was adapted by Richar Marvin Fernández Vílchez
%% Demo_1_Impacts.m
%%% This demo shows how to utilize Circlet Transform (CT) to find desired
%%% cirles in a given image. In this demo, it is assumed that the exact
%%% radius and number of the desire circles are known and we just need to
%%% find the location of the circles.
clear all; close all; clc;
%% read the image
I = imread('images/test001.png');
f = rgb2gray(I);
f = im2double(f);
%% apply the circlet transform
N = 3;                      % number of filters
r0 = [8:26]/2;      % radius range of the candidate circles
[Cl, Gkw] = fdct2(f, N, r0, 'abs');
%% analysis the CT coefficients to find the center of the desire circles
% just the band corresponding to k=2 is analized for lower computations;
% other bands could be used!
pts = maxcoeff(Cl, 2);
pts(:,3) = r0';
%% show the image and draw detected circles on it
figure,imshow(I); hold on;
for r = 1 : size(pts, 1)
    rectangle('Position',[pts(r,1) - pts(r,3), pts(r,2) - pts(r,3),...
        2*pts(r,3), 2*pts(r,3)],...
        'Curvature', [1,1], 'edgecolor', 'b', 'linewidth', 2);
end
hold off;