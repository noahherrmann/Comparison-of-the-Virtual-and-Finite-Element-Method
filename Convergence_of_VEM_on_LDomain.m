%% VEM Test: L-Domain (Poisson, Homogene BCs)
clear; clc; close all;

% Set up of directories
root = pwd; 
libDir = 'Used Libraries'; 
addpath(genpath(fullfile(root, libDir, 'ifem')));
addpath(genpath(fullfile(root, libDir, 'vem2D')));
addpath(fullfile(root, 'mesh'));

demo(1,'l_shape_tri')
set(figure(1), 'Name', 'L Domain - L2 Error vs k');
set(figure(2), 'Name', 'L Domain - H1 Error vs k');
set(figure(3), 'Name', 'L Doamin - L2 Error vs h');

demo(2,'l_shape_tri')