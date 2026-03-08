%% VEM Convergence Test on Voronoi and Hexagons
clear; clc; close all;

% Set up of directories
root = pwd; 
libDir = 'Used Libraries'; 
addpath(genpath(fullfile(root, libDir, 'vem2D')));
addpath(fullfile(root, 'mesh'));

% Run on Voronoi 1,2,3,4 that keep getting smaller with k=2,...,5
demo(1, 'vor'); 
set(figure(1), 'Name', 'Voronoi - L2 Error vs k');
set(figure(2), 'Name', 'Voronoi - H1 Error vs k');
set(figure(3), 'Name', 'Voronoi - L2 Error vs h');

% Run Hexagons 
demo(1, 'hex');
set(figure(4), 'Name', 'Hexagon - L2 Error vs k');
set(figure(5), 'Name', 'Hexagon - H1 Error vs k');
set(figure(6), 'Name', 'Hexagon - L2 Error vs h');

hmax1 = plotVEMMesh("voronoi_01.mat")
hmax2 = plotVEMMesh("voronoi_02.mat")
hmax3 = plotVEMMesh("voronoi_03.mat")

demo(5, 'hex'); %Uses Helmholtz equation to check mass matrix
set(figure(10), 'Name', 'Helmholtz - L2 Error vs k');
set(figure(11), 'Name', 'Helmholtz - H1 Error vs k');
set(figure(12), 'Name', 'Helmholtz - L2 Error vs h');
