%% L-Domain from Triangles
clear; clc; close all;

load('./mesh/triangles_03.mat');
baseMesh = mesh;

v1 = baseMesh.verts;               % [0,1] x [0,1]
v2 = baseMesh.verts + [-1, 0];     % [-1,0] x [0,1]
v3 = baseMesh.verts + [-1, -1];    % [-1,0] x [-1,0]
allVerts = [v1; v2; v3];

n = size(baseMesh.verts, 1);
e1 = baseMesh.elems;
e2 = cellfun(@(x) x + n,   baseMesh.elems, 'UniformOutput', false);
e3 = cellfun(@(x) x + 2*n, baseMesh.elems, 'UniformOutput', false);
allElems = [e1; e2; e3]; 

% Combines nodes on top of each other
[LMesh.verts, ~, ic] = unique(allVerts, 'rows', 'stable');
LMesh.elems = cellfun(@(x) ic(x)', allElems, 'UniformOutput', false);

% Declare boundary nodes
x = LMesh.verts(:,1); y = LMesh.verts(:,2);
is_boundary = (x == 1) | (x == -1) | (y == 1) | (y == -1) | ...
              (x >= 0 & y == 0) | (x == 0 & y <= 0);
LMesh.bndry = find(is_boundary);
% Plot (so we can see if it works)
plotVEMMesh(LMesh);
hold on;
plot(LMesh.verts(LMesh.bndry,1), LMesh.verts(LMesh.bndry,2), 'ro');
title('L Domain from triangle mesh');
axis equal;

% save as .mat
mesh = LMesh;

fileName = './mesh/l_shape_tri_03.mat';
save(fileName, 'mesh');
fprintf('Success! You have saved the file.\n', fileName);