%% VEM Patch Test for k=2, 3, 4
clear; clc; close all;
addpath(genpath("meshes"));
load('voronoi_02.mat');
k_values = [2, 3, 4];

% Define u = x(1-x)y(1-y) = (x^2-x)(y^2-y)
% -Delta u = -[ (2)(y^2-y) + (x^2-x)(2) ] = -2(y^2 - y + x^2 - x)
u_ex = @(x,y) (x.^2 - x) .* (y.^2 - y);
f_rhs = @(x,y) -2 * (y.^2 - y + x.^2 - x);

fprintf('VEM Consistency/Patch Test\n');

for k = k_values
    mesh_k = meshSetup(mesh, k);
    
    % Solve VEM
    [out, mesh_k] = vem2d(mesh_k, k, f_rhs);
    
    K = out.A;
    l = out.b;
    
    % Apply Homogeneous Dirichlet BC(u=0 on boundary)
    K(mesh_k.bdDOF,:) = 0;
    K(:,mesh_k.bdDOF) = 0;
    K(mesh_k.bdDOF,mesh_k.bdDOF) = speye(length(mesh_k.bdDOF));
    l(mesh_k.bdDOF) = 0;
    
    % Solve
    out.u = K \ l;
    errL2 = getL2Error(out, u_ex);
    
    fprintf('Degree k = %d | L2 Error: %.2e\n', k, errL2);
end