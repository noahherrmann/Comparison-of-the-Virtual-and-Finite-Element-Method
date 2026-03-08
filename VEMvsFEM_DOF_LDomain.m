%% L2 Error Convegence for increasing DOFs: VEM vs. FEM on L-Domain
clear; close all; clc;

% Set up of directories
root = pwd; 
libDir = 'Used Libraries'; 
addpath(genpath(fullfile(root, libDir, 'ifem')));
addpath(genpath(fullfile(root, libDir, 'vem2D')));
addpath(fullfile(root, 'mesh'));

meshNames = {'l_shape_tri_01.mat', 'l_shape_tri_02.mat', 'l_shape_tri_03.mat', 'l_shape_tri_04.mat'};
nMeshes = length(meshNames);

res = struct('vem_dof', [], 'vem_err', [], 'fem_dof', [], 'fem_err', []);

% Problem Definition
% u(x,y) = sin(3*pi*x) * sin(pi*y)
% f(x,y) = 10*pi^2 * sin(3*pi*x) * sin(pi*y)
sol_vem = @(x,y) sin(3*pi*x) .* sin(pi*y);
sol_ifem = @(p) sin(3*pi*p(:,1)) .* sin(pi*p(:,2));
f_handle = @(x,y) 10*pi^2 * sin(3*pi*x) .* sin(pi*y);

for i = 1:nMeshes
    load(['./mesh/', meshNames{i}]);
    
    % FEM
    node = double(mesh.verts);
    elem = double(reshape(cell2mat(mesh.elems), 3, [])');
    elem = fixorder(node, elem);
    bdFlag = setboundary(node, elem, 'Dirichlet');
    pde_input.d = 1; 
    pde_input.f = @(p) 10*pi^2 * sin(3*pi*p(:,1)) .* sin(pi*p(:,2));
    pde_input.g_D = @(p) zeros(size(p,1), 1);
    
    [soln, ~, ~] = PoissonP2(node, elem, bdFlag, pde_input);
    if isstruct(soln)
        uFEM = soln.u; 
    else
        uFEM = soln;
    end
    
    res.fem_dof(i) = length(uFEM);
    res.fem_err(i) = getL2error(node, elem, sol_ifem, uFEM);
    
    % VEM
    [outVEM, meshVEM] = vem2d(mesh, 2, f_handle);
    A = outVEM.A; b = outVEM.b;
    A(meshVEM.bdDOF,:) = 0; A(:,meshVEM.bdDOF) = 0;
    A(meshVEM.bdDOF,meshVEM.bdDOF) = speye(length(meshVEM.bdDOF));
    b(meshVEM.bdDOF) = 0;
    uVEM = A \ b;
    outVEM.u = uVEM;
    
    res.vem_dof(i) = length(uVEM);
    res.vem_err(i) = getL2Error(outVEM, sol_vem);
end

slopeFEM = polyfit(log(res.fem_dof), log(res.fem_err), 1);
slopeVEM = polyfit(log(res.vem_dof), log(res.vem_err), 1);

figure('Color', 'w', 'Position', [100 100 800 600]);
loglog(res.fem_dof, res.fem_err, '-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
loglog(res.vem_dof, res.vem_err, '-s', 'LineWidth', 2, 'MarkerSize', 8);


grid on;
set(gca, 'FontSize', 14);
xlabel('Number of DOFs', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('L^2-Error', 'FontSize', 16, 'FontWeight', 'bold');
lgd = legend(['FEM (Rate: ', num2str(-slopeFEM(1)), ')'], ...
       ['VEM (Rate: ', num2str(-slopeVEM(1)), ')'], 'Location', 'southwest');
lgd.FontSize = 18;
lgd.FontWeight = 'bold';
title('Comparision of the L^2-Error per DOF for k = 2');
set(gca, 'LooseInset', get(gca, 'TightInset'));