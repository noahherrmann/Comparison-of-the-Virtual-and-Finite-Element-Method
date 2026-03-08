%% L2 Error Convegence for increasing DOFs: VEM vs. FEM on Voronoi mesh
clear; clc; close all;

% Set up of directories
root = pwd; 
libDir = 'Used Libraries'; 
addpath(genpath(fullfile(root, libDir, 'ifem')));
addpath(genpath(fullfile(root, libDir, 'vem2D')));
addpath(fullfile(root, 'mesh'));

meshNames = {'voronoi_02.mat', 'voronoi_03.mat'}; 
nMeshes = length(meshNames);

res_vem = zeros(nMeshes, 2);
res_fem = zeros(nMeshes, 2);

% Define Problem
sol_exact = @(x,y) sin(pi*x).*sin(pi*y);
sol_exact_ifem = @(p) sin(pi*p(:,1)) .* sin(pi*p(:,2));
f_handle = @(x,y) 2*pi^2 * sin(pi*x) .* sin(pi*y);

for m = 1:nMeshes
    load(['./mesh/', meshNames{m}]);
    
    % VEM
    tic;
    [outVEM, meshVEM] = vem2d(mesh, 2, f_handle);
    A_vem = outVEM.A; b_vem = outVEM.b;
    A_vem(meshVEM.bdDOF,:) = 0; 
    A_vem(:,meshVEM.bdDOF) = 0;
    A_vem(meshVEM.bdDOF,meshVEM.bdDOF) = speye(length(meshVEM.bdDOF));
    b_vem(meshVEM.bdDOF) = 0;
    outVEM.u = A_vem \ b_vem;
    timeVEM = toc;
    errL2_VEM = getL2Error(outVEM, sol_exact);
    res_vem(m, :) = [length(outVEM.u), errL2_VEM];

    % FEM
    tic;
    [nodeFEM, elemFEM] = triangulateVoronoi(mesh.verts, mesh.elems);
    bdFlagFEM = setboundary(nodeFEM, elemFEM, 'Dirichlet');
    pde_input.f = @(p) 2*pi^2 * sin(pi*p(:,1)) .* sin(pi*p(:,2));
    pde_input.d = 1; pde_input.g_D = @(p) zeros(size(p,1),1);
    
    [uFEM_sol, ~, ~] = PoissonP2(nodeFEM, elemFEM, bdFlagFEM, pde_input);
    uFEM = uFEM_sol.u;
    timeFEM = toc;
    errL2_FEM = getL2error(nodeFEM, elemFEM, sol_exact_ifem, uFEM);
    res_fem(m, :) = [length(uFEM), errL2_FEM];
    
    fprintf('VEM: %d DOFs | Error: %.2e | Time: %.4fs\n', res_vem(m,1), errL2_VEM, timeVEM);
    fprintf('FEM: %d DOFs | Error: %.2e | Time: %.4fs\n', res_fem(m,1), errL2_FEM, timeFEM);
end

% Plot
figure('Color', 'w', 'Name', 'Comparision of the L2-Error per DOF for k = 2', 'Position', [100, 100, 800, 600]);
loglog(res_fem(:,1), res_fem(:,2), '-ok', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'k'); hold on;
loglog(res_vem(:,1), res_vem(:,2), '-sr', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');

grid on;
set(gca, 'FontSize', 14);

xlabel('Number of DOFs', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('L^2-Error', 'FontSize', 16, 'FontWeight', 'bold');

lgd = legend('FEM', 'VEM', 'Location', 'southwest');
lgd.FontSize = 18;
lgd.FontWeight = 'bold';

title('Comparision of L^2-Error per DOF for k = 2', 'FontSize', 18);

set(gca, 'LooseInset', get(gca, 'TightInset'));

%% Helper Function
function [nodes, triples] = triangulateVoronoi(verts, elems)
    triples = [];
    nodes = verts;
    for i = 1:length(elems)
        poly = elems{i};
        for j = 2:length(poly)-1
            triples = [triples; poly(1), poly(j), poly(j+1)];
        end
    end
end