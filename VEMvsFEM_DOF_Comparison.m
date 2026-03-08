%% L2 Error Convegence for increasing DOFs: VEM vs. FEM on Voronoi and Quad
clear; clc; close all;

% Set up of directories
root = pwd; 
libDir = 'Used Libraries'; 
addpath(genpath(fullfile(root, libDir, 'ifem')));
addpath(genpath(fullfile(root, libDir, 'vem2D')));
addpath(fullfile(root, 'mesh'));

meshDir = './mesh'; 
voronoiNames = {'voronoi_01.mat', 'voronoi_02.mat'}; 
quadNames    = {'quad_10.mat', 'quad_20.mat'};

res_vem_vor = zeros(length(voronoiNames), 2);
res_fem_vor = zeros(length(voronoiNames), 2);
res_vem_qua = zeros(length(quadNames), 2);
res_fem_qua = zeros(length(quadNames), 2);

% Define Problem
sol_exact = @(x,y) sin(pi*x).*sin(pi*y);
sol_exact_ifem = @(p) sin(pi*p(:,1)) .* sin(pi*p(:,2));
f_handle = @(x,y) 2*pi^2 * sin(pi*x) .* sin(pi*y);

% Voronoi mesh
for m = 1:length(voronoiNames)
    load(fullfile(meshDir, voronoiNames{m})); 
    
    % VEM
    [outVEM, meshVEM] = vem2d(mesh, 3, f_handle);
    u_sol = zeros(length(outVEM.b), 1);
    free = setdiff(1:length(outVEM.b), meshVEM.bdDOF);
    u_sol(free) = outVEM.A(free,free) \ (outVEM.b(free)); 
    outVEM.u = u_sol;
    res_vem_vor(m,:) = [length(outVEM.u), getL2Error(outVEM, sol_exact)];
    
    % FEM
    [nodeF, elemF] = triangulatePolygons(mesh.verts, mesh.elems);
    bdFlag = setboundary(nodeF, elemF, 'Dirichlet');
    pde_in.f = @(p) 2*pi^2 * sin(pi*p(:,1)) .* sin(pi*p(:,2));
    pde_in.d = 1; pde_in.g_D = @(p) zeros(size(p,1),1);
    [uF_sol, ~, ~] = PoissonP3(nodeF, elemF, bdFlag, pde_in);
    res_fem_vor(m,:) = [length(uF_sol.u), getL2error(nodeF, elemF, sol_exact_ifem, uF_sol.u)];
end

% Quad mesh
for m = 1:length(quadNames)
    load(fullfile(meshDir, quadNames{m}));

    % VEM
    [outVEM, meshVEM] = vem2d(mesh, 3, f_handle);
    u_sol = zeros(length(outVEM.b), 1);
    free = setdiff(1:length(outVEM.b), meshVEM.bdDOF);
    u_sol(free) = outVEM.A(free,free) \ (outVEM.b(free)); 
    outVEM.u = u_sol;
    res_vem_qua(m,:) = [length(outVEM.u), getL2Error(outVEM, sol_exact)];
    
    % FEM
    [nodeF, elemF] = triangulatePolygons(mesh.verts, mesh.elems);
    bdFlag = setboundary(nodeF, elemF, 'Dirichlet');
    [uF_sol, ~, ~] = PoissonP3(nodeF, elemF, bdFlag, pde_in);
    res_fem_qua(m,:) = [length(uF_sol.u), getL2error(nodeF, elemF, sol_exact_ifem, uF_sol.u)];
end

% Plot
figure('Color','w','Position',[100 100 850 650]);
% Voronoi
loglog(res_fem_vor(:,1), res_fem_vor(:,2), '-ok', 'LineWidth', 2, 'MarkerFaceColor', 'k'); hold on;
loglog(res_vem_vor(:,1), res_vem_vor(:,2), '-sr', 'LineWidth', 2, 'MarkerFaceColor', 'r');
% Quad
loglog(res_fem_qua(:,1), res_fem_qua(:,2), '--ok', 'LineWidth', 2, 'MarkerSize', 8); 
loglog(res_vem_qua(:,1), res_vem_qua(:,2), '--sr', 'LineWidth', 2, 'MarkerSize', 8);

grid on; set(gca, 'FontSize', 14);
xlabel('DOFs'); ylabel('L2-Error');
legend('FEM (Voronoi)', 'VEM (Voronoi)', 'FEM (Quad)', 'VEM (Quad)', 'Location', 'southwest', 'FontSize', 16);
title('Comparision of L^2 Error per DOF for k = 3');

function [nodes, triples] = triangulatePolygons(verts, elems)
    triples = []; nodes = verts;
    for i = 1:length(elems)
        poly = elems{i};
        for j = 2:length(poly)-1
            triples = [triples; poly(1), poly(j), poly(j+1)];
        end
    end
end