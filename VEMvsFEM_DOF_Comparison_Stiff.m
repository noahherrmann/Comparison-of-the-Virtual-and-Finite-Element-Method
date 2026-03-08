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

res_vem_vor = zeros(length(voronoiNames), 3);
res_fem_vor = zeros(length(voronoiNames), 3);
res_vem_qua = zeros(length(quadNames), 3);
res_fem_qua = zeros(length(quadNames), 3);

% Define Problem
eps_stiff = 0.01; 
u_func = @(x) x - (exp((x-1)/eps_stiff) - exp(-1/eps_stiff))/(1 - exp(-1/eps_stiff));
du_func = @(x) 1 - ( (1/eps_stiff)*exp((x-1)/eps_stiff) ) / (1 - exp(-1/eps_stiff));
u_deriv2 = @(x) -(1/eps_stiff^2 * exp((x-1)/eps_stiff)) / (1 - exp(-1/eps_stiff));

sol_exact = @(x,y) u_func(x) .* u_func(y);
sol_exact_ifem = @(p) u_func(p(:,1)) .* u_func(p(:,2));
du_ex = @(x,y) deal(du_func(x).*u_func(y), u_func(x).*du_func(y));
du_ifem = @(p) [du_func(p(:,1)).*u_func(p(:,2)), u_func(p(:,1)).*du_func(p(:,2))];

f_handle = @(x,y) -(u_deriv2(x).*u_func(y) + u_func(x).*u_deriv2(y));
f_ifem = @(p) -(u_deriv2(p(:,1)).*u_func(p(:,2)) + u_func(p(:,1)).*u_deriv2(p(:,2)));

% Voronoi mesh
for m = 1:length(voronoiNames)
    load(fullfile(meshDir, voronoiNames{m})); 
    [outVEM, meshVEM] = vem2d(mesh, 3, f_handle);
    u_sol = zeros(length(outVEM.b), 1);
    free = setdiff(1:length(outVEM.b), meshVEM.bdDOF);
    u_sol(free) = outVEM.A(free,free) \ (outVEM.b(free)); 
    outVEM.u = u_sol;
    res_vem_vor(m,:) = [length(u_sol), getL2Error(outVEM, sol_exact), getH1Error(outVEM, du_ex)];
    
    [nodeF, elemF] = triangulatePolygons(mesh.verts, mesh.elems);
    nodeF = double(nodeF); elemF = double(elemF);
    bdFlag = setboundary(nodeF, elemF, 'Dirichlet');
    pde_in.f = f_ifem; pde_in.d = 1; pde_in.g_D = @(p) zeros(size(p,1),1);
    [uF_sol, ~, ~] = PoissonP3(nodeF, elemF, bdFlag, pde_in);
    res_fem_vor(m,:) = [length(uF_sol.u), getL2error(nodeF, elemF, sol_exact_ifem, uF_sol.u), ...
                                          getH1error(nodeF, elemF, du_ifem, uF_sol.u)];
end

% Quad mesh
for m = 1:length(quadNames)
    load(fullfile(meshDir, quadNames{m}));
    [outVEM, meshVEM] = vem2d(mesh, 3, f_handle);
    u_sol = zeros(length(outVEM.b), 1);
    free = setdiff(1:length(outVEM.b), meshVEM.bdDOF);
    u_sol(free) = outVEM.A(free,free) \ (outVEM.b(free)); 
    outVEM.u = u_sol;
    res_vem_qua(m,:) = [length(u_sol), getL2Error(outVEM, sol_exact), getH1Error(outVEM, du_ex)];
    
    [nodeF, elemF] = triangulatePolygons(mesh.verts, mesh.elems);
    nodeF = double(nodeF); elemF = double(elemF);
    bdFlag = setboundary(nodeF, elemF, 'Dirichlet');
    [uF_sol, ~, ~] = PoissonP3(nodeF, elemF, bdFlag, pde_in);
    res_fem_qua(m,:) = [length(uF_sol.u), getL2error(nodeF, elemF, sol_exact_ifem, uF_sol.u), ...
                                          getH1error(nodeF, elemF, du_ifem, uF_sol.u)];
end

% Error Comparison and Plot
fprintf('\n================== Error Comparison (eps = %.4f) ==================\n', eps_stiff);
fprintf('%-15s | %-8s | %-12s | %-12s\n', 'Mesh/Method', 'DOFs', 'L2-Error', 'H1-Error');
fprintf('--------------------------------------------------------------------\n');
for i = 1:length(voronoiNames)
    fprintf('Voronoi VEM [%d] | %-8d | %-12.2e | %-12.2e\n', i, res_vem_vor(i,1), res_vem_vor(i,2), res_vem_vor(i,3));
    fprintf('Voronoi FEM [%d] | %-8d | %-12.2e | %-12.2e\n', i, res_fem_vor(i,1), res_fem_vor(i,2), res_fem_vor(i,3));
end
fprintf('--------------------------------------------------------------------\n');
for i = 1:length(quadNames)
    fprintf('Quad VEM [%d]    | %-8d | %-12.2e | %-12.2e\n', i, res_vem_qua(i,1), res_vem_qua(i,2), res_vem_qua(i,3));
    fprintf('Quad FEM [%d]    | %-8d | %-12.2e | %-12.2e\n', i, res_fem_qua(i,1), res_fem_qua(i,2), res_fem_qua(i,3));
end

figure('Color','w','Position',[100 100 850 650]);
loglog(res_fem_vor(:,1), res_fem_vor(:,2), '-ok', 'LineWidth', 2, 'MarkerFaceColor', 'k'); hold on;
loglog(res_vem_vor(:,1), res_vem_vor(:,2), '-sr', 'LineWidth', 2, 'MarkerFaceColor', 'r');
loglog(res_fem_qua(:,1), res_fem_qua(:,2), '--ok', 'LineWidth', 2, 'MarkerSize', 8); 
loglog(res_vem_qua(:,1), res_vem_qua(:,2), '--sr', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
set(gca, 'FontSize', 14);

xlabel('Number of DOFs', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('L^2-Error', 'FontSize', 16, 'FontWeight', 'bold');

lgd = legend('FEM (Voronoi)', 'VEM (Voronoi)', 'FEM (Quad)', 'VEM (Quad)', 'southwest');
lgd.FontSize = 18;
lgd.FontWeight = 'bold';

title('Comparision of L^2-Error per DOF for k = 3', 'FontSize', 18);

set(gca, 'LooseInset', get(gca, 'TightInset'));


function [nodes, triples] = triangulatePolygons(verts, elems)
    triples = []; nodes = verts;
    for i = 1:length(elems)
        poly = elems{i};
        for j = 2:length(poly)-1
            triples = [triples; poly(1), poly(j), poly(j+1)];
        end
    end
end