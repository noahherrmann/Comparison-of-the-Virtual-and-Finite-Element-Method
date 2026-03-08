%% L2 Error Convegence for increasing DOFs: VEM vs. FEM on Voronoi mesh
clear; close all; clc;

% Set up of directories
root = pwd; 
libDir = 'Used Libraries'; 
addpath(genpath(fullfile(root, libDir, 'ifem')));
addpath(genpath(fullfile(root, libDir, 'vem2D')));
addpath(fullfile(root, 'mesh'));

meshNames = {'voronoi_01.mat', 'voronoi_02.mat', 'voronoi_03.mat', 'voronoi_04.mat'};
nMeshes = length(meshNames);

res = struct('vem_dof', zeros(1,nMeshes), 'vem_err', zeros(1,nMeshes), ...
             'fem_dof', zeros(1,nMeshes), 'fem_err', zeros(1,nMeshes));

% Define Problem
% u(x,y) = |x - 0.5|^2.5
alpha = 2.5;
sol_vem  = @(x,y) abs(x - 0.5).^alpha;
sol_ifem = @(p) abs(p(:,1) - 0.5).^alpha;
f_handle = @(x,y) -(alpha*(alpha-1)) * abs(x - 0.5).^(alpha-2);

for i = 1:nMeshes
    load(['./mesh/', meshNames{i}]); 
    
    % FEM
    [node, elem] = triangulatePolygons(double(mesh.verts), mesh.elems);
    bdFlag = setboundary(node, elem, 'Dirichlet');
    pde_input.f = @(p) -(alpha*(alpha-1)) * abs(p(:,1) - 0.5).^(alpha-2);
    pde_input.d = 1; 
    pde_input.g_D = sol_ifem;
    
    [soln, ~, ~] = PoissonP3(node, elem, bdFlag, pde_input);
    uFEM = soln.u;
    
    res.fem_dof(i) = length(uFEM);
    res.fem_err(i) = getL2error(node, elem, sol_ifem, uFEM);
    
   % VEM
    warning('off')
    k_order = 3;
    quad = quadInfoGL(k_order);
    outVEM = vem2d(mesh, k_order, f_handle); 
    mV = outVEM.mesh; 
    A = outVEM.A;
    b = outVEM.b;
    
    g_D = zeros(length(mV.bdDOF), 1);
    nVerts = size(mesh.verts,1);
    nodes_int = quad.nodes(2:end-1); 
    
    for j = 1:length(mV.bdDOF)
        id = mV.bdDOF(j);
        
        if id <= nVerts
            g_D(j) = sol_vem(mesh.verts(id,1), mesh.verts(id,2));
        else
            offset = id - nVerts - 1;
            edge_idx = floor(offset / (k_order-1)) + 1;
            point_in_edge = mod(offset, k_order-1) + 1;

            v1 = mV.edges(edge_idx, 1);
            v2 = mV.edges(edge_idx, 2);
            P1 = mesh.verts(v1,:);
            P2 = mesh.verts(v2,:);
            xi = nodes_int(point_in_edge);
            P_phys = P1 + xi * (P2 - P1);
            
            g_D(j) = sol_vem(P_phys(1), P_phys(2));
        end
    end
    
    b = b - A(:, mV.bdDOF) * g_D; 
    A(mV.bdDOF,:) = 0; A(:,mV.bdDOF) = 0;
    A(mV.bdDOF,mV.bdDOF) = speye(length(mV.bdDOF));
    b(mV.bdDOF) = g_D;
    
    uVEM = A \ b;
    outVEM.u = uVEM;
    
    res.vem_dof(i) = length(uVEM);
    res.vem_err(i) = getL2Error(outVEM, sol_vem);

    warning('on');
end

% Plot
figure('Color', 'w', 'Position', [100 100 800 600]);
loglog(res.fem_dof, res.fem_err, '-ok', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'k'); hold on;
loglog(res.vem_dof, res.vem_err, '-sr', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'r');

grid on; 
set(gca, 'FontSize', 14);
xlabel('Number of DOFs', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('L^2-Error', 'FontSize', 16, 'FontWeight', 'bold');
lgd = legend('FEM', 'VEM', 'Location', 'southwest');
lgd.FontSize = 18;
lgd.FontWeight = 'bold';

title('Comparision of L^2-Error per DOF for k = 3', 'FontSize', 18);

set(gca, 'LooseInset', get(gca, 'TightInset'));

function [node, elem] = triangulatePolygons(verts, elems)
    node = verts; elem = [];
    for i = 1:length(elems)
        idx = elems{i};
        c = mean(verts(idx,:), 1);
        node = [node; c];
        c_idx = size(node, 1);
        for j = 1:length(idx)
            v1 = idx(j); v2 = idx(mod(j, length(idx)) + 1);
            elem = [elem; v1, v2, c_idx];
        end
    end
end