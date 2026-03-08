%% Gaussian function on L-Domain:
clear; close all; clc;

% Set up of directories
root = pwd; 
libDir = 'Used Libraries'; 
addpath(genpath(fullfile(root, libDir, 'ifem')));
addpath(genpath(fullfile(root, libDir, 'vem2D')));
addpath(fullfile(root, 'mesh'));

load('./mesh/l_shape_tri_04.mat'); 
node = double(mesh.verts);
elem = double(reshape(cell2mat(mesh.elems), 3, [])');

% Define Problem
sig = 10; x0 = 0.1; y0 = 0.1;
u_ex = @(x,y) exp(-sig * ((x-x0).^2 + (y-y0).^2));
du_ex = @(x,y) deal(-2*sig*(x-x0).*u_ex(x,y), -2*sig*(y-y0).*u_ex(x,y));
f_ex = @(x,y) 4*sig*u_ex(x,y) .* (1 - sig*((x-x0).^2 + (y-y0).^2));

u_ifem = @(p) exp(-sig * ((p(:,1)-x0).^2 + (p(:,2)-y0).^2));
du_ifem = @(p) [-2*sig*(p(:,1)-x0).*u_ifem(p), -2*sig*(p(:,2)-y0).*u_ifem(p)];
f_ifem = @(p) 4*sig*u_ifem(p) .* (1 - sig*((p(:,1)-x0).^2 + (p(:,2)-y0).^2));

[outVEM, meshVEM] = vem2d(mesh, 2, f_ex);
A_v = outVEM.A; b_v = outVEM.b;
A_v(meshVEM.bdDOF,:) = 0; A_v(:,meshVEM.bdDOF) = 0;
A_v(meshVEM.bdDOF,meshVEM.bdDOF) = speye(length(meshVEM.bdDOF));
b_v(meshVEM.bdDOF) = 0;
outVEM.u = A_v \ b_v;

pde.f = f_ifem; pde.g_D = @(p) zeros(size(p,1),1); pde.d = 1;
bdFlag = setboundary(node, elem, 'Dirichlet');
[soln, eqn, ~] = PoissonP2(node, elem, bdFlag, pde);
uFEM = soln.u;

errL2_VEM = getL2Error(outVEM, u_ex)
errH1_VEM = getH1Error(outVEM, du_ex)
errL2_FEM = getL2error(node, elem, u_ifem, uFEM)
errH1_FEM = getH1error(node, elem, du_ifem, uFEM)

% internal DOFs
freeVEM = setdiff(1:size(outVEM.A, 1), meshVEM.bdDOF);
condVEM = condest(outVEM.A(freeVEM, freeVEM));

% internal DOFs
fixedNodes = find(bdFlag);
allNodesFEM = 1:size(uFEM, 1);
freeFEM = setdiff(allNodesFEM, fixedNodes);
condFEM = condest(eqn.A(freeFEM, freeFEM));

% Plot
figure('Name', 'VEM vs FEM on Gaussian function', 'Color', 'w', 'Position', [100, 100, 1200, 500]);
subplot(1,2,1);
scatter3(meshVEM.verts(:,1), meshVEM.verts(:,2), outVEM.u(1:size(meshVEM.verts,1)), 20, outVEM.u(1:size(meshVEM.verts,1)), 'filled');
hold on;
title(sprintf('VEM L2-Err: %.2e)', errL2_VEM));
xlabel('x'); ylabel('y'); zlabel('u_h');
colormap parula; colorbar; view(3); grid on;


subplot(1,2,2);
N = size(node, 1);
trisurf(elem, node(:,1), node(:,2), uFEM(1:N));
title(sprintf('FEM L2-Err: %.2e)', errL2_FEM));
xlabel('x'); ylabel('y'); zlabel('u_h');
colorbar; view(3); shading interp; grid on;

sgtitle([' Gauss-Peak (\sigma = ', num2str(sig), ')']);