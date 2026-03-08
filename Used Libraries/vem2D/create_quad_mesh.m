%% Creates and plots Quad Meshes
clear; clc; close all;

% Refinement (10x10 and 20x20)
resolutions = [10, 20];

meshDir = './mesh';

for N = resolutions

    x = linspace(0, 1, N+1);
    [X, Y] = ndgrid(x, x);
    mesh.verts = [X(:), Y(:)];
    
    k_idx = 1;
    mesh.elems = cell(N*N, 1);
    totalArea = 0;
    
    for j = 1:N 
        for i = 1:N 
            % Enumerate vertices counter clockwise
            n1 = (j-1)*(N+1) + i;
            n2 = (j-1)*(N+1) + i + 1;
            n3 = j*(N+1) + i + 1;
            n4 = j*(N+1) + i;            
            nodes = [n1; n2; n3; n4];
            mesh.elems{k_idx} = nodes;
            k_idx = k_idx + 1;
        end
    end
    
    % Identify boundary
    eps_tol = 1e-10;
    v = mesh.verts;
    is_boundary = (v(:,1) < eps_tol | v(:,1) > 1-eps_tol | ...
                   v(:,2) < eps_tol | v(:,2) > 1-eps_tol);
    mesh.bndry = find(is_boundary);
    mesh.bndry = mesh.bndry(:); 
    
    save(fullfile(meshDir, sprintf('quad_%02d.mat', N)), 'mesh');
    plotVEMMesh(mesh); 
    title(sprintf('VEM Quad-Mesh %dx%d', N, N));
end