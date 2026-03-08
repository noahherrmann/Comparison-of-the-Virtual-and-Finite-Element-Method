%% VEM Fire Ball Time-lapse
clear; close all; clc;
load('./mesh/voronoi_03.mat'); 
f_handle = @(x,y) 0;
[K, ~, ~, ~, areas] = VEM(mesh, f_handle); 
ndof = size(mesh.verts, 1);

% We only have the values on the vertices.
% All we can do is distribute the area evenly on the vertices
M_diag = zeros(ndof, 1);
for el = 1:length(mesh.elems)
    v_ids = mesh.elems{el};
    n_sides = length(v_ids);
    M_diag(v_ids) = M_diag(v_ids) + areas(el) / n_sides;
end
M = sparse(1:ndof, 1:ndof, M_diag);

dt = 0.01;
T_final = 1.0;
steps = round(T_final/dt);

% Hot starting point
u = zeros(ndof, 1);
dist = sqrt((mesh.verts(:,1) - 0.5).^2 + (mesh.verts(:,2) - 0.5).^2);
u(dist < 0.2) = 100; 

% Cold Boundary set < 0.01 for other direction
isBdNode = (abs(mesh.verts(:,1)) > 0.99 | abs(mesh.verts(:,2)) > 0.99);
freeNodes = find(~isBdNode);

figure('Color', 'w');
h = patch('Faces', [], 'Vertices', mesh.verts, 'FaceVertexCData', u, ...
          'FaceColor', 'interp', 'EdgeColor', 'none');
colormap(hot); colorbar; clim([0 100]); axis equal; axis([-1 1 -1 1]);
title('VEM Heat Propagation');

for n = 1:steps
    % Implicit Euler:
    lhs = M + dt*K;
    rhs = M*u;
    
    u_next = zeros(ndof, 1);
    u_next(freeNodes) = lhs(freeNodes, freeNodes) \ rhs(freeNodes);
    u = u_next;
    
    % Triangluation for easier plots
    if n == 1
        tri = [];
        for el = 1:length(mesh.elems)
            v = mesh.elems{el};
            for j = 2:length(v)-1
                tri = [tri; v(1), v(j), v(j+1)];
            end
        end
    end
    
    trisurf(tri, mesh.verts(:,1), mesh.verts(:,2), u, 'EdgeColor', 'none');
    view(2); drawnow;
end

function [K,F,PP,stab,areas,diameters,out] = VEM(mesh,f,varargin)
% This code is a modification of the code introduced in 
% Sutton, O.J., "The virtual element method in 50 lines of Matlab", 2016.
% (c) Oliver J. Sutton, University of Leicester, 2016.
if(size(varargin,2)==1)
    coef = varargin{1};
else
    coef = ones(size(mesh.elems,1),1);
end
ndof = size(mesh.verts, 1); n_polys = 3; % Method has 1 degree of freedom per vertex
% vectors to store sparse matrix
numStore = cellfun(@numel,mesh.elems);
numStore = sum(numStore.^2);
proj = cell(size(mesh.elems,1),1); % projections for error estimates
ii = nan(numStore,1); jj = nan(numStore,1); ss = nan(numStore,1); pt = 0;
s2 = nan(numStore,1);
F  = zeros(ndof, 1); % Forcing vector
linear_polynomials = {[0,0], [1,0], [0,1]}; % Impose an ordering on the linear polynomials
mod_wrap = @(x, a) mod(x-1, a) + 1; % A utility function for wrapping around a vector
PP = cell(length(mesh.elems),1);
areas = nan(length(mesh.elems),1);
diameters = nan(length(mesh.elems),1);
for el_id = 1:length(mesh.elems) %
    vert_ids = mesh.elems{el_id}; % Global IDs of the vertices of this element
    verts = mesh.verts(vert_ids, :); % Coordinates of the vertices of this element
    n_sides = length(vert_ids); % Start computing the geometric information
    area_components = verts(:,1) .* verts([2:end,1],2) - verts([2:end,1],1) .* verts(:,2);
    area = 0.5 * abs(sum(area_components));
    areas(el_id) = area;
    centroid = sum((verts + verts([2:end,1],:)) .* repmat(area_components,1,2))/(6*area);
    diameter = 0; % Compute the diameter by looking at every pair of vertices
    for i = 1:(n_sides-1)
        for j = (i+1):n_sides
            diameter = max(diameter, norm(verts(i,:)-verts(j,:)));
        end
    end
    diameters(el_id) = diameter;
    D = zeros(n_sides, n_polys); D(:, 1) = 1;
    B = zeros(n_polys, n_sides); B(1, :) = 1/n_sides;
    for vertex_id = 1:n_sides
        vert = verts(vertex_id, :); % This vertex and its neighbours
        prev = verts(mod_wrap(vertex_id - 1, n_sides), :);
        next = verts(mod_wrap(vertex_id + 1, n_sides), :);
        vertex_normal = [next(2) - prev(2), prev(1) - next(1)]; % Average of normals on edges
        for poly_id = 2:n_polys % Only need to loop over non-constant polynomials
            poly_degree = linear_polynomials{poly_id};
            monomial_grad = poly_degree / diameter; % Gradient of a linear polynomial is constant
            D(vertex_id, poly_id) = dot((vert - centroid), poly_degree) / diameter;
            B(poly_id, vertex_id) = 0.5 * dot(monomial_grad, vertex_normal);
        end
    end
    projector = (B*D)\B; % Compute the local Ritz projector to polynomials
    proj{el_id} = projector;
    stabilising_term = (eye(n_sides) - D * projector)' * (eye(n_sides) - D * projector);
    G = B*D; G(1, :) = 0;
    if(isa(coef,'function_handle'))
        rhoV = coef(centroid(1),centroid(2));
    else
        rhoV = coef(el_id);
    end
    PP{el_id} = projector;
    local_stiffness = (projector' * G * projector + stabilising_term)*rhoV;
    [xx,yy] = meshgrid(vert_ids);
    adv = numel(local_stiffness);
    ii(pt+1:pt+adv) = xx(:);
    jj(pt+1:pt+adv) = yy(:);
    ss(pt+1:pt+adv) = local_stiffness(:);
    s2(pt+1:pt+adv) = stabilising_term(:);
    pt = pt + adv;
    F(vert_ids) = F(vert_ids) + f(centroid(1),centroid(2)) * area/n_sides;
end
K = sparse(ii,jj,ss,ndof,ndof);
stab = sparse(ii,jj,s2,ndof,ndof);

% save struct with info for L2-error
out = struct('mesh',mesh,'proj',{proj},'polys',linear_polynomials);
end