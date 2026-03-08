function [sol] = linear_SOLVER_VEM(geom, mat, BC, k, type, tau)
% linear_SOLVER_VEM solves the linear problem using the Virtual Element Method.
%
% INPUTS:
%   geom  - Structure with nodal coordinates (geom.coord) and connectivity (geom.conn)
%   mat   - Material properties structure or variable
%   BC    - Boundary condition structure containing:
%             BC.diric  - Dirichlet boundary conditions
%             BC.neum_c - Concentrated Neumann boundary conditions (point loads)
%             BC.neum_d - Distributed Neumann boundary conditions (applied on edges)
%   k     - Indicator for the type of element edges to use (e.g., k==1 for quadrilateral)
%   type  - Solver type: 'gauss' for direct or 'conj' for conjugate gradient solver
%   tau   - Stabilization parameter (if not provided, tau is set to 1)
%
% OUTPUT:
%   sol   - Structure containing the computed displacements (U), reactions (R),
%           separated x and y components, total elastic energy, and reduced stiffness matrix.

%% PREPROCESS
disp('--------------------------------------------------');
disp('          LINEAR VEM SOLVER STARTING              ');
disp('--------------------------------------------------');

DEBUG = 0;  % Set to 1 to plot element normals

coord = geom.coord;
conn = geom.conn;

ndof = 2;                   % Degrees of freedom per node (x and y)
npoin = size(coord, 1);     % Total number of nodes
nel = size(conn, 1);        % Total number of elements
globalDof = npoin * ndof;

disp(['[PREPROCESS] Total number of elements: ', num2str(nel)]);
disp('--------------------------------------------------');

% Initialize global stiffness matrix (sparse) for all nodes
K = sparse(globalDof, globalDof);

% Build the global force vector F
F = zeros(globalDof, 1);
% Apply concentrated (point) Neumann loads
if ~isempty(BC.neum_c)
    disp('[PREPROCESS] Applying concentrated Neumann loads...');
    for i = 1:npoin
        idxX = (i-1)*ndof + 1;
        idxY = idxX + 1;
        F(idxX) = BC.neum_c(i, 2);
        F(idxY) = BC.neum_c(i, 3);
    end
end

% Apply distributed Neumann loads on edges
if ~isempty(BC.neum_d)
    disp('[PREPROCESS] Applying distributed Neumann loads on edges...');
    for ibc = 1:size(BC.neum_d, 1)
        elem = BC.neum_d(ibc, 1);
        if k == 1
            nedgesEl_temp = conn(elem, 1);
            connEl_temp = conn(elem, 2:nedgesEl_temp+1);
            sides = edgeExtract(nedgesEl_temp, connEl_temp);
        elseif k == 2
            % Alternative element type handling can be added here
        else
            error('Polynomial order k > 2 is not supported!');
        end
        side_apply = sides(BC.neum_d(ibc, 2), :);
        Pn = BC.neum_d(ibc, 3);
        Pt = BC.neum_d(ibc, 4);
        Fedge = neumannAssign(coord, side_apply, Pn, Pt);
        F = F + Fedge;
    end
end

% Define prescribed (Dirichlet) degrees of freedom
disp('[PREPROCESS] Defining prescribed degrees of freedom...');
index = find(BC.diric(:, 2) ~= 0);
ipresc = BC.diric(index, :);
nfix = size(ipresc, 1);
ifdoffix = zeros(globalDof, 1);
valuedoffix = zeros(globalDof, 1);
for i_fix = 1:nfix
    ipoin = ipresc(i_fix, 1);
    bcType = ipresc(i_fix, 2);
    if bcType == 11
        ifdoffix(ipoin*ndof-1) = 1;
        ifdoffix(ipoin*ndof)   = 1;
        valuedoffix(ipoin*ndof-1) = ipresc(i_fix, 3);
        valuedoffix(ipoin*ndof)   = ipresc(i_fix, 4);
    elseif bcType == 1
        ifdoffix(ipoin*ndof) = 1;
        valuedoffix(ipoin*ndof) = ipresc(i_fix, 4);
    elseif bcType == 10
        ifdoffix(ipoin*ndof-1) = 1;
        valuedoffix(ipoin*ndof-1) = ipresc(i_fix, 3);
    elseif bcType == 2
        ifdoffix(ipoin*ndof-1) = 1;
        ifdoffix(ipoin*ndof) = 1;
        valuedoffix(ipoin*ndof-1) = ipresc(i_fix, 3);
        valuedoffix(ipoin*ndof) = ipresc(i_fix, 4);
    end
end
disp(['[PREPROCESS] Number of fixed nodes: ', num2str(nfix)]);
disp('--------------------------------------------------');

%% ANALYSIS
disp('[ANALYSIS] Assembling global stiffness matrix using VEM...');
for ielem = 1:nel
    % Get number of edges and element connectivity for the current element
    nedgesEl = conn(ielem, 1);
    connEl = conn(ielem, 2:nedgesEl+1);
    
    % Extract element edges and calculate outward normals
    edgesEl = edgeExtract(nedgesEl, connEl);
    normalsEl = calculateNormals(coord, edgesEl);
    
    % Optionally plot normals for debugging
    if DEBUG
        plotNormals(coord, edgesEl, normalsEl);
    end
    
    % Set tau if not provided
    if nargin < 6 || isempty(tau)
        tau = 1;
    end
    
    % Compute local stiffness matrix using VEM formulation
    ke = localK_VEM(coord, connEl, mat, nedgesEl, normalsEl, edgesEl, ndof, tau);
    
    % Assemble element contribution into global stiffness matrix K
    for inode = 1:nedgesEl
        for idof = 1:ndof
            for jnode = 1:nedgesEl
                for jdof = 1:ndof
                    lrow = (inode-1)*ndof + idof;
                    lcol = (jnode-1)*ndof + jdof;
                    globalRow = (connEl(1, inode)-1)*ndof + idof;
                    globalCol = (connEl(1, jnode)-1)*ndof + jdof;
                    K(globalRow, globalCol) = K(globalRow, globalCol) + ke(lrow, lcol);
                end
            end
        end
    end
end
disp('--------------------------------------------------');

% Build reduced system (Kstar, Fstar)
freedoftable = find(ifdoffix == 0);
nfreedof = length(freedoftable);
Kstar = K(freedoftable, freedoftable);
Fstar = F(freedoftable);
for igdof = 1:globalDof
    if ifdoffix(igdof) == 1
        Fstar = Fstar - K(freedoftable, igdof) * valuedoffix(igdof);
    end
end
disp(['[ANALYSIS] Reduced system size: ', num2str(nfreedof), ' DOFs']);
disp('--------------------------------------------------');

% Solve the reduced system
disp('----------------- SOLVING SYSTEM -----------------');
if strcmp(type, 'gauss')
    Ustar = Kstar \ Fstar;
elseif strcmp(type, 'conj')
    Ustar = pcg(Kstar, Fstar);
else
    error('Unknown solver type. Use ''gauss'' or ''conj''.');
end
disp(['[ANALYSIS] System solved using ', type, ' solver.']);
disp('--------------------------------------------------');

% Assemble full global displacement vector
U = zeros(globalDof, 1);
nonFixedIndices = find(ifdoffix == 0);
fixedIndices = find(ifdoffix == 1);
U(nonFixedIndices) = Ustar;
U(fixedIndices) = valuedoffix(fixedIndices);
sol.U = U;

%% POSTPROCESS
disp('[POSTPROCESS] Computing reaction forces at fixed DOFs...');
R = zeros(globalDof, 1);
for igdof = 1:globalDof
    if ifdoffix(igdof) == 1
        R(igdof) = K(igdof, :) * U - F(igdof);
    end
end

Ux = U(1:2:end);
Uy = U(2:2:end);
Rx = R(1:2:end);
Ry = R(2:2:end);
sol.Rx = Rx;
sol.Ry = Ry;
sol.Ux = Ux;
sol.Uy = Uy;
sol.UU = sqrt(Ux.^2 + Uy.^2);

disp(['[POSTPROCESS] Total reaction force in x: ', num2str(sum(Rx))]);
disp(['[POSTPROCESS] Total reaction force in y: ', num2str(sum(Ry))]);
disp(['[POSTPROCESS] Maximum displacement in x: ', num2str(max(abs(Ux)))]);
disp(['[POSTPROCESS] Maximum displacement in y: ', num2str(max(abs(Uy)))]);
disp(['[POSTPROCESS] Maximum overall displacement: ', num2str(max(abs(U)))]);
disp('--------------------------------------------------');

Energy = -0.5 * U' * K * U;
disp(['[POSTPROCESS] Total elastic energy: ', num2str(Energy)]);
sol.Energy = Energy;
sol.K = Kstar;

disp('--------------------------------------------------');
disp('          LINEAR VEM SOLVER FINISHED              ');
disp('--------------------------------------------------');
disp(' ');
disp(' ');
end
