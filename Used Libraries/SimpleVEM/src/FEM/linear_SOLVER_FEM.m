function sol = linear_SOLVER_FEM(geom, mat, BC, p, type)
% linear_SOLVER_FEM solves the linear finite element problem using the
% provided geometry, material, and boundary condition information.
%
% INPUTS:
%   geom   - Structure containing nodal coordinates (geom.coord) and 
%            connectivity (geom.conn)
%   mat    - Material properties (structure or variable as needed)
%   BC     - Boundary condition structure containing Dirichlet (BC.diric)
%            and Neumann (BC.neum_c and BC.neum_d) data.
%   p      - Polynomial order (1 or 2; p > 2 not supported)
%   type   - Type of solver ('gauss' for direct or 'conj' for conjugate gradient)
%
% OUTPUT:
%   sol    - Structure containing computed displacements (U), reactions (R),
%            separated x and y components, total elastic energy, etc.

%% PREPROCESS
disp('--------------------------------------------------');
disp('            LINEAR FEM SOLVER STARTING            ');
disp('--------------------------------------------------');

coord = geom.coord;
connectivity = geom.conn;

% Determine number of unique nodes and degrees of freedom per node
tmp = unique(coord, 'row');
npoin = size(tmp, 1);
ndofn = 2;  % two degrees of freedom (x and y)
globalDof = npoin * ndofn;

% Initialize global stiffness matrix and force vector
K = sparse(globalDof, globalDof);
F = zeros(globalDof, 1);

% Initialize total area accumulator
Atot = 0;

% For 2D quadrilateral elements, number of nodes per element is 4
nnode = 4;

% Number of elements in the mesh
nelem = size(connectivity, 1);
disp(['[PREPROCESS] Total number of elements: ', num2str(nelem)]);
disp('--------------------------------------------------');

% Apply point loads (concentrated Neumann BCs)
if ~isempty(BC.neum_c)
    disp('[PREPROCESS] Applying concentrated point loads...');
    for i = 1:npoin
        idxX = (i - 1) * ndofn + 1;
        idxY = idxX + 1;
        F(idxX) = BC.neum_c(i, 2);
        F(idxY) = BC.neum_c(i, 3);
    end
end

% Apply distributed loads (Neumann BCs on edges)
if ~isempty(BC.neum_d)
    disp('[PREPROCESS] Applying distributed loads on edges...');
    for ibc = 1:size(BC.neum_d, 1)
        elem = BC.neum_d(ibc, 1);
        if p == 1
            sides = [connectivity(elem, [1, 2]); 
                     connectivity(elem, [2, 3]); 
                     connectivity(elem, [3, 4]); 
                     connectivity(elem, [4, 1])];
        elseif p == 2
            sides = [connectivity(elem, [1, 2, 3]); 
                     connectivity(elem, [3, 5, 8]); 
                     connectivity(elem, [8, 7, 6]); 
                     connectivity(elem, [6, 4, 1])];
        else
            error('Polynomial order p > 2 is not supported!')
        end
        side_apply = sides(BC.neum_d(ibc, 2), :);
        Pn = BC.neum_d(ibc, 3);
        Pt = BC.neum_d(ibc, 4);
        Fedge = neumannAssign(coord, side_apply, Pn, Pt);
        F = F + Fedge;
    end
end

% Define fixed (prescribed) degrees of freedom
disp('[PREPROCESS] Defining fixed degrees of freedom...');
index = find(BC.diric(:, 2) ~= 0);
ipresc = BC.diric(index, :);
nfix = size(ipresc, 1);
ifdoffix = zeros(globalDof, 1);
valuedoffix = zeros(globalDof, 1);
for i_fix = 1:nfix
    ipoin = ipresc(i_fix, 1);
    bcType = ipresc(i_fix, 2);
    if bcType == 11
        ifdoffix(ipoin * ndofn - 1) = 1;
        ifdoffix(ipoin * ndofn) = 1;
        valuedoffix(ipoin * ndofn - 1) = ipresc(i_fix, 3);
        valuedoffix(ipoin * ndofn) = ipresc(i_fix, 4);
    elseif bcType == 1
        ifdoffix(ipoin * ndofn) = 1;
        valuedoffix(ipoin * ndofn) = ipresc(i_fix, 4);
    elseif bcType == 10
        ifdoffix(ipoin * ndofn - 1) = 1;
        valuedoffix(ipoin * ndofn - 1) = ipresc(i_fix, 3);
    elseif bcType == 2
        ifdoffix(ipoin * ndofn - 1) = 1;
        ifdoffix(ipoin * ndofn) = 1;
        valuedoffix(ipoin * ndofn - 1) = ipresc(i_fix, 3);
        valuedoffix(ipoin * ndofn) = ipresc(i_fix, 4);
    end
end
disp(['[PREPROCESS] Number of fixed nodes: ', num2str(nfix)]);
disp('--------------------------------------------------');

% Display total applied Neumann forces
indOdd = mod(1:length(F), 2) == 1;
indEven = mod(1:length(F), 2) == 0;
Fx_total = sum(F(indOdd));
Fy_total = sum(F(indEven));
disp(['[PREPROCESS] Total Neumann force in x: ', num2str(Fx_total)]);
disp(['[PREPROCESS] Total Neumann force in y: ', num2str(Fy_total)]);
disp('--------------------------------------------------');

%% ANALYSIS
% (Assembly and reduction steps are performed but not logged in detail)
disp('[ANALYSIS] Assembling global stiffness matrix and building reduced system...');
for ielem = 1:nelem
    % Compute element stiffness matrix and area
    [ke, A] = stiff(geom, mat, ielem);
    Atot = Atot + A;
    % Assemble element contribution into global stiffness matrix
    for inode = 1:nnode
        for idofn = 1:ndofn
            for jnode = 1:nnode
                for jdofn = 1:ndofn
                    lrow = (inode - 1) * ndofn + idofn;
                    lcol = (jnode - 1) * ndofn + jdofn;
                    globalRow = (connectivity(ielem, inode) - 1) * ndofn + idofn;
                    globalCol = (connectivity(ielem, jnode) - 1) * ndofn + jdofn;
                    K(globalRow, globalCol) = K(globalRow, globalCol) + ke(lrow, lcol);
                end
            end
        end
    end
end
disp('--------------------------------------------------');
disp(['[ANALYSIS] Total area of elements: ', num2str(Atot)]);
disp('--------------------------------------------------');

% Build reduced system (Kstar and Fstar) [summarized]
ngdof = globalDof;
freedoftable = find(ifdoffix == 0);
nfreedof = length(freedoftable);
Kstar = K(freedoftable, freedoftable);
Fstar = F(freedoftable);
% Add contributions from prescribed displacements
for igdof = 1:ngdof
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
Ustar_full = zeros(globalDof, 1);
nonFixedIdx = find(ifdoffix == 0);
fixedIdx = find(ifdoffix == 1);
Ustar_full(nonFixedIdx) = Ustar;
Ustar_full(fixedIdx) = valuedoffix(fixedIdx);
U = Ustar_full;
sol.U = U;

%% POSTPROCESS
disp('[POSTPROCESS] Computing reaction forces at fixed DOFs...');
R = zeros(globalDof, 1);
for igdof = 1:globalDof
    if ifdoffix(igdof) == 1
        R(igdof) = K(igdof, :) * U - F(igdof);
    end
end

% Separate displacements and reactions into x and y components
Ux = U(1:2:end);
Uy = U(2:2:end);
Rx = R(1:2:end);
Ry = R(2:2:end);

disp(['[POSTPROCESS] Total reaction force in x: ', num2str(sum(Rx))]);
disp(['[POSTPROCESS] Total reaction force in y: ', num2str(sum(Ry))]);
disp(['[POSTPROCESS] Maximum displacement in x: ', num2str(max(abs(Ux)))]);
disp(['[POSTPROCESS] Maximum displacement in y: ', num2str(max(abs(Uy)))]);
disp('--------------------------------------------------');

sol.Rx = Rx;
sol.Ry = Ry;
sol.Ux = Ux;
sol.Uy = Uy;
sol.UU = sqrt(Ux.^2 + Uy.^2);

disp(['[POSTPROCESS] Total elastic energy: ', num2str(-0.5 * U' * K * U)]);
sol.Energy = -0.5 * U' * K * U;

disp('--------------------------------------------------');
disp('            LINEAR FEM SOLVER FINISHED            ');
disp('--------------------------------------------------');

% Add extra spacing to separate solutions in the log
disp(' ');
disp(' ');
end
