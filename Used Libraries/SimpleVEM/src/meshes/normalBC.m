function [diric, neum_c, neum_d] = normalBC(coord, conn, options)
% normalBC applies boundary conditions (both Dirichlet and Neumann)
% for seven types of problems:
%   1. Tension bar
%   2. Cantilever beam with concentrated load at the free end
%   3. Cantilever beam with distributed load
%   4. Simply supported beam with distributed load
%   5. Tube under internal and external pressure
%   6. Perforated plate under tension
%   7. Cook's membrane (default case)
%
% It also displays the BCs on the current figure. Dirichlet nodes are shown
% as filled red markers. Distributed Neumann BCs are shown by drawing a red 
% line over the edge where the load is applied.
%
% INPUTS:
%   coord   - Nodal coordinates (Nx2 array)
%   conn    - Connectivity matrix (first column: number of nodes in element,
%             then the node numbers)
%   options - Structure containing:
%               options.problemtype : integer in {1,...,7}
%               options.data        : array with problem data (e.g., dimensions)
%               options.meshType    : (optional) indicator for mesh type (used in case 5)
%
% OUTPUTS:
%   diric   - Dirichlet BC matrix [node, BC identifier, prescribed ux, prescribed uy]
%   neum_c  - Neumann BC matrix applied at nodes [node, prescribed Fx, prescribed Fy]
%   neum_d  - Neumann BC matrix applied on edges [element, side index, Fx, Fy]

% Small offset for numerical robustness
ds = 1e-6;

% List of node indices
listnode = (1:size(coord,1))';

% Initialize outputs
diric   = zeros(length(listnode), 4); % [node, bcID, ux, uy]
neum_c  = [];  % For concentrated Neumann BC (per node)
neum_d  = [];  % For distributed Neumann BC (per edge)

% Determine the problem type
ind = options.problemtype;

switch ind
    case 1  % Tension bar
        % Problem data: length and height
        l = options.data(1);
        h = options.data(2);
        
        % DIRICHLET BC: left side (fixed)
        polygon = [-ds, -ds; ds, -ds; ds, h+ds; -ds, h+ds];
        in = findInside(polygon(:,1), polygon(:,2), coord);
        diric(:,1) = listnode;
        % Use an identifier (here 11) to indicate Dirichlet nodes
        diric(in==1,2) = 11;
        diric(in==1,3:4) = 0;  % Prescribed displacements
        
        % NEUMANN BC: right side (traction)
        polygon = [l-ds, -ds; l+ds, -ds; l+ds, h+ds; l-ds, h+ds];
        selected_edges_info = [];
        for i = 1:size(conn, 1)
            nedgesEl = conn(i,1);
            connEl = conn(i,2:nedgesEl+1);
            sides = edgeExtract(nedgesEl, connEl);
            for j = 1:size(sides,1)
                pt1 = coord(sides(j,1),:);
                pt2 = coord(sides(j,2),:);
                ptMid = (pt1 + pt2) / 2;
                t = 0.5;
                repPoint = (1-t)^2 * pt1 + 2*(1-t)*t * ptMid + t^2 * pt2;
                if inpolygon(repPoint(1), repPoint(2), polygon(:,1), polygon(:,2))
                    selected_edges_info = [selected_edges_info; i, j];
                end
            end
        end
        % Define Neumann BC arrow vector (e.g., traction value)
        value_x = 1000;
        value_y = 0;
        for i = 1:size(selected_edges_info,1)
            element_number = selected_edges_info(i,1);
            side_number    = selected_edges_info(i,2);
            neum_d = [neum_d; element_number, side_number, value_x, value_y];
        end
        
    case 2  % Cantilever beam with concentrated load at free end
        l = options.data(1);
        h = options.data(2);
        
        % DIRICHLET BC: fixed support at left side
        polygon = [-ds, -ds; ds, -ds; ds, h+ds; -ds, h+ds];
        in = findInside(polygon(:,1), polygon(:,2), coord);
        diric(:,1) = listnode;
        diric(in==1,2) = 11;
        diric(in==1,3:4) = 0;
        
        % NEUMANN BC: concentrated load at free end (applied at nodes)
        neum_c = zeros(length(listnode), 3);
        neum_c(:,1) = listnode;
        polygon = [l-ds, h-ds; l+ds, h-ds; l+ds, h+ds; l-ds, h+ds];
        in = findInside(polygon(:,1), polygon(:,2), coord);
        value_x = 0;
        value_y = -1000;
        neum_c(in==1,2) = value_x;
        neum_c(in==1,3) = value_y;
        
    case 3  % Cantilever beam with distributed load
        l = options.data(1);
        h = options.data(2);
        
        % DIRICHLET BC: fixed support at left side
        polygon = [-ds, -ds; ds, -ds; ds, h+ds; -ds, h+ds];
        in = findInside(polygon(:,1), polygon(:,2), coord);
        diric(:,1) = listnode;
        diric(in==1,2) = 11;
        diric(in==1,3:4) = 0;
        
        % NEUMANN BC: distributed load on right edge
        polygon = [-ds, h-ds; l+ds, h-ds; l+ds, h+ds; -ds, h+ds];
        selected_edges_info = [];
        for i = 1:size(conn, 1)
            nedgesEl = conn(i,1);
            connEl = conn(i,2:nedgesEl+1);
            sides = edgeExtract(nedgesEl, connEl);
            for j = 1:size(sides,1)
                pt1 = coord(sides(j,1),:);
                pt2 = coord(sides(j,2),:);
                ptMid = (pt1 + pt2) / 2;
                t = 0.5;
                repPoint = (1-t)^2 * pt1 + 2*(1-t)*t * ptMid + t^2 * pt2;
                if inpolygon(repPoint(1), repPoint(2), polygon(:,1), polygon(:,2))
                    selected_edges_info = [selected_edges_info; i, j];
                end
            end
        end
        value_x = -10;
        value_y = 0;
        for i = 1:size(selected_edges_info,1)
            element_number = selected_edges_info(i,1);
            side_number    = selected_edges_info(i,2);
            neum_d = [neum_d; element_number, side_number, value_x, value_y];
        end
        
    case 4  % Simply supported beam with distributed load
        l = options.data(1);
        h = options.data(2);
        
        % DIRICHLET BC: first support (left)
        % Simply supported
        polygon = [-ds, -ds; ds, -ds; ds, +ds; -ds, +ds];
        % Double Fix
        % polygon = [-ds, -ds; ds, -ds; ds, 2+ds; -ds, 2+ds];
        in = findInside(polygon(:,1), polygon(:,2), coord);
        diric(:,1) = listnode;
        diric(in==1,2) = 11;
        diric(in==1,3:4) = 0;
        % DIRICHLET BC: second support (right)
        % Simply supported:
        polygon = [15-ds, -ds; 15+ds, -ds; 15+ds, +ds; 15-ds, +ds];
        % Double Fix
        % polygon = [l-ds, -ds; l+ds, -ds; l+ds, 2+ds; l-ds, 2+ds];
        in = findInside(polygon(:,1), polygon(:,2), coord);
        diric(in==1,2) = 11;
        diric(in==1,3:4) = 0;
        
        % NEUMANN BC: distributed load on the top edge
        polygon = [-ds, h-ds; l+ds, h-ds; l+ds, h+ds; -ds, h+ds];
        selected_edges_info = [];
        for i = 1:size(conn, 1)
            nedgesEl = conn(i,1);
            connEl = conn(i,2:nedgesEl+1);
            sides = edgeExtract(nedgesEl, connEl);
            for j = 1:size(sides,1)
                pt1 = coord(sides(j,1),:);
                pt2 = coord(sides(j,2),:);
                ptMid = (pt1 + pt2) / 2;
                t = 0.5;
                repPoint = (1-t)^2 * pt1 + 2*(1-t)*t * ptMid + t^2 * pt2;
                if inpolygon(repPoint(1), repPoint(2), polygon(:,1), polygon(:,2))
                    selected_edges_info = [selected_edges_info; i, j];
                end
            end
        end
        value_x = -100;
        value_y = 0;
        for i = 1:size(selected_edges_info,1)
            element_number = selected_edges_info(i,1);
            side_number    = selected_edges_info(i,2);
            neum_d = [neum_d; element_number, side_number, value_x, value_y];
        end
        
    case 5  % Tube under internal and external pressure
        % Check mesh type (Voronoi or not)
        if isempty(options.meshType)
            rr = input('Is the mesh generated with Voronoi? (1-yes, 2-no):\n');
        else
            rr = options.meshType;
        end
        if rr ~= 1
            ds = 1e-2;
        else
            ds = 5e-3;
        end
        r_int = options.data(1);
        r_ext = options.data(2);
        
        % DIRICHLET BC: inner boundary (first support)
        polygon = [-ds, r_int-ds; ds, r_int-ds; ds, r_ext+ds; -ds, r_ext+ds];
        in = findInside(polygon(:,1), polygon(:,2), coord);
        diric(:,1) = listnode;
        diric(in==1,2) = 10;
        diric(in==1,3:4) = 0;
        % DIRICHLET BC: lower boundary (second support)
        polygon = [r_int-ds, -ds; r_ext+ds, -ds; r_ext+ds, ds; r_int-ds, ds];
        in = findInside(polygon(:,1), polygon(:,2), coord);
        diric(in==1,2) = 1;
        diric(in==1,3:4) = 0;
        
        % NEUMANN BC: distributed pressure on curved boundaries
        selected_edges_info = [];
        % First curved segment
        th = deg2rad(linspace(-1,91,100))';
        polygon = [(r_int/2)*cos(th), (r_int/2)*sin(th); (r_int+ds)*cos(flip(th)), (r_int+ds)*sin(flip(th))];
        for i = 1:size(conn,1)
            nedgesEl = conn(i,1);
            connEl = conn(i,2:nedgesEl+1);
            sides = edgeExtract(nedgesEl, connEl);
            for j = 1:size(sides,1)
                pt1 = coord(sides(j,1),:);
                pt2 = coord(sides(j,2),:);
                ptMid = (pt1 + pt2) / 2;
                t = 0.5;
                repPoint = (1-t)^2 * pt1 + 2*(1-t)*t * ptMid + t^2 * pt2;
                if inpolygon(repPoint(1), repPoint(2), polygon(:,1), polygon(:,2))
                    selected_edges_info = [selected_edges_info; i, j];
                end
            end
        end
        % Second curved segment
        % polygon = [(r_ext-ds)*cos(th), (r_ext-ds)*sin(th); (r_ext+ds)*cos(flip(th)), (r_ext+ds)*sin(flip(th))];
        % for i = 1:size(conn,1)
        %     nedgesEl = conn(i,1);
        %     connEl = conn(i,2:nedgesEl+1);
        %     sides = edgeExtract(nedgesEl, connEl);
        %     for j = 1:size(sides,1)
        %         pt1 = coord(sides(j,1),:);
        %         pt2 = coord(sides(j,2),:);
        %         ptMid = (pt1 + pt2) / 2;
        %         t = 0.5;
        %         repPoint = (1-t)^2 * pt1 + 2*(1-t)*t * ptMid + t^2 * pt2;
        %         if inpolygon(repPoint(1), repPoint(2), polygon(:,1), polygon(:,2))
        %             selected_edges_info = [selected_edges_info; i, j];
        %         end
        %     end
        % end
        value_x = -10000;
        value_y = 0;
        for i = 1:size(selected_edges_info,1)
            element_number = selected_edges_info(i,1);
            side_number    = selected_edges_info(i,2);
            neum_d = [neum_d; element_number, side_number, value_x, value_y];
        end
        
    case 6  % Perforated plate under tension
        l = options.data(1);
        h = options.data(2);
        r = options.data(3);
        
        % DIRICHLET BC: one side support
        polygon = [-ds, r-ds; ds, r-ds; ds, h+ds; -ds, h+ds];
        in = findInside(polygon(:,1), polygon(:,2), coord);
        diric(:,1) = listnode;
        diric(in==1,2) = 10;
        diric(in==1,3:4) = 0;
        % DIRICHLET BC: other side support
        polygon = [r-ds, -ds; l+ds, -ds; l+ds, ds; r-ds, ds];
        in = findInside(polygon(:,1), polygon(:,2), coord);
        diric(in==1,2) = 1;
        diric(in==1,3:4) = 0;
        
        % NEUMANN BC: distributed load on the opposite side
        polygon = [l-ds, -ds; l+ds, -ds; l+ds, h+ds; l-ds, h+ds];
        selected_edges_info = [];
        for i = 1:size(conn,1)
            nedgesEl = conn(i,1);
            connEl = conn(i,2:nedgesEl+1);
            sides = edgeExtract(nedgesEl, connEl);
            for j = 1:size(sides,1)
                pt1 = coord(sides(j,1),:);
                pt2 = coord(sides(j,2),:);
                ptMid = (pt1 + pt2) / 2;
                t = 0.5;
                repPoint = (1-t)^2 * pt1 + 2*(1-t)*t * ptMid + t^2 * pt2;
                if inpolygon(repPoint(1), repPoint(2), polygon(:,1), polygon(:,2))
                    selected_edges_info = [selected_edges_info; i, j];
                end
            end
        end
        value_x = 10000;
        value_y = 0;
        for i = 1:size(selected_edges_info,1)
            element_number = selected_edges_info(i,1);
            side_number    = selected_edges_info(i,2);
            neum_d = [neum_d; element_number, side_number, value_x, value_y];
        end
        
    otherwise  % Cook's membrane (default case)
        H1 = options.data(1);
        H2 = options.data(2);
        H3 = options.data(3);
        L_val = options.data(4);
        
        % DIRICHLET BC: lower boundary
        polygon = [-ds, -ds; ds, -ds; ds, H1+ds; -ds, H1+ds];
        in = findInside(polygon(:,1), polygon(:,2), coord);
        diric(:,1) = listnode;
        diric(in==1,2) = 11;
        diric(in==1,3:4) = 0;
        
        % NEUMANN BC: applied on the top boundary
        polygon = [L_val-ds, (H1+H2-H3)-ds; L_val+ds, (H1+H2-H3)-ds; L_val+ds, (H1+H2)+ds; L_val-ds, (H1+H2)+ds];
        selected_edges_info = [];
        for i = 1:size(conn,1)
            nedgesEl = conn(i,1);
            connEl = conn(i,2:nedgesEl+1);
            sides = edgeExtract(nedgesEl, connEl);
            for j = 1:size(sides,1)
                pt1 = coord(sides(j,1),:);
                pt2 = coord(sides(j,2),:);
                ptMid = (pt1 + pt2) / 2;
                t = 0.5;
                repPoint = (1-t)^2 * pt1 + 2*(1-t)*t * ptMid + t^2 * pt2;
                if inpolygon(repPoint(1), repPoint(2), polygon(:,1), polygon(:,2))
                    selected_edges_info = [selected_edges_info; i, j];
                end
            end
        end
        value_x = 0;
        value_y = 1000;
        for i = 1:size(selected_edges_info,1)
            element_number = selected_edges_info(i,1);
            side_number    = selected_edges_info(i,2);
            neum_d = [neum_d; element_number, side_number, value_x, value_y];
        end
end

%% --- Plotting of Boundary Conditions ---
hold on;
% Plot Dirichlet nodes (colored markers)
dirichletNodes = find(diric(:,2) ~= 0);
if ~isempty(dirichletNodes)
    scatter(coord(dirichletNodes,1), coord(dirichletNodes,2), 50, 'r', 'filled');
end

% For concentrated Neumann BCs (neum_c), you might choose to display them
% with a different marker if desired. (Here, we do not plot them as red lines.)
if ~isempty(neum_c)
    for i = 1:size(neum_c,1)
        if neum_c(i,2) ~= 0 || neum_c(i,3) ~= 0
            % For example, mark the node with a blue circle.
            plot(coord(neum_c(i,1),1), coord(neum_c(i,1),2), 'bo', 'MarkerSize', 8, 'LineWidth', 2);
        end
    end
end

% Plot Neumann BCs applied on edges (neum_d) as red lines.
% For each entry, retrieve the edge endpoints and plot a red line over the edge.
if ~isempty(neum_d)
    for i = 1:size(neum_d,1)
        elem = neum_d(i,1);
        sideIdx = neum_d(i,2);
        nedgesEl = conn(elem,1);
        connEl = conn(elem,2:nedgesEl+1);
        sides = edgeExtract(nedgesEl, connEl);
        pt1 = coord(sides(sideIdx,1),:);
        pt2 = coord(sides(sideIdx,2),:);
        plot([pt1(1) pt2(1)], [pt1(2) pt2(2)], 'r-', 'LineWidth', 2);
    end
end

end
