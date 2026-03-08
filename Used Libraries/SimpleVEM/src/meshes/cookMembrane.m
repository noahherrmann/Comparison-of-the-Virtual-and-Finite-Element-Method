% DESCRIPTION:
% This function returns the structure "geom", composed by:
% 1) coord: Matrix that contain the coordinations of each point of the
%           mesh (every raw of this matrix define a point)
% 2) conn: Matrix that contain the indicies of each point that compose an
%          element (every raw difine an element)
% The Mesh proposed is formed by only qudrilateral elements.
% The function also return the contour of the figure. 
% 
% PARAMETERS:
% H1: boundary height 
% H2: distance from P and the uppest point of the boundary on y
% H3: height of the free vertical face
% L: membrane length
% nx: number of point on axis x
% ny: number of point on axis y

function [geom, contour, options] = cookMembrane(H1, H2, H3, L, nx, ny)

x = linspace(0,L,nx)';
y1 = linspace(0,H1,ny)';
y2 = linspace(H1+H2-H3,H1+H2,ny)';

coord = [];
for i = 1:ny
    coord = [coord; x,linspace(y1(i),y2(i),nx)']; %#ok<AGROW>
end

n_el = (nx-1)*(ny-1);
conn = zeros(n_el,5);
conn(:,1) = 4;
j = 1;
for i = 1:n_el
    if mod(j,nx) == 0
        j = j+1;
    end
    conn(i, 2:end) = [j, j+1, j+nx+1, j+nx];
    j = j+1;
end

geom.conn = conn;
geom.coord = coord;
contour = [0, 0; L, H1+H2-H3; L, H1+H2; 0, H1];

options.problemtype = 7;
options.data = [H1, H2, H3, L];
end