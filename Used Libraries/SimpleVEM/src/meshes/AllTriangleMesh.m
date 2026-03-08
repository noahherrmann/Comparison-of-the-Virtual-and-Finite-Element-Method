% DESCRIPTION:
% This function returns the structure "geom", composed by:
% 1) coord: Matrix that contain the coordinations of each point of the
%           mesh (every raw of this matrix define a point)
% 2) conn: Matrix that contain the indicies of each point that compose an
%          element (every raw difine an element)
% The Mesh proposed is formed by only triangular elements. Every element
% is equal to the others.
% 
% PARAMETERS:
% x: length beam
% y: width beam
% nx: number of nodes on x axis
% ny: number of nodes in y axis

function [geom, contour] = AllTriangleMesh(x,y,nx,ny)
n = nx*ny;
Points_ind = linspace(1,n,n);
x_Points = linspace(0,x,nx)';
y_Points = linspace(0,y,ny)';

Points_coord = zeros(n,2);
for i = 1:ny
    for j = 1:nx
        Points_coord((i-1)*nx+j, 1) = x_Points(j);
        Points_coord((i-1)*nx+j, 2) = y_Points(i);
    end
end
n_el = (nx-1)*(ny-1)*2;
Element_conn = zeros(n_el,4);
Element_conn(:,1) = 4;
j = 1;
k = 1;
for i = 1:n_el/2
    if mod(j,nx) == 0
        j = j+1;
    end
    Element_conn(k:k+1,2:5) = [j,j+nx,j+1,j;
                               j+1,j+nx+1,j+nx,j+1];

    k = k+2;
    j = j+1;
end

geom.coord = Points_coord;
geom.conn = Element_conn;
contour = [0, 0; x,0; x, y; 0, y];
end
