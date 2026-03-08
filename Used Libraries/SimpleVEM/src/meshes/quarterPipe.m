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
% r_int: internal radius 
% r_ext: external radius
% nr: number of nodes in radius direction
% nth: number of nodes in theta direction
% nn: number of points to discretize the hole contour

function [geom, contour, options] = quarterPipe(r_int,r_ext,nr,nth,nn)

r = linspace(r_int, r_ext, nr);
th = linspace(0, pi/2, nth)';

coord = [];
for i = 1:nr
    coord = [coord; r(i)*cos(th), r(i)*sin(th)]; %#ok<AGROW>
end

n_el = (nr-1)*(nth-1);
conn = zeros(n_el, 5);
conn(:,1) = 4;
j = 1;
for i = 1:n_el
    if mod(j,nth) == 0
        j = j+1;
    end
    conn(i, 2:end) = [j, j+nth, j+nth+1, j+1];
    j = j+1;
end

geom.conn = conn;
geom.coord = coord;

teta1 = linspace(0, pi/2, nn)';
teta2 = flip(teta1);
contour = [r_int*cos(teta1),r_int*sin(teta1);r_ext*cos(teta2),r_ext*sin(teta2)];

options.problemtype = 5;
options.data = [r_int,r_ext];
end