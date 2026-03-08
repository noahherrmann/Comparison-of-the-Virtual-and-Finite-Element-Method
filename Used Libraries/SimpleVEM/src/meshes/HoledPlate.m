% DESCRIPTION:
% This function returns the structure "geom", composed by:
% 1) coord: Matrix that contain the coordinations of each point of the
%           mesh (every raw of this matrix define a point)
% 2) conn: Matrix that contain the indicies of each point that compose an
%          element (every raw difine an element)
% The function also return the contour of the figure. 
% The Mesh proposed is formed by only qudrilateral elements.
% 
% PARAMETERS:
% x: length 
% y: width
% r: hole radius
% nx: number of nodes on x axis in the second part of the plate
% nth: number of nodes in theta direction
% nr: number of nodes in radius direction
% nn: number of points to discretize the hole contour
% 
% NOTES:
% It's very important to define "nth" as an odd number (nunmero dispari), otherwise the
% discretization won't be what you expect

function [geom, contour, options] = HoledPlate(x,y,r,nx,nr,nth,nn)

th = linspace(pi/2,0,nth);
r_vect = linspace(r,y,nr);
x_rett = linspace(y,x,nx);
y_rett = zeros(1,round(nth/2));

n = nth*nr+(nx-1)*round(nth/2);
n1 = nth*nr;
n_el = (nth-1)*(nr-1)+(round(nth/2)-1)*(nx-1);
n_el1 = (nth-1)*(nr-1);
n_el2 = (round(nth/2)-1)*(nx-1);

Points_coord = zeros(n,2);
Points_conn = zeros(n_el,5);
Points_conn(:,1) = 4;

for i = 1:nr-1
    for j = 1:nth
        Points_coord((i-1)*nth+j,:) = [r_vect(i)*cos(th(j)),r_vect(i)*sin(th(j))];
    end
end

j = 1;
for i = 1:nth
    if i < round(nth/2)
        Points_coord((nr-1)*nth+i,:) = [y*tan(pi/2-th(i)),y];
    else 
        Points_coord((nr-1)*nth+i,:) = [y,y*tan(th(i))];
        y_rett(j) = y*tan(th(i));
        j = j+1;
    end
end

for i = 1:nx-1
    for j =1:round(nth/2)
        Points_coord((i-1)*round(nth/2)+n1+j,:) = [x_rett(i+1),y_rett(j)];
    end
end

j = 1;
for i = 1:n_el1
    if mod(j,nth) == 0
        j = j+1;
    end
    Points_conn(i,2:5) = [j,j+1,j+nth+1,j+nth];
    j = j+1;
end

k = 1;
ind = j+round(nth/2);
for i = 1:n_el2
    if mod(k,round(nth/2)) == 0
        ind = ind+1;
        k = k+1;
    end
    Points_conn(n_el1+i,2:5) = [ind,ind+1,ind+round(nth/2)+1,ind+round(nth/2)];
    ind = ind+1;
    k = k+1;
end

geom.coord = Points_coord;
geom.conn = Points_conn;

contour = zeros(nn+3,2);
teta = linspace(0,deg2rad(90),nn);
contour(1:nn,1) = r*cos(teta);
contour(1:nn,2) = r*sin(teta);
contour(nn+1,:) = [0, y];
contour(nn+2,:) = [x, y];
contour(nn+3,:) = [x, 0];

options.problemtype = 6;
options.data = [x, y, r];
end
