function [NN,DNNx,DNNy] = LinearShapeFunction(pts,evaluatePts)
% Nn = 1*C1 + x*C2 + y*C3 + x*y*C4
% evaluatePts are usually the GaussPts

p1 = pts(1,:);
p2 = pts(2,:);
p3 = pts(3,:);
p4 = pts(4,:);

n = [ 1 p1(1) p1(2) p1(1)*p1(2) ; ...
    1 p2(1) p2(2) p2(1)*p2(2) ; ...
    1 p3(1) p3(2) p3(1)*p3(2) ; ...
    1 p4(1) p4(2) p4(1)*p4(2) ];


x = evaluatePts(1);
y = evaluatePts(2);

if rank(n) < length(n)
    N = pinv(n);
    NN=[1 x y x*y]*N;
else
    N = inv(n);
    NN=[1 x y x*y]/n;
end

Nx = [N(2,:); N(4,:)];
Ny = [N(3,:); N(4,:)];

DNNx = [1 y]*Nx;
DNNy = [1 x]*Ny;