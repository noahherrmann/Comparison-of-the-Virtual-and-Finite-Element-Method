%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                 %
%                 CALCULATE STIFFNESS MATRIX                      %
%                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ke, Afem] = stiff(geom, mat, ielem)

% Recover necessary quantities
coord = geom.coord;
connectivity = geom.conn;
thick = mat.thick;
ngaus = 4;
nnode = 4;
ndofn = 2;
D = mat.D;

Afem = 0;

% Recover nodes coordinates
x = coord(connectivity(ielem,:),1);
y = coord(connectivity(ielem,:),2);

% Gauss Points:
ke = 0;

[pg w] = legzo(sqrt(ngaus),-1,1);
pg = flip(pg);
w = flip(w);

for i=1:length(pg)
    eta = pg(i);
    w_eta = w(i);

    for j=1:length(pg)
        xi = pg(j);
        w_xi = w(j);

        % Jacobian:
        J = jacFEM(xi, eta, x, y);

        % Global B matrix (Bxy)
        Bxy = BmatrixFEM(xi, eta, J, nnode, ndofn);

        % Area
        Afem = Afem + det(J)*w_xi*w_eta;

        % Element stiffness
        ke = ke + (Bxy'*D*Bxy)*thick*det(J)*w_xi*w_eta;
        
    end
end
% DNDx = [Bxy(1,1) Bxy(1,3) Bxy(1,5) Bxy(1,7)];
% DNDy = [Bxy(2,2) Bxy(2,4) Bxy(2,6) Bxy(2,8)];
% quiver(x,y,DNDx',DNDy', 0, 'r')
end
