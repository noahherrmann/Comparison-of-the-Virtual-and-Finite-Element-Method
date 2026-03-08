%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                 %
%                       CALCULATE AREA FEM                        %
%                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Afem] = area_FEM(coord, conn, ielem)

% Recover necessary quantities
ngaus = 4;
nnode = 4;
ndofn = 2;

% Recover nodes coordinates
x = coord(conn(ielem,:),1);
y = coord(conn(ielem,:),2);

Afem = 0;

% Recover nodes coordinates
x = coord(conn(ielem,:),1);
y = coord(conn(ielem,:),2);

% Gauss Points:
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
    end
end

end
