function [e] = HPComparisonUU(n,tau)

% MATERIALE 
mat = linearElastic2D_Steel;

% PIASTRA FORATA
p = 1;

[geom, ~,options] = HoledPlate(1,0.3,0.1,n(1),n(2),n(3),20);
geom_VEM = geom;
geom_FEM = geom;
geom_FEM.conn = geom_FEM.conn(:,2:end);

options.p = p;

[BC.diric, BC.neum_c, BC.neum_d] = normalBC(geom_VEM.coord, geom_VEM.conn, options);

sol_VEM = linear_SOLVER_VEM(geom_VEM, mat, BC, options.p, "gauss",tau);
sol_FEM = linear_SOLVER_FEM(geom_FEM, mat, BC, options.p, "gauss"); 

UV = sol_VEM.UU;
UF = sol_FEM.UU;
EV = sol_VEM.Energy;
EF = sol_FEM.Energy;
% e = sum((UF-UV).^2)/length(UV);
e = abs(EF-EV);

end