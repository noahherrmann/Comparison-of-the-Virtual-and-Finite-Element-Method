function matprop = linearElastic2D_Steel
    
% young module
matprop.young = 200000;
% poisson ratio
matprop.poiss = 0.3;
% thermic conducibility
matprop.k = 1;
% element thickness
matprop.thick = 1; % element thickness

% problem type (1 = plane stress ; 2 = plane strain)
matprop.ntype = 1;

matprop.D = constMatrix2D(matprop);
    
end