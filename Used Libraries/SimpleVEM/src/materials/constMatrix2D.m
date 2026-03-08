function [D_out] = constMatrix2D(matprop)

ntype = matprop.ntype;
young = matprop.young;
poiss = matprop.poiss;

if  ntype == 1                   % plane stress
    const = young/(1-poiss^2);
    D_out = const*[1 poiss 0 ; poiss 1 0 ; 0 0 (1-poiss)/2 ];
elseif  ntype == 2               % Plane strain
    matprop.thick = 1;
    const = young/((1+poiss)*(1-2*poiss));
    D_out = const*[1-poiss poiss 0 ; poiss 1-poiss 0 ; 0 0 (1-2*poiss)/2 ];
end

end