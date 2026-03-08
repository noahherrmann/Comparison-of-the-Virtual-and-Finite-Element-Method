%sf_dsf_1D
% Calculate shape function and derivatives for 1D problem from p = 1 to p =
% 4 in a [-1 1] domain
function [N,DN] = sf_dsf_1D(p,xi)
switch p
    case 1
        % Shape Functions
        N(1) = (1-xi)/2;
        N(2) = (1+xi)/2;
        % Derivative of Shape Functions
        DN(1) = -1/2;
        DN(2) = 1/2;
    case 2
        % Shape Functions
        N(1) = xi/2*(xi-1);
        N(2) = (1+xi)*(1-xi);
        N(3) = xi/2*(xi+1);
        % Derivative of Shape Functions
        DN(1) = xi-1/2;
        DN(2) = -2*xi;
        DN(3) = xi+1/2;
    case 3
        % Shape Functions
        N(1) = -1/16 + 1/16*xi + 9/16*xi^2 - 9/16*xi^3;
        N(2) = 9/16 - 27/16*xi - 9/16*xi^2 + 27/16*xi^3;
        N(3) = 9/16 + 27/16*xi - 9/16*xi^2 - 27/16*xi^3;
        N(4) = -1/16 - 1/16*xi + 9/16*xi^2 + 9/16*xi^3;
        % Derivative of Shape Functions
        DN(1) = 1/16 + 9/8*xi - 27/16*xi^2;
        DN(2) = -27/16 - 9/8*xi + 81/16*xi^2;
        DN(3) = 27/16 - 9/8*xi - 81/16*xi^2;
        DN(4) = -1/16 + 9/8*xi + 27/16*xi^2;
    case 4
        % Shape Functions
        N(1) = 1/6*xi - 1/6*xi^2 - 2/3*xi^3 + 2/3*xi^4;
        N(2) = -4/3*xi + 8/3*xi^2 + 4/3*xi^3 - 8/3*xi^4;
        N(3) = 1 - 5*xi^2 + 4*xi^4;
        N(4) = 4/3*xi + 8/3*xi^2 - 4/3*xi^3 - 8/3*xi^4;
        N(5) = -1/6*xi - 1/6*xi^2 + 2/3*xi^3 + 2/3*xi^4;
        % Derivative of Shape Functions
        DN(1) = 1/6 - 1/3*xi - 2*xi^2 + 8/3*xi^3;
        DN(2) = -4/3 + 16/3*xi + 4*xi^2 -32/3*xi^3;
        DN(3) = -10*xi + 16*xi^3;
        DN(4) = 4/3 + 16/3*xi - 4*xi^2 -32/3*xi^3;
        DN(5) = -1/6 - 1/3*xi + 2*xi^2 + 8/3*xi^3;
end
end