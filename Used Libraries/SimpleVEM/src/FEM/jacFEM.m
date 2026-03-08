function J = jacFEM(xi, eta, x, y)
J = zeros(2,2);
DNNxi = s2424_DNNx([-1 -1 ; 1 -1 ; 1 1 ; -1 1],[xi,eta]);
DNNeta = s2424_DNNy([-1 -1 ; 1 -1 ; 1 1 ; -1 1],[xi,eta]);
% Jacobian elements:
J(1,1) = DNNxi*x;
J(1,2) = DNNeta*x;
J(2,1) = DNNxi*y;
J(2,2) = DNNeta*y;
end