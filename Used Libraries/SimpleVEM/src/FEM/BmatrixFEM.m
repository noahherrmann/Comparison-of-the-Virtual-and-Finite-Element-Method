function Bxy = BmatrixFEM(xi, eta, J, nnode, ndofn)
DNNxi = s2424_DNNx([-1 -1 ; 1 -1 ; 1 1 ; -1 1],[xi,eta]);
DNNeta = s2424_DNNy([-1 -1 ; 1 -1 ; 1 1 ; -1 1],[xi,eta]);

% POSSIBILE ERRORE GRAVE!
for inode=1:nnode
    BB(:,inode) = J'\[DNNxi(inode) ; DNNeta(inode)];
end

% Bi = [DNNxi ; DNNeta];
% BB = J\Bi;

% Assemble of B matrix
c1 = 0;
c2 = 0;
for j=1:nnode*ndofn
    if mod(j,2) ~= 0
        c1 = c1 + 1;
        Bxy(1,j) = BB(1,c1);
        Bxy(3,j) = BB(2,c1);
    else
        c2 = c2 + 1;
        Bxy(2,j) = BB(2,c2);
        Bxy(3,j) = BB(1,c2);
    end
end
end