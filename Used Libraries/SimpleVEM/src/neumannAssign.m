function Fneumann = neumannAssign(coord,side,Pn,Pt)

% ngaus_edge = sqrt(elgrp.ngaus);
ndofn = 2;
p = 1;
ngaus = 2;

% gauss point defined in [-1 1]
[gp, w] = legzo(ngaus,-1,1);
gp = flip(gp);
w = flip(w);

tmp = unique(coord,'row');
npoin = size(tmp,1);
Fneumann = zeros(ndofn*npoin,1);
% Recover coordinates of Neumann nodes
x = coord(side(:),1);
y = coord(side(:),2);

pGxyz = 0;
pGxyz_n = 0;
pGxyz_t = 0;
ll = 0;
llx = 0;
lly = 0;
% loop over Gauss points
for j=1:length(gp)
    
    xi = gp(j);
    
    % Calculate shape function and derivatives
    [N, DN] = sf_dsf_1D(p,xi);
    
    % Jacobian
    % NOTA: c'é un valore assoluto perché questa "lunghezza" non
    % deve essere negativa (non avrebbe senso), il segno deve
    % essere deciso dal carico applicato (???)
    Jx = DN*x;
    Jy = DN*y;
    J = norm([Jx Jy],2);
    
    xx = N*x;
    yy = N*y;
    
    % length of the edge (DEBUG LINE)
    ll = ll + J*w(j);
    if ll<0
        disp(j)
    end
    llx = llx + Jx*w(j);
    lly = lly + Jy*w(j);
    
    % Neumann condition on nodes
%     pq_n = BV(typeOfedge).pn(xx,yy);
%     pq_t = BV(typeOfedge).pt(xx,yy);
    pn = [];
    pt = [];
    for k=1:length(x)
        pn = [pn Pn];
        pt = [pt Pt];
    end

    % Load interpolated with shape function for each edge
    pq_n = pn*N';
    pq_t = pt*N';
    
    % Normal calculation, ref [The finite element method. Vol 1 - Zienkiewicz,
    % Taylor] pag. 211
    % in this case [0 0 1] means positive traction and negative
    % compression, otherwise [0 0 -1] means the opposit
    pGxyz_n = pGxyz_n + pq_n*N*w(j).*cross([Jx Jy 0],[0 0 1])';
    % Tangential pressure positive in anticlockwise direction
    pGxyz_t = pGxyz_t + pq_t*N*w(j).*[Jx ; Jy ; 0];
    
    pGxyz = pGxyz_n + pGxyz_t;
    
end
% length of edge (DEBUG LINE)
llxy = sqrt(llx^2+lly^2);
% Assemble global vector
Fneumann(2*side(:)-1) = pGxyz(1,:);
Fneumann(2*side(:)) = pGxyz(2,:);

end