function Kloc = localK_VEM(coord, connEl, mat, nedges, normals, edges, ndof,tau)

DEBUG = 0;

C = mat.D;

Np = eye(3);

% Delaunay Triangulation - !!! NON FUNZIONA CON POLIGONI CONVESSI !!!
DT = delaunayTriangulation(coord(connEl,:));
triangles = connEl(DT.ConnectivityList);

% Nuova funzione per poligoni convessi:
% [~, triangles] = triangulateNonConvexPolygon(coord, connEl);

trianglesCheck = [zeros(size(triangles,1),1)+3,triangles];

if DEBUG
    checkOrientation(coord, trianglesCheck)
end

% Plot Triangulation
if DEBUG
    hold on;
    patch('Faces',triangles,'Vertices',coord,...
        'FaceColor', [0.8 0.8 0.8],'EdgeColor','k', 'LineStyle', '-');
    axis equal;
    plot(coord(triangles,1),coord(triangles,2),'ks');
end

% Number of triangles on which integrate
nTri = size(triangles,1);

% Calculate G
G = zeros(3,3);
ngpsTri = 1;
% Loop over triangles
for i=1:nTri
    % Loop over gauss points
    [pt,w] = GLTrisQuad(coord(triangles(i,:),:), ngpsTri);

    if DEBUG
        hold on;
        plot(pt(:,1),pt(:,2),'xr')
    end

    for j=1:ngpsTri
        G = G + Np'*Np*w(j);
    end
end

if det(G) < 1e-10
    warning('La matrice G è singolare')
end

% Calculate Gt
Gt = zeros(3,3);
% Loop over triangles
for i=1:nTri
    % Loop over gauss points
    [pt,w] = GLTrisQuad(coord(triangles(i,:),:), ngpsTri);
    for j=1:ngpsTri
        Gt = Gt + Np'*C*Np*w(j);
    end
end

if det(Gt) < 1e-10
    disp('La matrice Gt è singolare')
end

% Calculate B
ngpsEdge = 2;
B = zeros(3,ndof*nedges);
% Loop over edges
for i=1:nedges

    nx = normals(i,1);
    ny = normals(i,2);
    
    if isnan(nx) || isnan(ny)
        warning('Attenzione: nx o ny contiene NaN.');
    end

    Ne = [nx 0 ny ; 0 ny nx];

    w = GLLineQuad(coord(edges(i,1),:), coord(edges(i,2),:));
    xi = [-1 1];
    Nv = zeros(2,ndof*nedges); % to be generalised

    % It works for k=1
    edgesL = edgesLocal(nedges);

    % Loop over gauss points (works for k=1)
    for j=1:ngpsEdge
        N1 = (1-xi(j))/2;
        N2 = (1+xi(j))/2;

        % Matrix Assembly
        Nv(1,2*edgesL(i,1)-1) = N1;
        Nv(2,2*edgesL(i,1)) = N1;

        Nv(1,2*edgesL(i,2)-1) = N2;
        Nv(2,2*edgesL(i,2)) = N2;

        B = B + (Ne*Np)'*Nv*w(j);
    end
end

% Calculate Kc
Kc = B'*(G\Gt/G)*B;

% Calculate D
Di = zeros(ndof*nedges,6);

% to be generalised over all elements (works for k=1)
for i=1:nedges
    Di(2*i-1,1) = 1;
    Di(2*i,2) = 1;

    Di(2*i-1,3) = coord(connEl(i),1);
    Di(2*i,4) = coord(connEl(i),1);

    Di(2*i-1,5) = coord(connEl(i),2);
    Di(2*i,6) = coord(connEl(i),2);
end

% Normalization
invDtD = inv(Di'*Di);

% Stabilization parameter
% if isempty(tau) == 1
%     tau = 1;
% end

% Ks
Ks = tau*trace(Kc)*(eye(ndof*nedges)-Di*invDtD*Di');

% Finale Stiffness Matrix
Kloc = Kc + Ks;

end