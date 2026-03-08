function [allTriangles] = VEM2FEM_delaunayTriangulation(geom_VEM)
    
coord = geom_VEM.coord;
conn = geom_VEM.conn;

% Numero di elementi
nel = size(conn, 1);

% Array per memorizzare tutte le triangolazioni e i rispettivi valori di UU
allTriangles = [];

% Loop sugli elementi
for ielem = 1:nel
    nedgesEl = conn(ielem, 1); % Numero di nodi dell'elemento
    connEl = conn(ielem, 2:nedgesEl+1);  % Indici dei nodi dell'elemento

    % Esegui la triangolazione di Delaunay per l'elemento
    DT = delaunayTriangulation(coord(connEl, :));
    triangles = connEl(DT.ConnectivityList);  % Ottieni i triangoli

    % Memorizza tutti i triangoli
    allTriangles = [allTriangles; triangles]; %#ok<AGROW>
end

allTriangles = [allTriangles, allTriangles(:,1)];

end