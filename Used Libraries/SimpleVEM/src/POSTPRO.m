function POSTPRO(geom, sol)

    % Estrazione delle variabili dai dati in input
    coord = geom.coord;
    conn = geom.conn;
    U = sol.U;
    UU = sol.UU;
    Ux = sol.Ux;
    Uy = sol.Uy;

    % Calcolo della distanza massima e del fattore di amplificazione
    distances = pdist(coord, 'euclidean');
    max_distance = max(distances);
    max_U = max(abs(U));

    amp = (0.05 * max_distance) / max_U;
    % disp(['Amplification factor = ', num2str(amp)])

    % Creazione delle coordinate deformate
    coordDef = coord;
    coordDef(:,1) = coordDef(:,1) + Ux * amp;
    coordDef(:,2) = coordDef(:,2) + Uy * amp;

    % Inizializzazione della figura
    h = NewPrettyFigure;

    % Plot della mesh deformata (senza campo scalare)
    Plot_MeshFEM(geom.conn, coordDef, 1, '-', 1, h);

    % Numero di elementi
    nel = size(conn, 1);

    % Array per memorizzare tutte le triangolazioni e i rispettivi valori di UU
    allTriangles = [];
    
    % Loop sugli elementi
    for ielem = 1:nel
        nedgesEl = conn(ielem, 1); % Numero di nodi dell'elemento
        connEl = conn(ielem, 2:nedgesEl+1);  % Indici dei nodi dell'elemento

        % Esegui la triangolazione di Delaunay per l'elemento
        DT = delaunayTriangulation(coordDef(connEl, :));
        triangles = connEl(DT.ConnectivityList);  % Ottieni i triangoli
        
        % [~, triangles] = triangulateNonConvexPolygon(coord, connEl);

        % Memorizza tutti i triangoli
        allTriangles = [allTriangles; triangles]; %#ok<AGROW>
    end

    % Plot di tutti i triangoli insieme interpolando il campo scalare UU
    Plot_InterpoledColormap(coordDef, allTriangles, UU, h);

    % Colora il grafico finale con la barra colorata
    PrettifyColorbar;
    axis equal;
    axis tight;

    % Plot finale della mesh deformata sovrapposta
    Plot_MeshFEM(geom.conn, coordDef, 1, '-', 1, h);

end
