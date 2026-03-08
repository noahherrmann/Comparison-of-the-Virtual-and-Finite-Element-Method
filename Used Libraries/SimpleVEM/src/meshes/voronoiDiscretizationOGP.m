function [geom] = voronoiDiscretizationOGP(geometry, n, m)

% Creazione di un oggetto polyshape dalla matrice
geomShape = polyshape(geometry(:,1), geometry(:,2));

% Numero di punti da generare

% Generazione di punti all'interno e attorno alla geometria
padding = 1; 
minX = min(geometry(:,1)) - padding;
maxX = max(geometry(:,1)) + padding;
minY = min(geometry(:,2)) - padding;
maxY = max(geometry(:,2)) + padding;

[x, y] = meshgrid(linspace(minX,maxX,n), linspace(minY,maxY,m));  
points = [x(:), y(:)];  

% Calcolo delle celle di Voronoi
[V, C] = voronoin(points);

% Inizializzazione elementi di memorizzazione
allNodes = []; % Memorizza i nodi
elementNodeCoord = {}; % Memorizza gli indici dei nodi per ciascun elemento

% Plot delle celle di Voronoi ritagliate all'interno della geometria
% figure;
% plot(geomShape, 'FaceColor', 'none', 'EdgeColor', 'k');
% hold on;
% axis equal;

for i = 1:length(C)
    if all(C{i}~=1)  % Escludi le celle con vertici all'infinito

        % Coordinate x,y dei nodi della singola cella    
        cellNodes = V(C{i}, :); 

        % Genera una cella a partire dalle coordinate x,y dei vertici
        voronoiCell = polyshape(cellNodes(:,1), cellNodes(:,2)); 

        % Intersezione tra la cella di Voronoi e la geometria
        clippedCell = intersect(voronoiCell, geomShape); 
        
        % Salvataggio vertici delle celle ritagliate
        if ~isempty(clippedCell.Vertices)
            
            % Coordinate dei nodi di tutti gli elementi
            allNodes = [allNodes; clippedCell.Vertices]; %#ok<AGROW>
            
            % Visualizzazione della cella ritagliata
            % fill(clippedCell.Vertices(:,1), clippedCell.Vertices(:,2), rand(1,3), 'EdgeColor', 'none'); % Colora casualmente la cella

            % Ad ogni cella corrispondono le coordinate dei nodi dei singoli elementi
            elementNodeCoord{end+1} = clippedCell.Vertices; %#ok<AGROW>
        end
    end
end

% Rimozione dei nodi duplicati e aggiornamento della matrice allNodes
allNodes = round(allNodes,4);
[allNodes, ~, ~] = unique(allNodes, 'rows');

% Salvataggio di cooordinate dei nodi e connessioni tra di essi
numElements = length(elementNodeCoord); % Numero elementi
maxNodesPerElement = max(cellfun(@(x) size(x, 1), elementNodeCoord)); % Massimo numero di nodi 

% Inizializza la matrice degli indici con zeri (il numero di colonne dipende dal numero massimo di nodi per elemento)
elementConnections = zeros(numElements, maxNodesPerElement + 1); 

% Popola la matrice degli indici
for i = 1:numElements
    elementVertices = elementNodeCoord{i}; % Estrazione verticici dei singoli elementi
    elementVertices = round(elementVertices,4);
    numNodes = size(elementVertices, 1); % Numerod dei nodi dell'elemento
    
    % Ricerca dell'indice di ciascun vertice dell'elemento in allNodes
    for j = 1:maxNodesPerElement
        if j<=numNodes
            [~, nodeIndex] = ismember(elementVertices(j, :), allNodes, 'rows');
            elementConnections(i, j + 1) = nodeIndex;  % Memorizza l'indice del nodo
        else
            elementConnections(i, j + 1) = nodeIndex(end);
        end
    end
    
    % Memorizza il numero di nodi dell'elemento nella prima colonna
    elementConnections(i, 1) = numNodes;
end

geom.conn = elementConnections;
geom.coord = allNodes;

end