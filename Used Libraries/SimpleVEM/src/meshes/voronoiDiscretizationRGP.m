function [geom] = voronoiDiscretizationRGP(geometry, numPoints)
% DESCRIZIONE:
% Questa funzione produce una mesh composta da poligoni con numero di lati
% arbitrario, prendendo come input il contorno della geometria da meshare
% ed il numero di punti di generazione. L'algoritmo di generazione è
% quello di voronoi.
% 
% Input:
% - geometry: contorno della geometria
% - numPoints: numero di punti di generazione

% --- Sezione di posizionamento figures (commentata per velocizzare) ---
% screenSize = get(0, 'ScreenSize');
% figWidth = screenSize(3)/2; 
% figHeight = screenSize(4)/2-100; 
% xPosUP = screenSize(3) - figWidth; 
% yPosUP = screenSize(4) - figHeight-80; 
% xPosDOWN = screenSize(3) - figWidth; 
% yPosDOWN = 40; 

% Creazione di un oggetto polyshape corrispondente alla geometria in esame
geomShape = polyshape(geometry(:,1), geometry(:,2));

% Generazione di punti all'interno e attorno alla geometria
padding = 1; 
minX = min(geometry(:,1)) - padding;
maxX = max(geometry(:,1)) + padding;
minY = min(geometry(:,2)) - padding;
maxY = max(geometry(:,2)) + padding;

xContour = [minX; maxX; maxX; minX];
yContour = [minY; minY; maxY; maxY];
contour = polyshape(xContour, yContour);

pointsX = (maxX - minX) * rand(1, numPoints) + minX;
pointsY = (maxY - minY) * rand(1, numPoints) + minY;
pointsX = pointsX';
pointsY = pointsY';

while true
    % --- Plot iniziale (commentato per velocizzare) ---
    % plotPoints = figure;
    % set(plotPoints, 'Position', [xPosUP, yPosUP, figWidth, figHeight]);
    % plot(geomShape, 'FaceColor', 'none', 'EdgeColor', 'k','LineWidth',2);
    % hold on;
    % plot(pointsX, pointsY, 'r.', 'MarkerSize', 10);
    % plot(contour, 'FaceColor', 'none', 'EdgeColor', 'k','LineWidth',0.5);
    % axis equal;
    
    % Calcolo delle celle di Voronoi
    [V, C] = voronoin([pointsX pointsY]);
    
    % Inizializzazione elementi di memorizzazione
    allNodes = [];           % Memorizza i nodi
    elementNodeCoord = {};   % Memorizza gli indici dei nodi per ciascun elemento
    
    % --- Plot delle celle di Voronoi (commentato per velocizzare) ---
    % plotMesh = figure;
    % set(plotMesh, 'Position', [xPosDOWN, yPosDOWN, figWidth, figHeight]);
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
            indNaN = any(isnan(clippedCell.Vertices),2);
            VV = clippedCell.Vertices(indNaN==0,:);
            clippedCell.Vertices = VV;
            
            % Salvataggio vertici delle celle ritagliate
            if ~isempty(clippedCell.Vertices)
                % Aggiunge le coordinate dei nodi di tutti gli elementi
                allNodes = [allNodes; VV]; %#ok<AGROW>
                
                % --- Visualizzazione della cella (commentata per velocizzare) ---
                % fill(VV(:,1), VV(:,2), rand(1,3), 'EdgeColor', 'none');
                
                % Salva le coordinate dei nodi dell'elemento corrente
                elementNodeCoord{end+1} = VV; %#ok<AGROW>
            end
        end
    end
    
    % Se non si desidera raffinare ulteriormente la mesh
    rr = 2;
    % rr = input('Vuoi infittire la mesh?\n1 - yes\n2 - no\n');
    % warning('La funzione di refine è stata DISATTIVATA')
    
    if rr == 1
        % --- Sezione di raffinatezza (eventuale, commentata in parte per velocizzare) ---
        % figure(plotPoints)
        newPoints = meshRefine_VEM();
        pointsX = [pointsX; newPoints(:,1)]; %#ok<AGROW>
        pointsY = [pointsY; newPoints(:,2)]; %#ok<AGROW>
        % close(plotMesh)
        % close(plotPoints)
    else
        % close(plotMesh)
        % close(plotPoints)
        break;
    end
end

% Rimozione dei nodi duplicati e aggiornamento della matrice allNodes
allNodes = round(allNodes,4);
allNodes = unique(allNodes, 'rows');

% Salvataggio delle coordinate dei nodi e delle connessioni tra di essi
numElements = length(elementNodeCoord); 
maxNodesPerElement = max(cellfun(@(x) size(x, 1), elementNodeCoord));

% Inizializza la matrice degli indici con zeri (la prima colonna indica il numero di nodi per elemento)
elementConnections = zeros(numElements, maxNodesPerElement + 1); 

% Popola la matrice degli indici
for i = 1:numElements
    elementVertices = elementNodeCoord{i};
    elementVertices = round(elementVertices,4);
    elementVertices(2:end,1) = flip(elementVertices(2:end,1));
    elementVertices(2:end,2) = flip(elementVertices(2:end,2));
    numNodes = size(elementVertices, 1);
    
    % Ricerca dell'indice di ciascun vertice dell'elemento in allNodes
    for j = 1:maxNodesPerElement
        if j <= numNodes
            [~, nodeIndex] = ismember(elementVertices(j, :), allNodes, 'rows');
            elementConnections(i, j + 1) = nodeIndex;
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
