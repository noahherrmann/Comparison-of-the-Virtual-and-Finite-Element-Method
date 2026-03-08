function [points, triangles] = triangulateNonConvexPolygon(coord, connEl)

DEBUG = 0;

vertices = coord(connEl,:);

% Calcola l'area target del poligono
targetArea = polyarea(vertices(:,1), vertices(:,2)); % Area del poligono originale

% Inizializza i parametri per il tuning di alpha
alpha = 0; % Valore di alpha iniziale
stepSize = 1e-3; % Variazione incrementale per alpha
tolerance = 1e-5; % Tolleranza di area
maxiter = 1e5;

% Loop per ottimizzare alpha
for i = 1:maxiter % Numero massimo di iterazioni per il tuning
    shp = alphaShape(vertices(:,1), vertices(:,2), alpha); % Genera l'alphaShape
    computedArea = area(shp); % Calcola l'area dell'alphaShape
    
    % Controlla la differenza dell'area rispetto al target
    areaDiff = abs(computedArea - targetArea);
    
    if areaDiff < tolerance
        if DEBUG
            disp(['Alpha ottimizzato: ', num2str(alpha)]);
        end

        break;
    end
    
    % Incrementa alpha per la prossima iterazione
    alpha = alpha + stepSize;
end

if i == maxiter
    warning(['iterazioni Massime Raggiunte per poligoni convessi: ', num2str(i)]);
end

if DEBUG
    % Visualizzazione del risultato
    plot(shp, 'FaceColor', 'cyan'); hold on;
    title(['Triangolazione con alpha ottimizzato: ', num2str(alpha)]);
    hold off;
end

% Ottieni la connettività dei triangoli
rawTriangles  = alphaTriangulation(shp);
points = shp.Points;

% Mappa `points` su `coord` per riscrivere `triangles`
triangles = zeros(size(rawTriangles));
for j = 1:size(rawTriangles, 1)
    for k = 1:3
        % Trova l'indice di `points` che corrisponde a `coord`
        [~, idx] = min(vecnorm(coord - points(rawTriangles(j, k), :), 2, 2));
        triangles(j, k) = idx;
    end
end

end