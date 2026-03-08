function orientation = checkOrientation(coord, conn)
    numElements = size(conn, 1);
    orientation = zeros(numElements, 1);
    for i = 1:numElements
        elementNodes = conn(i, 2:conn(i, 1) + 1); % Ottieni i nodi dell'elemento
        elementCoords = coord(elementNodes, :); % Ottieni le coordinate dei nodi
        % Calcola l'area con segno usando la formula del poligono
        signedArea = 0;
        numNodes = length(elementNodes);
        for j = 1:numNodes
            x1 = elementCoords(j, 1);
            y1 = elementCoords(j, 2);
            x2 = elementCoords(mod(j, numNodes) + 1, 1);
            y2 = elementCoords(mod(j, numNodes) + 1, 2);
            signedArea = signedArea + (x1 * y2 - x2 * y1);
        end
        signedArea = signedArea / 2;
        % Determina l'orientamento in base all'area con segno
        if signedArea > 0
            orientation(i) = 1;
        else
            orientation(i) = 0;
        end
    end
end