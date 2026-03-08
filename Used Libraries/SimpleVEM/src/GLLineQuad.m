%Gauss-Lobatto Line Segment
function transformedWeights = GLLineQuad(point1, point2)
    % Nodi e pesi per Gauss-Lobatto con 2 punti
    nodes = [-1, 1];
    weights = [1, 1];

    % Calcola la lunghezza del segmento
    length = norm([point2(1) - point1(1), point2(2) - point1(2)]);

    % Trasforma i pesi
    transformedWeights = weights * (length / 2);

    % Ritorna i pesi trasformati
end