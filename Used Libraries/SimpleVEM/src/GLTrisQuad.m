%Gauss-Legendre Triangle Quadrature
function [transformedPoints, transformedWeights] = GLTrisQuad(vertices, nPoints)
    % Estrai i vertici
    x1 = vertices(1,1);
    y1 = vertices(1,2);
    x2 = vertices(2,1);
    y2 = vertices(2,2);
    x3 = vertices(3,1);
    y3 = vertices(3,2);

    % Definizione dei punti e dei pesi in base al numero di punti desiderato
    switch nPoints
        case 1
            weights = 1/2;
            points = [1/3, 1/3];
        case 3
            weights = [1/6, 1/6, 1/6];
            points = [1/2, 0; 1/2, 1/2; 0, 1/2];
        otherwise
            error('Numero di punti non supportato');
    end

    % Calcola il determinante del Jacobiano della trasformazione affine dal triangolo di riferimento al triangolo specifico
    jacobianDet = abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));

    % Trasforma i punti nel nuovo triangolo
    numPoints = size(points, 1);
    transformedPoints = zeros(numPoints, 2);
    for i = 1:numPoints
        transformedPoints(i,1) = x1 + points(i,1) * (x2-x1) + points(i,2) * (x3-x1);
        transformedPoints(i,2) = y1 + points(i,1) * (y2-y1) + points(i,2) * (y3-y1);
    end

    % Adegua i pesi moltiplicando per il determinante del Jacobiano
    transformedWeights = weights * jacobianDet;
end
