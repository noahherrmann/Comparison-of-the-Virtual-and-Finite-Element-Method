function [curve1, curve2] = equidistantCurves(originalCurveCoord, distance)

COORD = originalCurveCoord;
d = distance;
n = size(COORD, 1);

% Inizializza le due curve equidistanti e simmetriche
curve1 = zeros(n, 2);
curve2 = zeros(n, 2);

for i = 1:n-1
    % Segmento tra i punti i e i+1
    p1 = COORD(i, :);
    p2 = COORD(i+1, :);
    
    % Vettore del segmento
    vec = p2 - p1;
    
    % Direzione normale
    normal = [-vec(2), vec(1)];
    
    % Normalizzazione direzione normale
    normal = normal / norm(normal);
    
    % Sposta i punti lungo la normale di una distanza d
    curve1(i, :) = p1 + d * normal;
    curve2(i, :) = p1 - d * normal;
end

% Ultimo punto
p1 = COORD(end-1, :);
p2 = COORD(end, :);
vec = p2 - p1;
normal = [-vec(2), vec(1)];
normal = normal / norm(normal);
curve1(end, :) = p1 + d * normal;
curve2(end, :) = p1 - d * normal; 