function normals = calculateNormals(vertices, edges)
    % Numero di lati
    n = size(edges, 1);
    
    % Inizializzazione dell'array delle normali
    normals = zeros(n, 2);  % Cambia la dimensione se i tuoi punti sono 3D
    
    % Calcolo dei vettori e delle normali
    for i = 1:n
        % Indici dei vertici che formano il lato
        v1 = edges(i, 1);
        v2 = edges(i, 2);
        
        % Vettore dal primo al secondo vertice
        vector = vertices(v2, :) - vertices(v1, :);
        
        % Calcolo della normale (perpendicolare al vettore)
        normals(i, :) = [vector(2), -vector(1)];  % Ruota di 90 gradi
    end
    
    % Normalizzazione delle normali (opzionale)
    for i = 1:n
        normals(i, :) = normals(i, :) / norm(normals(i, :));
        if isnan(normals(i, :))
            normals(i, :) = [0 , 0];
        end
    end
end
