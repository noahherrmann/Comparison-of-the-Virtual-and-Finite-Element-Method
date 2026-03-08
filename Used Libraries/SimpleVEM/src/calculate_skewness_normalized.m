function skewness_val_norm = calculate_skewness_normalized(element_coords)
    % Calcolo degli angoli interni del poligono
    num_vertices = size(element_coords, 1);
    angles = zeros(num_vertices, 1);

    for i = 1:num_vertices
        % Coordinate dei vertici coinvolti
        p1 = element_coords(i, :);
        p2 = element_coords(mod(i, num_vertices) + 1, :);
        p3 = element_coords(mod(i+1, num_vertices) + 1, :);
        
        % Vettori tra i vertici
        v1 = p2 - p1;
        v2 = p3 - p2;
        
        % Calcolo del coseno dell'angolo tra i due vettori
        cos_theta = dot(v1, v2) / (norm(v1) * norm(v2));
        angles(i) = acos(cos_theta); % Angolo in radianti
    end

    % Angolo ideale per un poligono regolare (dipende dal numero di lati)
    ideal_angle = (num_vertices - 2) * pi / num_vertices;

    % Calcolo della skewness come deviazione media degli angoli
    skewness_val = mean(abs(angles - ideal_angle));

    % Calcolo della skewness massima (approssimazione: un angolo prossimo a 0 o 180 gradi)
    skewness_max = pi - ideal_angle;  % Limite massimo quando uno degli angoli è prossimo a pi (180 gradi)

    % Normalizzazione della skewness tra 0 e 1
    skewness_val_norm = skewness_val / skewness_max;
end
