function skewness_val = calculate_skewness(element_coords)
    % Calcolo degli angoli interni del poligono
    num_vertices = size(element_coords, 1);
    angles = zeros(num_vertices, 1);

    for i = 1:num_vertices
        p1 = element_coords(i, :);
        p2 = element_coords(mod(i, num_vertices) + 1, :);
        p3 = element_coords(mod(i+1, num_vertices) + 1, :);
        
        v1 = p2 - p1;
        v2 = p3 - p2;
        
        cos_theta = dot(v1, v2) / (norm(v1) * norm(v2));
        angles(i) = acos(cos_theta); % Angolo in radianti
    end

    % Angolo ideale per un poligono regolare (dipende dal numero di lati)
    ideal_angle = (num_vertices - 2) * pi / num_vertices;

    % Calcolo della skewness come deviazione media degli angoli
    skewness_val = mean(abs(angles - ideal_angle));

    % Calcolo della skewness come deviazione massima degli angoli
    % skewness_val = max(abs(angles - ideal_angle));
end