function [skewness_el] = check_skewness(coord, conn_poly)
num_elements = size(conn_poly, 1);
skewness_el = [];

for elem = 1:num_elements
    num_sides = conn_poly(elem, 1); % Numero di lati del poligono
    element_nodes = conn_poly(elem, 2:num_sides+1); % Nodi del poligono
    element_coords = coord(element_nodes, :); % Coordinate del poligono

    % Calcolo della skewness per l'elemento corrente
    skewness_val = calculate_skewness(element_coords);
    skewness_el = [skewness_el, skewness_val];
end
end