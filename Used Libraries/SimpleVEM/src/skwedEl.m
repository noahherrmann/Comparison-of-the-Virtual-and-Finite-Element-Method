function skewness_vector = skwedEl(skewness_el, threshold)
num_elements = length(skewness_el);
skewness_vector = zeros(num_elements, 1);

for elem = 1:num_elements
    % Confronto con la soglia
    skewness_val = skewness_el(elem);
    if skewness_val > threshold
        skewness_vector(elem) = 1; % Skewness elevata
    else
        skewness_vector(elem) = 0; % Skewness ok
    end
end
end