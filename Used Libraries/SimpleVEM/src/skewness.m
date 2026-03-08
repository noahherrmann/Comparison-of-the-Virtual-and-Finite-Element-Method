function skewness(coord, conn_poly)

skew_treshold = 0.8;

skewness_el = check_skewness(coord, conn_poly);

skewness_el_norm = (skewness_el - min(skewness_el)) / (max(skewness_el) ...
    - min(skewness_el));
skewness_vector = skwedEl(skewness_el_norm, skew_treshold);

max_skewness = max(skewness_el);
max_skewness_norm = max(skewness_el_norm);

fprintf('Max Skewness real: %f\n',max_skewness);
fprintf('Max Skewness normalized: %f\n',max_skewness_norm);

hold on;
Plot_skew(conn_poly(skewness_vector==1,:), coord, 'red');
hold on;
Plot_skew(conn_poly(skewness_vector==0,:), coord, 'green');

end