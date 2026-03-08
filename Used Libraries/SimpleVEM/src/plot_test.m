% Definire i vertici del poligono (x, y) per ogni elemento
vertices = [0 0; 1 0; 1 1; 0.5 2; 0 1]; % Esempio di quadrilatero
solution_values = [0.1, 0.5, 0.3, 0.8, 0.]; % Valori della soluzione nodale

% Creare un grafico colorato con 'patch'
figure;
patch('Vertices', vertices, 'Faces', [1 2 3 4 5], 'FaceVertexCData', solution_values', ...
      'FaceColor', 'interp', 'EdgeColor', 'k');
colorbar; % Aggiungi una barra dei colori per visualizzare i valori
axis equal;
title('Soluzione nodale interpolata su un poligono');
