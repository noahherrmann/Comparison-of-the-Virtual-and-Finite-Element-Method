% DESCRIZIONE:
% Questa funzione Permette un facile remeshing dando la possibilità di inserire, in tre modalità, 
% un numero di punti arbitrario all'interno della matrice newPoints, tramite la selezione 
% di puntio aree di interesse direttamente su una figure precedentemente creata.


function [newPoints] = meshRefine_VEM()

refineType = input('Che tipo di refine vuoi usare?\n1 - Puntual Refine\n2 - Rectangular area Refine\n3 - Circular area Refine\n');
% clc

newPoints = [];


if refineType == 1

        while true
            % Seleziona un punto con il mouse
            % figure(plotPoints)
            [x_point, y_point, button] = ginput(1);
            
            % Se l'utente preme Invio (button vuoto), esci dal ciclo
            if isempty(button)
                break;
            end
            
            % Salvataggio del punto selezionato
            newPoints = [newPoints; x_point, y_point]; %#ok<AGROW>
    
            % Plotta tutti i punti selezionati fino ad ora
            plot(newPoints(:,1),newPoints(:,2), 'b.', 'MarkerSize', 15, 'LineWidth', 2);    
        end

elseif refineType == 2 

    % selezione area rettangolare di interesse
    disp('Clicca e trascina per disegnare un rettangolo di selezione');
    rect = getrect; 
    % clc

    % Estrai i limiti del rettangolo selezionato
    x_min = rect(1);
    y_min = rect(2);
    width = rect(3);
    height = rect(4);
    x_max = x_min + width;
    y_max = y_min + height;

    rectX = [x_min; x_max; x_max; x_min];
    rectY = [y_min; y_min; y_max; y_max];
    rectangle = polyshape(rectX', rectY');
    plot(rectangle, 'EdgeColor', 'b','LineWidth',0.5); 
    
    % Numero di punti da generare
    numPoints = input('Quanti punti vuoi che siano aggiunti?\n');
    % clc
  
    % Genera punti casuali all'interno del rettangolo selezionato
    x_point = x_min + (x_max - x_min) * rand(numPoints, 1);
    y_point = y_min + (y_max - y_min) * rand(numPoints, 1);
    
    % Salvataggio dei punti selezionati
    newPoints = [newPoints; x_point, y_point]; 
    % Salvataggio dei punti selezionati
    plot(newPoints(:,1),newPoints(:,2), 'b.', 'MarkerSize', 15, 'LineWidth', 2); 
    pause(2)
else

    % Seleziona il centro del cerchio interattivamente
    disp('Seleziona il centro del cerchio\n');
    center = ginput(1); % Ottieni le coordinate (x,y) del centro
    % clc

    % Seleziona un punto sulla circonferenza per definire il raggio
    disp('Seleziona un punto sulla circonferenza del cerchio per definire il valore del raggio\n');
    circPoint = ginput(1); % Ottieni il punto sulla circonferenza
    % clc

    % Calcola il raggio del cerchio
    radius = sqrt((circPoint(1) - center(1))^2 + (circPoint(2) - center(2))^2);
    
    % Disegna il cerchio selezionato
    theta = linspace(0, 2*pi, 100); % Parametro angolare per creare il cerchio
    x_circle = center(1) + radius * cos(theta);
    y_circle = center(2) + radius * sin(theta);
    plot(x_circle, y_circle, 'r-', 'LineWidth', 1.5); % Plot del cerchio
    
    % Numero di punti da generare
    numPoints = input('Quanti punti vuoi che siano aggiunti?\n');
    % clc

    % Genera raggi casuali scalati (distribuzione uniforme nel cerchio)
    r = radius * sqrt(rand(numPoints, 1));
    
    % Genera angoli casuali tra 0 e 2pi
    theta = 2 * pi * rand(numPoints, 1);
    
    % Coordinate dei punti casuali all'interno del cerchio
    x_point = center(1) + r .* cos(theta);
    y_point = center(2) + r .* sin(theta);

    % Salvataggio dei punti selezionati
    newPoints = [newPoints; x_point, y_point]; 
    % Plot dei punti generati
    plot(x_point, y_point, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k');
    pause(2)
end

end
