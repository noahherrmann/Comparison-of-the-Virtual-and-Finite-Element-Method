function [diric, neum_c, neum_d] = imposeBC(coord, conn, options)

% Enlarge limits:
% Get current x and y limits
x_limits = xlim;
y_limits = ylim;

p = options.p;

percentage = 10;
% Calculate enlargement
x_enlargement = diff(x_limits) * percentage / 100;
y_enlargement = diff(y_limits) * percentage / 100;

% Set new limits
xlim([x_limits(1) - x_enlargement, x_limits(2) + x_enlargement]);
ylim([y_limits(1) - y_enlargement, y_limits(2) + y_enlargement]);

fprintf('DEFINE CONSTRAINTS:\n')
listnode = linspace(1,size(coord,1),size(coord,1))';

%% IMPOSING DIRICHLET
fprintf('DIRICHLET:\n')
tmp = 1;
diric = zeros(length(listnode),3);
diric(:,1) = listnode;
while tmp == 1
    % Get the polygon vertices from the user
    h = impoly;
    wait(h); % Wait for the user to double-click
    polygon = getPosition(h);
    in = findInside(polygon(:,1), polygon(:,2), coord);
    hold on;
    % plot(coord(in==1,1),coord(in==1,2),'sr');
    plotNodes([coord(in==1,1),coord(in==1,2)],'red')
    cc = input('1) Constraint in x\n2) Constraint in y\n3) Constraint in xy\n4) Imposed displacement\n\n');

    if cc == 1
        diric(in==1,2) = 10;
        diric(in==1,3) = 0;
        diric(in==1,4) = 0;
    elseif cc == 2
        diric(in==1,2) = 01;
        diric(in==1,3) = 0;
        diric(in==1,4) = 0;
    elseif cc == 3
        diric(in==1,2) = 11;
        diric(in==1,3) = 0;
        diric(in==1,4) = 0;
    elseif cc == 4
        val_x = input('Displacement over x: ');
        val_y = input('\nDisplacement over y: ');
        diric(in==1,2) = 2;
        diric(in==1,3) = val_x;
        diric(in==1,4) = val_y;
    else
        warning('No options!')
    end
    
    tmp = input('Assign new constraint?\n1 - yes\n2 - no\n');
    clc;
    % splashScreen();
end

fprintf('DIRICHLET ASSIGNED!\n\n')

%% IMPOSING NEUMANN
fprintf('NEUMANN:\n')
tmp = 1;
neum_c = [];
conc = input('Puntual Load\n1 - yes\n2 - no\n');
if conc == 1
    neum_c = zeros(length(listnode),3);
    neum_c(:,1) = listnode;
    while tmp == 1
        % Carichi concentrati
        % Get the polygon vertices from the user
        h = impoly;
        wait(h); % Wait for the user to double-click
        polygon = getPosition(h);
        in = findInside(polygon(:,1), polygon(:,2), coord);
        iw = input('is it x or y direction?\n1 - x\n2 - y\n3 - xy\n');
        if iw == 3
            value_x = input('Value x?: ');
            value_y = input('Value y?: ');
        else
            value = input('Value?: ');
        end

        % Fx
        if iw == 1
            drawArrow(in, coord, value, 'x')
            neum_c(in==1,2) = value;
        % Fy
        elseif iw == 2
            drawArrow(in, coord, value, 'y')
            neum_c(in==1,3) = value;
        elseif iw == 3
            drawArrow(in, coord, value_x, 'x')
            drawArrow(in, coord, value_y, 'y')
            neum_c(in==1,2) = value_x;
            neum_c(in==1,3) = value_y;
        else
            warning('No options!')
        end
        hold on;
        plot(coord(in==1,1),coord(in==1,2),'ob');
        % estrai nodi e applica (chiedi valore)
        tmp = input('Add another Puntual Load?\n1 - yes\n2 - no\n');
        clc;
        % splashScreen();
    end
end

clc;
% splashScreen();

fprintf('DIRICHLET ASSIGNED!\n\n')
fprintf('NEUMANN:\n')
% Carichi Distribuiti
dist = input('Distributed Load\n1 - yes\n2 - no\n');
tmp = 1;
neum_d = [];
if dist == 1
    while tmp == 1
        title('Click near edges to select them. Press Enter when done.');

        % Get the polygon vertices from the user
        h = impoly;
        wait(h); % Wait for the user to double-click
        polygon = getPosition(h);

        selected_edges_info = []; % Inizializza un array per tenere traccia di elementi e lati selezionati

        for i = 1:size(conn, 1)
            % Definizione dei lati in base ai nodi
            if p == 1
                nedgesEl = conn(i,1);
                connEl = conn(i,2:nedgesEl+1);
                sides = edgeExtract(nedgesEl, connEl);
            elseif p == 2

            else
                error('no p > 2');
            end

            for j = 1:size(sides, 1)

                if p == 1
                    pt1 = coord(sides(j, 1), :);
                    pt2 = coord(sides(j, 2), :);
                    ptMid = (pt2+pt1)/2;
                elseif p == 2

                else
                    error('no p > 2');
                end
                
                % Calcolo del punto rappresentativo usando la curva di Bézier
                t = 0.5;
                representativePoint = (1 - t)^2 * pt1 + 2*(1 - t)*t * ptMid + t^2 * pt2;

                % Verifica se il punto rappresentativo è all'interno del poligono
                if inpolygon(representativePoint(1), representativePoint(2), polygon(:, 1), polygon(:, 2))
                    plot([pt1(1), ptMid(1), pt2(1)], [pt1(2), ptMid(2), pt2(2)], 'r-', 'LineWidth', 2);
                    selected_edges_info = [selected_edges_info; i, j]; % Salva l'elemento e l'indice del lato
                end
            end
        end

        value_x = input('Normal pressure (+) traction = ');
        value_y = input('Tangential pressure = ');
        for i = 1:size(selected_edges_info, 1)
            element_number = selected_edges_info(i, 1);
            side_number = selected_edges_info(i, 2);
            neum_d = [neum_d; element_number, side_number, value_x, value_y];
        end

        tmp = input('Add another Distributed Load?\n1 - yes\n2 - no\n');
        clc;
        % splashScreen();
    end
end

clc;
% splashScreen();

end