function drawArrow(in, coord, value, dir)

x_max = max(coord(:,1)) - min(coord(:,1));
y_max = max(coord(:,2)) - min(coord(:,2));
% z_max = max(coord(:,3)) - min(coord(:,3));
z_max = 0;

dimensione_massima = max([x_max, y_max, z_max]);
fs = 0.05 * dimensione_massima;

for i = 1:length(in)
    if dir == 'x'
        if in(i) == 1
            x1 = coord(i,1) - sign(value)*fs;                      % Multiply for the proportion of max lentgh?
            x2 = coord(i,1);
            y1 = coord(i,2);
            y2 = coord(i,2);
            dx = x2-x1;
            dy = y2-y1;
            hold on;
            quiver(x1, y1, dx, dy, 'MaxHeadSize', 0.5, 'LineWidth', 2, 'Color', 'r');
        end
    elseif dir == 'y'
        if in(i) == 1
            x1 = coord(i,1);
            x2 = coord(i,1);
            y1 = coord(i,2) - sign(value)*fs;                      % Multiply for the proportion of max lentgh?;
            y2 = coord(i,2);
            dx = x2-x1;
            dy = y2-y1;
            hold on;
            quiver(x1, y1, dx, dy, 'MaxHeadSize', 0.5, 'LineWidth', 2, 'Color', 'r');
        end
    else
        warning('NO x OR y DIRECTION FOR ARRORW DRAW!')
    end
end
end