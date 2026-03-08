function geoCheckQuads(coord, conn)

nel = size(conn, 1);
for i = 1:nel
    delta = distortionQL1(conn(i, :), coord);

    one = conn(i, 1);
    two = conn(i, 2);
    three = conn(i, 3);
    four = conn(i, 4);

    x1 = coord(one, 1);
    x2 = coord(two, 1);
    x3 = coord(three, 1);
    x4 = coord(four, 1);

    y1 = coord(one, 2);
    y2 = coord(two, 2);
    y3 = coord(three, 2);
    y4 = coord(four, 2);

    xx = [x1 x2 x3 x4 x1];
    yy = [y1 y2 y3 y4 y1];

    AreaGauss = area_FEM(coord, conn, 1);

    if AreaGauss > 0
        if delta >= 0.2
            disp(['WARNING: elemento ', num2str(i), ' sembra essere distorto!']);
            fill(xx, yy, 'r');
            hold on;
            plot(xx, yy, 'k-', 'LineWidth', 1);
        else
            fill(xx, yy, 'g');
            hold on;
            plot(xx, yy, 'k-', 'LineWidth', 1);
        end
    else
        disp(['WARNING: elemento ', num2str(i), ' ha area negativa!']);
        fill(xx, yy, 'b');
        plot(xx, yy, 'k-', 'LineWidth', 1);
        hold on;
    end
end
hold off;
end
