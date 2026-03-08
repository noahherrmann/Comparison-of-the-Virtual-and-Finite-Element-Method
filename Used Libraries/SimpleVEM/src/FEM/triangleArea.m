function A = triangleArea(conn, coord)

    % Calcolo della lunghezza del primo lato 1-2
    x1 = coord(conn(1), 1);
    x2 = coord(conn(2), 1);
    y1 = coord(conn(1), 2);
    y2 = coord(conn(2), 2);
    a = sqrt((x1 - x2)^2 + (y1 - y2)^2);

    % Calcolo della lunghezza del secondo lato 2-3
    x1 = coord(conn(2), 1);
    x2 = coord(conn(3), 1);
    y1 = coord(conn(2), 2);
    y2 = coord(conn(3), 2);
    b = sqrt((x1 - x2)^2 + (y1 - y2)^2);

    % Calcolo della lunghezza del terzo lato 3-1
    x1 = coord(conn(3), 1);
    x2 = coord(conn(1), 1);
    y1 = coord(conn(3), 2);
    y2 = coord(conn(1), 2);
    c = sqrt((x1 - x2)^2 + (y1 - y2)^2);

    % Calcolo del semiperimetro
    p = (a + b + c) / 2;

    % Calcolo dell'area del triangolo con la formula di Erone
    A = sqrt(p * (p - a) * (p - b) * (p - c));

end
