function qA = quadrangleArea(conn, coord)

    qA = 0.0;

    for i = 1:2

        if i == 1
            one = 1;
            two = 2;
            three = 3;
        else
            one = 1;
            two = 3;
            three = 4;
        end

        % Calcolo della lunghezza del primo lato 1-2
        x1 = coord(conn(one), 1);
        x2 = coord(conn(two), 1);
        y1 = coord(conn(one), 2);
        y2 = coord(conn(two), 2);
        a = sqrt((x1 - x2)^2 + (y1 - y2)^2);

        % Calcolo della lunghezza del secondo lato 2-3
        x1 = coord(conn(two), 1);
        x2 = coord(conn(three), 1);
        y1 = coord(conn(two), 2);
        y2 = coord(conn(three), 2);
        b = sqrt((x1 - x2)^2 + (y1 - y2)^2);

        % Calcolo della lunghezza del terzo lato 3-1
        x1 = coord(conn(three), 1);
        x2 = coord(conn(one), 1);
        y1 = coord(conn(three), 2);
        y2 = coord(conn(one), 2);
        c = sqrt((x1 - x2)^2 + (y1 - y2)^2);

        % Calcolo del semiperimetro
        p = (a + b + c) / 2;

        % Calcolo dell'area del triangolo con la formula di Erone
        tA = sqrt(p * (p - a) * (p - b) * (p - c));

        % Somma dell'area del triangolo alla quadrangolare
        qA = qA + tA;

    end

end
