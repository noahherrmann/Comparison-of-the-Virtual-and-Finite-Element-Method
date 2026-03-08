function geoCheckTris(coord, conn)

VALORE_SOGLIA = 10;

nel = size(conn, 1);
for i = 1:nel
    x1 = coord(conn(i,1), 1);
    x2 = coord(conn(i,2), 1);
    x3 = coord(conn(i,3), 1);

    y1 = coord(conn(i,1), 2);
    y2 = coord(conn(i,2), 2);
    y3 = coord(conn(i,3), 2);

    b1 = y3 - y2;
    b2 = y1 - y3;
    b3 = y1 - y2;

    c1 = x2 - x3;
    c2 = x1 - x3;
    c3 = x1 - x2;

    % Calcolo del Jacobiano
    J = [-c3, -b3;
        -c2, -b2];

    % Calcolo dell'area
    A = det(J) / 2;

    % Matrice Bs
    Bs = [-1, 1, 0;
        -1, 0, 1];

    % Matrice B nel sistema cartesiano
    Bi = J \ Bs;

    % Matrice di rigidezza
    K = A * (Bi') * Bi;

    % Calcolo degli autovalori
    eigv = eig(K);

    % Rimuovere gli autovalori pari a 0
    eps = 1e-10;
    eigv2 = eigv(eigv > eps);

    % Fattore di condizionamento normalizzato
    cond = (max(eigv2) / min(eigv2)) / 3;

    xx = [x1 x2 x3 x1];
    yy = [y1 y2 y3 y1];

    if A > 0
        if cond >= VALORE_SOGLIA
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
