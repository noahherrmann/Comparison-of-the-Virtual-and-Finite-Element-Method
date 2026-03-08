function [conn_tris, conn_quads, conn_poly] = separate_conn(conn)

    % Inizializzazione delle matrici vuote per i triangoli e i quadrilateri
    conn_tris = [];
    conn_quads = [];
    conn_poly = [];

    nel = size(conn, 1);
    for i = 1:nel
        if conn(i, 1) == 3
            conn_tris = [conn_tris; conn(i, 2:4)];
        elseif conn(i, 1) == 4
            conn_quads = [conn_quads; conn(i, 2:5)];
        else
            conn_poly = [conn_poly; conn(i, 1:end)];
        end
    end

end
