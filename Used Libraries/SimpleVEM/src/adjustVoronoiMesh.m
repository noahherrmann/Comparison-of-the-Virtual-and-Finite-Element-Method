function geom_VEM_voro = adjustVoronoiMesh(geometry, nelem_fem, n, max_iter)
% adjustVoronoiMesh: Esegue la discretizzazione Voronoi fino a quando il
% numero di elementi ottenuti è entro un range di ±5%% rispetto a nelem_fem,
% utilizzando una ricerca binaria per ottimizzare il parametro n.
%
% INPUT:
%   geometry  - Contorno della geometria da meshare.
%   nelem_fem - Numero di elementi target.
%   n         - Valore iniziale del parametro per la discretizzazione.
%   max_iter  - Numero massimo di iterazioni per evitare loop infiniti.
%
% OUTPUT:
%   geom_VEM_voro - Struttura della mesh Voronoi ottenuta, contenente ad esempio il campo .conn.
%
% NOTA:
%   Si assume che la relazione tra n e il numero di elementi sia (quasi) monotona.

    tolerance = 0.05; % ±5%

    % Prima valutazione con il valore iniziale
    geom_VEM_voro = voronoiDiscretizationRGP(geometry, n);
    nelem_voro = size(geom_VEM_voro.conn, 1);
    
    % Se il valore iniziale rientra già nella tolleranza, restituisce subito
    if nelem_voro >= (1-tolerance)*nelem_fem && nelem_voro <= (1+tolerance)*nelem_fem
        fprintf('Match iniziale: n = %d, nelem_voro = %d\n', n, nelem_voro);
        return;
    end
    
    % --- Fase di bracketing ---
    if nelem_voro < nelem_fem
        % Se la mesh ha meno elementi del target, aumentiamo n
        n_low = n;
        n_high = n;
        for i = 1:max_iter
            n_high = n_high * 2;  % raddoppia n_high
            geom_VEM_voro = voronoiDiscretizationRGP(geometry, n_high);
            nelem_voro = size(geom_VEM_voro.conn, 1);
            if nelem_voro >= nelem_fem
                break;
            end
        end
    else
        % Se la mesh ha più elementi del target, riduciamo n
        n_high = n;
        n_low = n;
        for i = 1:max_iter
            n_low = max(1, floor(n_low / 2));  % dimezza n_low (non scende sotto 1)
            geom_VEM_voro = voronoiDiscretizationRGP(geometry, n_low);
            nelem_voro = size(geom_VEM_voro.conn, 1);
            if nelem_voro <= nelem_fem
                break;
            end
        end
    end

    % --- Fase di ricerca binaria ---
    for iter = 1:max_iter
        n_mid = round((n_low + n_high) / 2);
        geom_VEM_voro = voronoiDiscretizationRGP(geometry, n_mid);
        nelem_voro = size(geom_VEM_voro.conn, 1);
        fprintf('Iterazione %d: n_mid = %d, nelem_voro = %d\n', iter, n_mid, nelem_voro);
        
        % Verifica se il numero di elementi è entro la tolleranza
        if nelem_voro >= (1-tolerance)*nelem_fem && nelem_voro <= (1+tolerance)*nelem_fem
            fprintf('Match trovato: n_mid = %d, nelem_voro = %d\n', n_mid, nelem_voro);
            return;
        elseif nelem_voro < nelem_fem
            n_low = n_mid + 1;
        else
            n_high = n_mid - 1;
        end
        
        % Se l'intervallo si annulla, esce dal ciclo
        if n_low > n_high
            break;
        end
    end

    warning('Il numero di elementi desiderato non è stato raggiunto entro il numero massimo di iterazioni.');
end
