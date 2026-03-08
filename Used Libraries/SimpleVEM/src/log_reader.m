function T = log_reader()
% RUN ME INSIDE THE FOLDER WHERE YOU HAVE THE simulation_log.txt file!!

fid = fopen('simulation_log.txt','r');
if fid == -1
    error('Impossibile aprire simulation_log.txt');
end

% Inizializzazione degli array per i dati
cases = [];
methods = {};
elements = [];
dofs = [];
max_disp = [];
energy = [];

% Variabili per il blocco corrente
current_case = NaN;
current_method = '';
current_elements = NaN;
current_dofs = NaN;
current_energy = NaN;
disp_x = NaN;
disp_y = NaN;

% Funzione interna per salvare il blocco corrente (se valido)
    function storeBlock()
        if ~isempty(current_method)
            if ~isnan(disp_x) && ~isnan(disp_y)
                current_max_disp = sqrt(disp_x^2 + disp_y^2);
            else
                current_max_disp = NaN;
            end
            cases(end+1,1) = current_case;
            methods{end+1,1} = current_method;
            elements(end+1,1) = current_elements;
            dofs(end+1,1) = current_dofs;
            max_disp(end+1,1) = current_max_disp;
            energy(end+1,1) = current_energy;
        end
    end

while ~feof(fid)
    tline = fgetl(fid);
    if ischar(tline)
        % Rileva la riga "Processing Case X ..." per aggiornare il caso corrente
        tokens = regexp(tline, 'Processing Case\s+(\d+)', 'tokens');
        if ~isempty(tokens)
            current_case = str2double(tokens{1}{1});
            continue;
        end
        
        % Rileva la riga contenente il metodo (es. "-> FEM-quad-1")
        tokens = regexp(tline, '->\s*([A-Z]+-[a-z]+)-\d+', 'tokens');
        if ~isempty(tokens)
            % Salva il blocco precedente se esiste
            storeBlock();
            % Imposta il nuovo metodo e resetta le variabili
            current_method = tokens{1}{1};
            current_elements = NaN;
            current_dofs = NaN;
            current_energy = NaN;
            disp_x = NaN;
            disp_y = NaN;
            continue;
        end
        
        % Estrae il numero totale di elementi
        tokens = regexp(tline, 'Total number of elements:\s*([\d\.]+)', 'tokens');
        if ~isempty(tokens)
            current_elements = str2double(tokens{1}{1});
            continue;
        end
        
        % Estrae il numero di DOFs (Reduced system size)
        tokens = regexp(tline, 'Reduced system size:\s*([\d\.]+)', 'tokens');
        if ~isempty(tokens)
            current_dofs = str2double(tokens{1}{1});
            continue;
        end
        
        % Estrae la massima deformazione in x
        tokens = regexp(tline, 'Maximum displacement in x:\s*([\d\.Ee\-\+]+)', 'tokens');
        if ~isempty(tokens)
            disp_x = str2double(tokens{1}{1});
            continue;
        end
        
        % Estrae la massima deformazione in y
        tokens = regexp(tline, 'Maximum displacement in y:\s*([\d\.Ee\-\+]+)', 'tokens');
        if ~isempty(tokens)
            disp_y = str2double(tokens{1}{1});
            continue;
        end
        
        % Estrae l'energia elastica totale
        tokens = regexp(tline, 'Total elastic energy:\s*([\d\.\-Ee\+]+)', 'tokens');
        if ~isempty(tokens)
            current_energy = str2double(tokens{1}{1});
            continue;
        end
    end
end

% Salva l'ultimo blocco
storeBlock();
fclose(fid);

% Crea la tabella e limita a 9 righe se ce ne sono di più
T = table(cases, methods, elements, dofs, max_disp, energy, ...
    'VariableNames', {'Case','Method','Elements','DOFs','MaxDisplacement','Energy'});
if height(T) > 9
    T = T(1:9,:);
end

% Esporta la tabella in formato LaTeX
exportLatexTable(T, 'simulation_table.tex');

end

function exportLatexTable(T, filename)
% EXPORTLATEXTABLE esporta la tabella MATLAB T in formato LaTeX nel file specificato.
fid = fopen(filename, 'w');
if fid == -1
    error('Impossibile aprire il file %s per la scrittura', filename);
end
fprintf(fid, '\\begin{table}[ht]\\centering\n');
fprintf(fid, '\\begin{tabular}{cccccc}\n\\noalign{\\hrule height 1.5pt}\n');
fprintf(fid, '\\textbf{Case} & \\textbf{Method} & \\textbf{Elements} & \\textbf{DOFs} & \\textbf{Max Displacement} & \\textbf{Energy} \\\\ \n \\noalign{\\hrule height 1.5pt}\n');
for i = 1:height(T)
    fprintf(fid, '%d & %s & %d & %d & %.4f & %.4f \\\\ \\hline\n', ...
        T.Case(i), T.Method{i}, T.Elements(i), T.DOFs(i), T.MaxDisplacement(i), T.Energy(i));
end
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\caption{Simulation results.}\n');
fprintf(fid, '\\label{tab:simulation_results}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);
end
