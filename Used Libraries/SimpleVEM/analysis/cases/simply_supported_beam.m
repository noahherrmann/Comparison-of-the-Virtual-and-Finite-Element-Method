%simply supported beam

close all
clear
clc

% -------------------------------------------------------------------------
% 1) Suppress All Warnings
% -------------------------------------------------------------------------
warning('off', 'all')

% -------------------------------------------------------------------------
% 2) Prompt for Output Folder
% -------------------------------------------------------------------------
outputFolder = uigetdir('', 'Select the Output Folder for Saving Results');
if outputFolder == 0
    error('No folder selected. Exiting the script.');
end

% -------------------------------------------------------------------------
% 3) Parameters and Tau Values
% -------------------------------------------------------------------------
p = 1;
tau_opt = 0.13;              % Example additional tau value
tau_values = [1, tau_opt];     % We will run simulations for these tau values

% Create an array of folder names (one for each tau value)
tau_folder_names = cell(size(tau_values));
for i = 1:length(tau_values)
    tau_str = num2str(tau_values(i));
    % Replace the dot with an underscore for folder naming if necessary
    tau_str = strrep(tau_str, '.', '_');
    subfolder = fullfile(outputFolder, ['tau_' tau_str]);
    if ~exist(subfolder, 'dir')
        mkdir(subfolder);
    end
    tau_folder_names{i} = subfolder;
end

% -------------------------------------------------------------------------
% 4) Create the Geometries
% -------------------------------------------------------------------------
[~, contour, options] = AllQuadrilateralMesh(15, 2, 5, 3);
geometry = contour;

nx = 25;
ny = 7;

% Quadrilateral geometries (three refinements)
[geom1, contour1] = AllQuadrilateralMesh(15, 2, nx,   ny);
[geom2, contour2] = AllQuadrilateralMesh(15, 2, 2*nx, 2*ny);
[geom3, contour3] = AllQuadrilateralMesh(15, 2, 3*nx, 3*ny);

% FEM: remove the first column of connectivity
geom_FEM_quad_1 = geom1;  geom_FEM_quad_1.conn = geom_FEM_quad_1.conn(:, 2:end);
geom_FEM_quad_2 = geom2;  geom_FEM_quad_2.conn = geom_FEM_quad_2.conn(:, 2:end);
geom_FEM_quad_3 = geom3;  geom_FEM_quad_3.conn = geom_FEM_quad_3.conn(:, 2:end);

% VEM: same quadrilateral meshes
geom_VEM_quad_1 = geom1;
geom_VEM_quad_2 = geom2;
geom_VEM_quad_3 = geom3;

% Voronoi geometries for VEM
nelem_fem_1 = size(geom_FEM_quad_1.conn,1);
n1 = size(geom_FEM_quad_1.coord,1);
geom_VEM_voro_1 = adjustVoronoiMesh(geometry, nelem_fem_1, n1, 1000);

n2 = size(geom_FEM_quad_2.coord,1);
nelem_fem_2 = size(geom_FEM_quad_2.conn,1);
geom_VEM_voro_2 = adjustVoronoiMesh(geometry, nelem_fem_2, n2, 1000);

n3 = size(geom_FEM_quad_3.coord,1);
nelem_fem_3 = size(geom_FEM_quad_3.conn,1);
geom_VEM_voro_3 = adjustVoronoiMesh(geometry, nelem_fem_3, n3, 1000);
% geom_VEM_voro_1 = voronoiDiscretizationRGP(geometry,   nx*ny);
% geom_VEM_voro_2 = voronoiDiscretizationRGP(geometry, 4*nx*ny);
% geom_VEM_voro_3 = voronoiDiscretizationRGP(geometry, 9*nx*ny);

% (Optional) check orientation of Voronoi meshes
orientation1 = checkOrientation(geom_VEM_voro_1.coord, geom_VEM_voro_1.conn);
if any(orientation1 ~= 1)
    error('Stop calculation!');
end
orientation2 = checkOrientation(geom_VEM_voro_2.coord, geom_VEM_voro_2.conn);
if any(orientation2 ~= 1)
    error('Stop calculation!');
end
orientation3 = checkOrientation(geom_VEM_voro_3.coord, geom_VEM_voro_3.conn);
if any(orientation3 ~= 1)
    error('Stop calculation!');
end

% Additional options
options.p = p;
options.problemtype = 4;  % Example problem type

% -------------------------------------------------------------------------
% 5) Material Definition
% -------------------------------------------------------------------------
mat = linearElastic2D_Steel;  % Replace with your own material definition if needed

% -------------------------------------------------------------------------
% 6) MAIN LOOP over all tau values
% -------------------------------------------------------------------------
fprintf('[SUPPORTED BEAM]\n');
for tau_idx = 1:length(tau_values)
    current_tau = tau_values(tau_idx);
    resultsFolder = tau_folder_names{tau_idx};
    
    % 6.1) Start logging (diary) in the subfolder
    diaryFile = fullfile(resultsFolder, 'simulation_log.txt');
    diary(diaryFile);
    diary on;
    
    fprintf('>>> Starting simulation for tau = %f\n\n', current_tau);
    
    % ---------------------------------------------------------------------
    % 6.2) Pre-processing Plots and Boundary Conditions (Quad and Voronoi)
    % ---------------------------------------------------------------------
    % --- Quad 1
    Plot_MeshFEM(geom_VEM_quad_1.conn, geom_VEM_quad_1.coord, 1, '-', 0);
    [BC1_quad.diric, BC1_quad.neum_c, BC1_quad.neum_d] = ...
        normalBC(geom_VEM_quad_1.coord, geom_VEM_quad_1.conn, options);
    saveas(gcf, fullfile(resultsFolder, 'Prepro_Quad_1.png'));
    
    % --- Quad 2
    Plot_MeshFEM(geom_VEM_quad_2.conn, geom_VEM_quad_2.coord, 1, '-', 0);
    [BC2_quad.diric, BC2_quad.neum_c, BC2_quad.neum_d] = ...
        normalBC(geom_VEM_quad_2.coord, geom_VEM_quad_2.conn, options);
    saveas(gcf, fullfile(resultsFolder, 'Prepro_Quad_2.png'));
    
    % --- Quad 3
    Plot_MeshFEM(geom_VEM_quad_3.conn, geom_VEM_quad_3.coord, 1, '-', 0);
    [BC3_quad.diric, BC3_quad.neum_c, BC3_quad.neum_d] = ...
        normalBC(geom_VEM_quad_3.coord, geom_VEM_quad_3.conn, options);
    saveas(gcf, fullfile(resultsFolder, 'Prepro_Quad_3.png'));
    
    % --- Voro 1
    Plot_MeshFEM(geom_VEM_voro_1.conn, geom_VEM_voro_1.coord, 1, '-', 0);
    [BC1_voro.diric, BC1_voro.neum_c, BC1_voro.neum_d] = ...
        normalBC(geom_VEM_voro_1.coord, geom_VEM_voro_1.conn, options);
    saveas(gcf, fullfile(resultsFolder, 'Prepro_Voro_1.png'));
    
    % --- Voro 2
    Plot_MeshFEM(geom_VEM_voro_2.conn, geom_VEM_voro_2.coord, 1, '-', 0);
    [BC2_voro.diric, BC2_voro.neum_c, BC2_voro.neum_d] = ...
        normalBC(geom_VEM_voro_2.coord, geom_VEM_voro_2.conn, options);
    saveas(gcf, fullfile(resultsFolder, 'Prepro_Voro_2.png'));
    
    % --- Voro 3
    Plot_MeshFEM(geom_VEM_voro_3.conn, geom_VEM_voro_3.coord, 1, '-', 0);
    [BC3_voro.diric, BC3_voro.neum_c, BC3_voro.neum_d] = ...
        normalBC(geom_VEM_voro_3.coord, geom_VEM_voro_3.conn, options);
    saveas(gcf, fullfile(resultsFolder, 'Prepro_Voro_3.png'));
    
    % ---------------------------------------------------------------------
    % 6.3) Initialize Arrays for Energy Values
    % ---------------------------------------------------------------------
    Energy_FEM_quad = [];
    Energy_VEM_quad = [];
    Energy_VEM_voro = [];
    
    %% Initialize arrays for max displacements
    md_FEM_quad = [];
    md_VEM_quad = [];
    md_VEM_voro = [];
    
    % ---------------------------------------------------------------------
    % 6.4) Loop over the three mesh refinements (n = 1:3)
    % ---------------------------------------------------------------------
    for n = 1:3
        fprintf('\n-------------------\n');
        fprintf('Processing Case %d for tau = %f\n', n, current_tau);
        
        %% FEM - Quadrilateral
        fprintf('  -> FEM-quad-%d\n', n);
        sol_FEM_quad = linear_SOLVER_FEM( ...
            eval(sprintf('geom_FEM_quad_%d', n)), ...
            mat, ...
            eval(sprintf('BC%d_quad', n)), ...
            options.p, "gauss");
        
        %% VEM - Quadrilateral
        fprintf('  -> VEM-quad-%d (tau = %f)\n', n, current_tau);
        sol_VEM_quad = linear_SOLVER_VEM( ...
            eval(sprintf('geom_VEM_quad_%d', n)), ...
            mat, ...
            eval(sprintf('BC%d_quad', n)), ...
            options.p, "gauss", current_tau);
        
        %% VEM - Voronoi
        fprintf('  -> VEM-voro-%d (tau = %f)\n', n, current_tau);
        sol_VEM_voro = linear_SOLVER_VEM( ...
            eval(sprintf('geom_VEM_voro_%d', n)), ...
            mat, ...
            eval(sprintf('BC%d_voro', n)), ...
            options.p, "gauss", current_tau);
        
        %% Post-Processing
        % -- FEM Quadrilateral
        POSTPRO_FEM( ...
            eval(sprintf('geom_VEM_quad_%d', n)), ...
            eval(sprintf('geom_FEM_quad_%d', n)), ...
            sol_FEM_quad);
        title(sprintf('FEM-quad-%d', n), 'Interpreter','none');
        saveas(gcf, fullfile(resultsFolder, sprintf('Postpro_FEM_quad_%d.png', n)));
        
        % -- VEM Quadrilateral
        POSTPRO( ...
            eval(sprintf('geom_VEM_quad_%d', n)), ...
            sol_VEM_quad);
        title(sprintf('VEM-quad-%d', n), 'Interpreter','none');
        saveas(gcf, fullfile(resultsFolder, sprintf('Postpro_VEM_quad_%d.png', n)));
        
        % -- VEM Voronoi
        POSTPRO( ...
            eval(sprintf('geom_VEM_voro_%d', n)), ...
            sol_VEM_voro);
        title(sprintf('VEM-voro-%d', n), 'Interpreter','none');
        saveas(gcf, fullfile(resultsFolder, sprintf('Postpro_VEM_voro_%d.png', n)));
        
        pause(1); % Optional pause
        
        % -----------------------------------------------------------------
        % Store Energy Data
        % -----------------------------------------------------------------
        Energy_FEM_quad = [Energy_FEM_quad, sol_FEM_quad.Energy];
        Energy_VEM_quad = [Energy_VEM_quad, sol_VEM_quad.Energy];
        Energy_VEM_voro = [Energy_VEM_voro, sol_VEM_voro.Energy];

        md_FEM_quad = [md_FEM_quad, max(abs(sol_FEM_quad.UU))];
        md_VEM_quad = [md_VEM_quad, max(abs(sol_VEM_quad.UU))];
        md_VEM_voro = [md_VEM_voro, max(abs(sol_VEM_voro.UU))];
    end
    
    % ---------------------------------------------------------------------
    % 6.5) Plot Energy Comparison (Including Optional Analytical Lines)
    % ---------------------------------------------------------------------
    % Un = 50.131;  % Example reference
    Un = 136.71;
    figure;
    hold on;
    plot(abs(Energy_FEM_quad), '-ob', 'LineWidth', 1.5, ...
         'MarkerSize', 6, 'MarkerFaceColor', 'k');
    plot(abs(Energy_VEM_quad), '-or', 'LineWidth', 1.5, ...
         'MarkerSize', 6, 'MarkerFaceColor', 'r');
    plot(abs(Energy_VEM_voro), '-og', 'LineWidth', 1.5, ...
         'MarkerSize', 6, 'MarkerFaceColor', 'g');
    plot([Un, Un, Un], '-k', 'LineWidth', 1.5);
    
    % If you want to plot the analytical lines for each sample:
    % plot(abs([Ua Ua Ua]), '-k', 'LineWidth', 1.5);
    % plot(abs([Un Un Un]), '-k', 'LineWidth', 1.5);

    grid on;
    title('Energy Comparison', 'FontSize', 20, 'Interpreter', 'latex');
    legend('FEM-quad', 'VEM-quad', 'VEM-voro', 'Analytical', ...
           'FontSize', 14, 'Interpreter', 'latex', 'Location', 'best');
    xlabel('Sample Index', 'FontSize', 16, 'Interpreter', 'latex');
    ylabel('Energy Value', 'FontSize', 16, 'Interpreter', 'latex');
    set(gca, 'FontSize', 14);
    box on;
    
    % Save the comparative energy plot
    saveas(gcf, fullfile(resultsFolder, 'Energy_Comparison.png'));

    %% Plot Energy Error
    figure;
    hold on;
    plot(log10(abs(abs(Energy_FEM_quad)-Un)), '-ob', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'k');
    plot(log10(abs(abs(Energy_VEM_quad)-Un)), '-or', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    plot(log10(abs(abs(Energy_VEM_voro)-Un)), '-og', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'g');
    % Plot reference energies (replicated for each case)
    grid on;
    title('Energy Error', 'FontSize', 20, 'Interpreter', 'latex');
    legend('FEM-quad', 'VEM-quad', 'VEM-voro', 'Analytical (Ua)', 'Analytical (Un)', ...
           'FontSize', 14, 'Interpreter', 'latex', 'Location', 'best');
    xlabel('Sample Index', 'FontSize', 16, 'Interpreter', 'latex');
    ylabel('Energy Value', 'FontSize', 16, 'Interpreter', 'latex');
    set(gca, 'FontSize', 14);
    box on;
    saveas(gcf, fullfile(resultsFolder, 'Energy_Error.png'));

    %% Plot Displacement Comparisons
    % mdn = 0.11992;
    mdn = 0.29314;
    figure;
    hold on;
    plot(abs(md_FEM_quad), '-ob', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'k');
    plot(abs(md_VEM_quad), '-or', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    plot(abs(md_VEM_voro), '-og', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'g');
    % Plot reference energies (replicated for each case)
    plot([mdn, mdn, mdn], '-k', 'LineWidth', 1.5);
    grid on;
    title('Max Displacement Comparison', 'FontSize', 20, 'Interpreter', 'latex');
    legend('FEM-quad', 'VEM-quad', 'VEM-voro', 'Analytical (Ua)', 'Analytical (Un)', ...
           'FontSize', 14, 'Interpreter', 'latex', 'Location', 'best');
    xlabel('Sample Index', 'FontSize', 16, 'Interpreter', 'latex');
    ylabel('Energy Value', 'FontSize', 16, 'Interpreter', 'latex');
    set(gca, 'FontSize', 14);
    box on;
    saveas(gcf, fullfile(resultsFolder, 'Max_Displacement_Comparison.png'));


    fprintf('[SUPPORTED BEAM]\n');
    % Close the diary (log file) for this tau
    diary off;
end

fprintf('\n>>> All simulations have been completed successfully.\n');
fprintf('[SUPPORTED BEAM]\n');
close all
