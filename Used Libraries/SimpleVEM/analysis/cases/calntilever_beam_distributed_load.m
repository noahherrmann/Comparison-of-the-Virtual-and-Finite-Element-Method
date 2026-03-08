%cantilever beam - distributed load

close all
clear
clc

% Suppress all warnings
warning('off', 'all')

%% Prompt for Output Folder
outputFolder = uigetdir('', 'Select Output Folder for Saving Results');
if outputFolder == 0
    error('No folder selected. Exiting the script.');
end

%% Parameters and Tau values
p = 1;
tau_opt = 0.13;
tau_values = [1, tau_opt]; % Run simulation for these two tau values

% Create subfolders for each tau value
tau_folder_names = cell(size(tau_values));
for i = 1:length(tau_values)
    tau_str = num2str(tau_values(i));
    % Replace dot with underscore for folder naming if necessary
    tau_str = strrep(tau_str, '.', '_');
    subfolder = fullfile(outputFolder, ['tau_' tau_str]);
    if ~exist(subfolder, 'dir')
        mkdir(subfolder);
    end
    tau_folder_names{i} = subfolder;
end

%% Create Base Geometries
[~, contour, options] = AllQuadrilateralMesh(15, 2, 5, 3);
geometry = contour;

nx = 21;
ny = 7;

% Quadrilateral Geometries for FEM and VEM
[geom1, ~] = AllQuadrilateralMesh(15, 2, nx, ny);
[geom2, ~] = AllQuadrilateralMesh(15, 2, 2*nx, 2*ny);
[geom3, ~] = AllQuadrilateralMesh(15, 2, 3*nx, 3*ny);

% FEM geometries (remove first column of connectivity for FEM solver)
geom_FEM_quad_1 = geom1; 
geom_FEM_quad_2 = geom2; 
geom_FEM_quad_3 = geom3;
geom_FEM_quad_1.conn = geom_FEM_quad_1.conn(:, 2:end);
geom_FEM_quad_2.conn = geom_FEM_quad_2.conn(:, 2:end);
geom_FEM_quad_3.conn = geom_FEM_quad_3.conn(:, 2:end);

% VEM quadrilateral geometries (same as quadrilateral meshes)
geom_VEM_quad_1 = geom1;
geom_VEM_quad_2 = geom2;
geom_VEM_quad_3 = geom3;

%% Voronoi Geometries for VEM
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
% geom_VEM_voro_1 = voronoiDiscretizationRGP(geometry, nx * ny);
% geom_VEM_voro_2 = voronoiDiscretizationRGP(geometry, 4 * nx * ny);
% geom_VEM_voro_3 = voronoiDiscretizationRGP(geometry, 9 * nx * ny);

%% Check Orientations of Voronoi Meshes
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

%% Set Options
options.p = p;
options.problemtype = 3;

%% Material Properties
mat = linearElastic2D_Steel;

fprintf('[CANTILEVER BEAM - DISTRIBUTED LOAD]\n');

%% Loop over tau values
for tau_idx = 1:length(tau_values)
    current_tau = tau_values(tau_idx);
    resultsFolder = tau_folder_names{tau_idx};

    % Start logging to a diary file in the results subfolder
    diaryFile = fullfile(resultsFolder, 'simulation_log.txt');
    diary(diaryFile);
    diary on;
  
    fprintf('>>> Starting simulation for tau = %f\n\n', current_tau);
    
    %% Plot Pre-processing Figures for Quadrilateral and Voronoi meshes
    % Quadrilateral meshes pre-processing plot
    Plot_MeshFEM(geom_VEM_quad_1.conn, geom_VEM_quad_1.coord, 1, '-', 0);
    [BC1_quad.diric, BC1_quad.neum_c, BC1_quad.neum_d] = normalBC(geom_VEM_quad_1.coord, geom_VEM_quad_1.conn, options);
    % Save the quadrilateral pre-processing plot
    saveas(gcf, fullfile(resultsFolder, 'Prepro_Quad_1.png'));

    Plot_MeshFEM(geom_VEM_quad_2.conn, geom_VEM_quad_2.coord, 1, '-', 0);
    [BC2_quad.diric, BC2_quad.neum_c, BC2_quad.neum_d] = normalBC(geom_VEM_quad_2.coord, geom_VEM_quad_2.conn, options);
    % Save the quadrilateral pre-processing plot
    saveas(gcf, fullfile(resultsFolder, 'Prepro_Quad_2.png'));
    
    Plot_MeshFEM(geom_VEM_quad_3.conn, geom_VEM_quad_3.coord, 1, '-', 0);
    [BC3_quad.diric, BC3_quad.neum_c, BC3_quad.neum_d] = normalBC(geom_VEM_quad_3.coord, geom_VEM_quad_3.conn, options);
    % Save the quadrilateral pre-processing plot
    saveas(gcf, fullfile(resultsFolder, 'Prepro_Quad_3.png'));
    

    % Voronoi meshes pre-processing plot
    Plot_MeshFEM(geom_VEM_voro_1.conn, geom_VEM_voro_1.coord, 1, '-', 0);
    [BC1_voro.diric, BC1_voro.neum_c, BC1_voro.neum_d] = normalBC(geom_VEM_voro_1.coord, geom_VEM_voro_1.conn, options);
    % Save the Voronoi pre-processing plot
    saveas(gcf, fullfile(resultsFolder, 'Prepro_Voro_1.png'));

    Plot_MeshFEM(geom_VEM_voro_2.conn, geom_VEM_voro_2.coord, 1, '-', 0);
    [BC2_voro.diric, BC2_voro.neum_c, BC2_voro.neum_d] = normalBC(geom_VEM_voro_2.coord, geom_VEM_voro_2.conn, options);
    % Save the Voronoi pre-processing plot
    saveas(gcf, fullfile(resultsFolder, 'Prepro_Voro_2.png'));
    
    Plot_MeshFEM(geom_VEM_voro_3.conn, geom_VEM_voro_3.coord, 1, '-', 0);
    [BC3_voro.diric, BC3_voro.neum_c, BC3_voro.neum_d] = normalBC(geom_VEM_voro_3.coord, geom_VEM_voro_3.conn, options);
    % Save the Voronoi pre-processing plot
    saveas(gcf, fullfile(resultsFolder, 'Prepro_Voro_3.png'));


    %% Initialize arrays for energy results
    Energy_FEM_quad = [];
    Energy_VEM_quad = [];
    Energy_VEM_voro = [];

    %% Initialize arrays for DOFs
    DOFs_FEM_quad = [];
    DOFs_VEM_quad = [];
    DOFs_VEM_voro = [];

    %% Initialize arrays for max displacements
    md_FEM_quad = [];
    md_VEM_quad = [];
    md_VEM_voro = [];
    
    %% Loop over three test cases (different mesh resolutions)
    for n = 1:3
        fprintf('\n-------------------\n');
        fprintf('Processing Case %d for tau = %f\n', n, current_tau);
                
        %% FEM - Quadrilateral (solver works only for quadrilateral meshes)
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
        % Post-processing for FEM Quadrilateral
        POSTPRO_FEM( ...
            eval(sprintf('geom_VEM_quad_%d', n)), ...
            eval(sprintf('geom_FEM_quad_%d', n)), ...
            sol_FEM_quad);
        title(sprintf('FEM-quad-%d', n));
        saveas(gcf, fullfile(resultsFolder, sprintf('Postpro_FEM_quad_%d.png', n)));
        
        % Post-processing for VEM Quadrilateral
        POSTPRO( ...
            eval(sprintf('geom_VEM_quad_%d', n)), ...
            sol_VEM_quad);
        title(sprintf('VEM-quad-%d', n));
        saveas(gcf, fullfile(resultsFolder, sprintf('Postpro_VEM_quad_%d.png', n)));
        
        % Post-processing for VEM Voronoi
        POSTPRO( ...
            eval(sprintf('geom_VEM_voro_%d', n)), ...
            sol_VEM_voro);
        title(sprintf('VEM-voro-%d', n));
        saveas(gcf, fullfile(resultsFolder, sprintf('Postpro_VEM_voro_%d.png', n)));
        
        pause(1); % Optional pause between cases
        
        % Gather energy data
        Energy_FEM_quad = [Energy_FEM_quad, sol_FEM_quad.Energy];
        Energy_VEM_quad = [Energy_VEM_quad, sol_VEM_quad.Energy];
        Energy_VEM_voro = [Energy_VEM_voro, sol_VEM_voro.Energy];

        % Gather DOF data
        DOFs_FEM_quad = [DOFs_FEM_quad, length(sol_FEM_quad.UU)];
        DOFs_VEM_quad = [DOFs_VEM_quad, length(sol_VEM_quad.UU)];
        DOFs_VEM_voro = [DOFs_VEM_voro, length(sol_VEM_voro.UU)];

        md_FEM_quad = [md_FEM_quad, max(abs(sol_FEM_quad.UU))];
        md_VEM_quad = [md_VEM_quad, max(abs(sol_VEM_quad.UU))];
        md_VEM_voro = [md_VEM_voro, max(abs(sol_VEM_voro.UU))];
    end
    
    %% Plot Energy Comparisons
    Un = 14.664;
    figure;
    hold on;
    plot(abs(Energy_FEM_quad), '-ob', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'k');
    plot(abs(Energy_VEM_quad), '-or', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    plot(abs(Energy_VEM_voro), '-og', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'g');
    % Plot reference energies (replicated for each case)
    plot([Un, Un, Un], '-k', 'LineWidth', 1.5);
    grid on;
    title('Energy Comparison', 'FontSize', 20, 'Interpreter', 'latex');
    legend('FEM-quad', 'VEM-quad', 'VEM-voro', 'Analytical (Ua)', 'Analytical (Un)', ...
           'FontSize', 14, 'Interpreter', 'latex', 'Location', 'best');
    xlabel('Sample Index', 'FontSize', 16, 'Interpreter', 'latex');
    ylabel('Energy Value', 'FontSize', 16, 'Interpreter', 'latex');
    set(gca, 'FontSize', 14);
    box on;
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
    mdn = 0.48464;
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



    fprintf('[CANTILEVER BEAM - DISTRIBUTED LOAD]\n');
    diary off;

    %% Plot Error reduction per DOF (Convergence Plot)
    figure;
    err_FEM = abs(abs(Energy_FEM_quad) - Un);
    err_VEM_q = abs(abs(Energy_VEM_quad) - Un);
    err_VEM_v = abs(abs(Energy_VEM_voro) - Un);
    
    loglog(DOFs_FEM_quad, err_FEM, '-ob', 'LineWidth', 2, 'MarkerFaceColor', 'k'); hold on;
    loglog(DOFs_VEM_quad, err_VEM_q, '-or', 'LineWidth', 2, 'MarkerFaceColor', 'r');
    loglog(DOFs_VEM_voro, err_VEM_v, '-og', 'LineWidth', 2, 'MarkerFaceColor', 'g');
    
    grid on; grid minor;
    title(['Error vs DOFs ($\tau = ' num2str(current_tau) '$)'], 'FontSize', 18, 'Interpreter', 'latex');
    legend('FEM-quad', 'VEM-quad', 'VEM-voro', 'FontSize', 14, 'Interpreter', 'latex', 'Location', 'southwest');
    xlabel('Number of Degrees of Freedom (DOFs)', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('Absolute Energy Error $|U_{sim} - U_n|$', 'FontSize', 14, 'Interpreter', 'latex');
    set(gca, 'FontSize', 14);
    saveas(gcf, fullfile(resultsFolder, 'Error_vs_DOFs_LogLog.png'));
end

fprintf('\n>>> All simulations have been completed successfully.\n');
close all
fprintf('[CANTILEVER BEAM - DISTRIBUTED LOAD]\n');
