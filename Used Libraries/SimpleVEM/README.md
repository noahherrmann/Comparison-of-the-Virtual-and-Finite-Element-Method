# mySimpleVEM

![Testo alternativo](imgs/abstract.png)

 This MATLAB project is designed to perform 2D structural analysis using both the Finite Element Method (FEM) and the Virtual Element Method (VEM) using both structural and voronoi meshes.

## Project Overview

The simulation script carries out the following main tasks:

- **Mesh Generation:** Creates base geometries and generates multiple mesh resolutions using both quadrilateral and Voronoi discretizations.
- **Pre-Processing:** Sets up the problem by applying appropriate boundary conditions and checking mesh orientations.
- **Simulation Execution:** Runs simulations with different values of the stabilization parameter `tau` for both FEM and VEM solvers.
- **Post-Processing:** Compares simulation results by plotting energy, energy error, and maximum displacement for each mesh resolution.
- **Output Management:** Saves plots, logs, and simulation data in designated output folders for further analysis.

## Prerequisites

- **MATLAB:** Ensure you have a compatible MATLAB version installed.
- **Dependencies:** The following custom functions must be available in your MATLAB path:
  - `Curve Fitting Toolbox`
  - `Optimization Toolbox`

## How to Use the Simulation Script

1. **Open the Script:**
   
   - Launch MATLAB and run `main.m` to load the entire project.
   - Move to paper_analysis/cases, choose one of the 6 available cases (e.g. `calntilever_beam_distributed_load`) and run it

2. **Select Output Folder:**
   
   - Upon running the script, a dialog box will prompt you to select an output folder where all results (plots, logs, etc.) will be saved.
   - If no folder is selected, the script will terminate with an error message.

3. **Set Simulation Parameters:**
   
   - The script defines key parameters such as the stabilization parameter `tau`. Two values of `tau` (e.g., `1` and an optimal value like `0.13`) are used, and separate output subfolders are created for each value.
   - You can modify these values as needed for your specific analysis.

4. **Generate Meshes:**
   
   - The base geometry is generated and used to create different mesh resolutions:
     - **FEM Meshes:** Quadrilateral meshes are adjusted for the FEM solver.
     - **VEM Meshes:** Both quadrilateral and Voronoi meshes are generated for the VEM solver.
   - The script includes a check for optimal voronoi element number and correct orientation of the Voronoi meshes. If the orientation check fails, the simulation is halted.

5. **Run Simulations:**
   
   - The simulation loop iterates over each `tau` value and for three different mesh resolutions:
     - **FEM Quadrilateral:** The FEM solver is applied to the quadrilateral mesh.
     - **VEM Quadrilateral:** The VEM solver is applied to the same quadrilateral mesh.
     - **VEM Voronoi:** The VEM solver is applied to the Voronoi mesh.
   - For each simulation, energy values and maximum displacements are computed and stored.

6. **Post-Processing and Results:**
   
   - After each simulation, the script generates pre-processing and post-processing plots:
     - Mesh plots with boundary conditions.
     - Comparative plots for energy, energy error (using a reference analytical value), and maximum displacement.
   - All generated figures and logs are automatically saved in the corresponding `tau` subfolder.

7. **Logging:**
   
   - A log file (`simulation_log.txt`) is created in each output subfolder to record the simulation progress and details for each `tau` value.

## Customization and Adaptation

- **Geometry and Meshes:**  
  Modify the functions and parameters related to mesh generation if you wish to analyze a different geometry or use alternative discretizations.

- **Boundary Conditions:**  
  The script uses the `normalBC` function to set up boundary conditions. To adapt the simulation to different loading or support conditions, update this function or replace it with a custom implementation.

- **Material Properties:**  
  The material properties are defined via the `linearElastic2D_Steel` function. Change or extend this function to simulate other materials. 