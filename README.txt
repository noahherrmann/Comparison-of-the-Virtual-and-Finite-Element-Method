A COMPARISON OF THE VIRTUAL AND FINITE ELEMENT METHOD

This project conducts an extensive numerical comparison between the 
Finite Element Method (FEM) and the Virtual Element Method (VEM).
The study focuses on robustness and convergence rates, but most importantly,
it focuses on the error reduction per degree of freedom (DOF).
This work was developed as part of a Bachelor's Thesis (2026).

------------------------------------------------------------------------
DIRECTORY STRUCTURE
------------------------------------------------------------------------

/Comparison of the Virtual and Finite Element Method
│
├── Used Libraries/         % External frameworks and toolboxes
│   ├── ifem/               % iFEM Toolbox (Reference FEM solutions)
│   ├── vem2D/              % VEM library for higher-order (k > 1) tests
│   └── SimpleVEM/          % mySimpleVEM code from Santi et al. (2025)
│
├── mesh/                   % Directory for .mat mesh files
│
└── [Main Scripts]          % Analysis scripts

------------------------------------------------------------------------
MAIN SCRIPTS DESCRIPTION
------------------------------------------------------------------------

1. Convergence_Analysis.m
   A standard test for calculating L2 and H1 error norms and
   determining numerical convergence rates on Voronoi and hexagonal meshes.

2. Convergence_of_VEM_on_LDomain.m
   A standard test for calculating L2 and H1 error norms and
   determining numerical convergence rates on L-shaped domains on triangular meshes.

3. VEMvsFEM_DOF_Voronoi.m
   Compares the L^2-Error per Degree of Freedom (DOF) between the FEM and the VEM
   on two Voronoi meshes in very smooth problems.

4. VEMvsFEM_DOF_Comparison.m
   Compares the L^2-Error per Degree of Freedom (DOF) between the FEM and the VEM
   on two quadrilateral and two Voronoi meshes in very smooth problems.

5. VEMvsFEM_DOF_Comparison_Stiff.m
   Compares the L^2-Error per Degree of Freedom (DOF) between the FEM and the VEM
   on two quadrilateral and two Voronoi meshes in stiffer problems.

6. VEMvsFEM_DOF_absolute_value_function
   Compares the L^2-Error per Degree of Freedom (DOF) between the FEM and the VEM
   on four Voronoi meshes with increasing refinement in stiffer problems.

7. VEMvsFEM_DOF_LDomain
   Compares the L^2-Error per Degree of Freedom (DOF) between the FEM and the VEM
   on an L-shaped domain consisting of four triangular meshes with increasing
   refinement in very smooth problems.

8. Gaussian_function_Peak_Test.m
   Compares the FEM and the VEM for a Gaussian function on a L-shaped domain.

9. Fireball_Time_Lapse.m
   This simulates a "fireball" and should be viewed as a fun implementation and
   use of the VEM for k=1.


------------------------------------------------------------------------
HOW TO RUN
------------------------------------------------------------------------

1. Open MATLAB and set the root directory of this project as your 
   Current Folder.
2. Run any of the Main Scripts.
3. The scripts automatically initialize the necessary paths.

------------------------------------------------------------------------
LICENSE & CREDITS
------------------------------------------------------------------------

- The 'SimpleVEM' library was created by Gian Maria Santi, Daniela Francia,
  Enrico Stacchini, Francesco Cesari and can be downloaded from their GitHub:
  https://github.com/mySimpleProjects/SimpleVEM.
- The used FEM solvers are part of the iFEM framework by Long Chen and can
  be downloaded from their GitHub: https://github.com/lyc102/ifem.
- The used VEM solvers are part of the vem2D framework by Juan G. Calvo and
  can be downloaded from their GitHub: https://github.com/jgcalvo/vem2D.

========================================================================