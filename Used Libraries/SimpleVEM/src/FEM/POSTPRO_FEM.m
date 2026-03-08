function POSTPRO_FEM(geom_VEM, geom, sol)
DISPLACEMENTS = 1;
STRESSES = 0;

coord = geom.coord;

distances = pdist(coord, 'euclidean');
max_distance = max(distances);
max_U = max(abs(sol.U));

amp = (0.05*max_distance)/max_U;
% disp(['Amplification factor = ',num2str(amp)])
coordDef = coord;
coordDef(:,1) = coordDef(:,1) + sol.Ux*amp;
coordDef(:,2) = coordDef(:,2) + sol.Uy*amp;

h = NewPrettyFigure;
Plot_MeshFEM(geom_VEM.conn, coordDef, 1, '-', 1, h);

%% - - - - - - - - D I S P L A C E M E N T S - - - - - - - - -
if DISPLACEMENTS == 1

    Plot_InterpoledColormap_FEM(geom.conn, coordDef, sol.UU, h);
    PrettifyColorbar;
    axis equal, axis tight;
end
geom.conn = [linspace(1,size(geom.conn,1), size(geom.conn,1))', geom.conn];
Plot_MeshFEM(geom.conn, coordDef, 1, '-', 1, h);
geom.conn(:,1) = [];

% nomeFile = num2str(amp); % La variabile con il nome del file
% estensione = '.png'; % L'estensione del file
% nomeCompleto = [nomeFile, estensione]; % Concatenazione del nome del file e l'estensione
% 
% if options.exportGraphics
%     exportgraphics(dispPlot, nomeCompleto, 'Resolution', 300);
% end

%% - - - - - - - - S T R E S S E S - - - - - - - - -
if STRESSES == 1
    Sx = [];
    Sy = [];
    Txy = [];
    VonM = [];

    for i=1:size(elgrp_fem.conn,1)
        element = elgrp_fem.conn(i,:);
        x = coord(element,1);
        y = coord(element,2);
        UUx = Ux(element);
        UUy = Uy(element);
        [sigmaX, sigmaY, tauXY, vonMises] = stressRecovery_FEM(UUx,UUy,x,y,elgrp_fem);

        Sx = [Sx; sigmaX'];
        Sy = [Sy; sigmaY'];
        Txy = [Txy; tauXY'];
        VonM = [VonM; vonMises'];
    end
    [avrSigmaX, avrSigmaY, avrTauXY, avrVonMises] = stressAverage_FEM(size(coord,1),elgrp_fem.conn,Sx,Sy,Txy,VonM);

    sigX = PlotSigmaX(elgrp_fem,coord,elgrp_fem.conn,Ux,Uy,avrSigmaX,amplification,1);
    Plot_Sigma_NEFEM_p2(elgrp_nefem, elgrp_fem, coord, Ux,Uy, avrSigmaX, NURBS,amplification,1, sigX);

    sigy = PlotSigmaY(elgrp_fem,coord,elgrp_fem.conn,Ux,Uy,avrSigmaY,amplification,1);
    Plot_Sigma_NEFEM_p2(elgrp_nefem, elgrp_fem, coord, Ux,Uy, avrSigmaY, NURBS,amplification,1, sigy);

    taXY = PlotTauXY(elgrp_fem,coord,elgrp_fem.conn,Ux,Uy,avrTauXY,amplification,1);
    Plot_Sigma_NEFEM_p2(elgrp_nefem, elgrp_fem, coord, Ux,Uy, avrTauXY, NURBS,amplification,1, taXY);

    von = PlotVonMises(elgrp_fem,coord,elgrp_fem.conn,Ux,Uy,avrVonMises,amplification,1);
    Plot_Sigma_NEFEM_p2(elgrp_nefem, elgrp_fem, coord, Ux,Uy, avrVonMises, NURBS,amplification,1, von);
end
end