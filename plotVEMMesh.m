function hmax = plotVEMMesh(meshSource)
% PLOTVEMMESH Visualizes a VEM polygonal mesh from a file or struct
% USE: plotVEMMesh('mesh/voronoi_02.mat') or plotVEMMesh(meshStruct)
    if ischar(meshSource) || isstring(meshSource)
        data = load(meshSource);
        mesh = data.mesh;
    else
        mesh = meshSource;
    end

    numElements = length(mesh.elems);
    hmax = 0;

    figure('Color', 'w'); hold on; axis equal; grid on;

    for i = 1:numElements
        v_idx = mesh.elems{i};
        verts = mesh.verts(v_idx, :);

        [h, ~, ~] = geomElement(verts);
        hmax = max(hmax, h);
        
        patch('Faces', 1:length(v_idx), 'Vertices', verts, ...
              'FaceColor', [0.9 0.9 1], 'EdgeColor', 'k');
    end

    title(['Mesh: ', num2str(numElements), ' Elements, h_{max} = ', num2str(hmax)]);
    xlabel('x'); ylabel('y');
end