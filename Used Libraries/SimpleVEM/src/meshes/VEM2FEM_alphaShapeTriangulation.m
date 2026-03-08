function [allTriangles] = VEM2FEM_alphaShapeTriangulation(geom_VEM)

coord = geom_VEM.coord;
conn = geom_VEM.conn;

nel = size(conn,1);
allTriangles = [];

for ielem = 1:nel
   nedgesEl = conn(ielem,1);
   connEl = conn(ielem,2:nedgesEl+1);
    [~,triangles] = triangulateNonConvexPolygon(coord,connEl);
    allTriangles = [allTriangles;triangles]; %#ok<AGROW>
end

allTriangles = [allTriangles,allTriangles(:,1)];

end