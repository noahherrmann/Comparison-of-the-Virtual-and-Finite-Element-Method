function meshChecker(geom, isVEM)

if isVEM
    skewness(geom.coord, geom.conn)
else
    [conn_tris, conn_quads, conn_poly] = separate_conn(geom.conn);
    geoCheckTris(geom.coord, conn_tris)
    hold on;
    geoCheckQuads(geom.coord, conn_quads)
    hold on;
    skewness(geom.coord, conn_poly)
end
hold on;
plotNodes(geom.coord, 'black')

end