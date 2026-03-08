function edges = edgeExtract(nedges, conn)

connTmp = [conn,conn(1)];
edges = [];

for i=1:nedges
    edges = [edges; connTmp(i),connTmp(i+1)];
end

end