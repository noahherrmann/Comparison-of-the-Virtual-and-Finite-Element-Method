function edgeL = edgesLocal(nedges)

edgeL = zeros(nedges,2);

for i=1:nedges-1
    edgeL(i,1) = i;
    edgeL(i,2) = i+1;
end

edgeL(end,1) = nedges;
edgeL(end,2) = 1;

end