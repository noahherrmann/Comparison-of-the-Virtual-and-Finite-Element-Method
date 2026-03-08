function C = constitutive_eq(mat)

young = mat.young;
poiss = mat.poiss;
type = mat.type;

if type == "plane-stress"
    C = young/(1-poiss^2)*[1 poiss 0 ; poiss 1 0 ; 0 0 (1-poiss)/2];
end

if type == "plane-strain"
    
end

if type == "axisymm"

end

end