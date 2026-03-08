function delta = distortionQL1(conn, coord)

    A = 4 * quadrangleArea(conn, coord);
    
    conn123 = [conn(1) conn(2) conn(3)];
    conn124 = [conn(1) conn(2) conn(4)];
    
    B = 4 * (triangleArea(conn123, coord) - ...
             triangleArea(conn124, coord));

    conn134 = [conn(1) conn(3) conn(4)];
    
    C = 4 * (triangleArea(conn134, coord) - ...
             triangleArea(conn124, coord));

    delta = (sqrt(B^2 + C^2)) / A;

end
