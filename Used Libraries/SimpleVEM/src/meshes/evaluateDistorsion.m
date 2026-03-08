function [Distorsion,overDistortedElements] = evaluateDistorsion(coord,conn)

DEBUG = 0;

indL = [];
indA = [];
indAng = [];
nel = size(conn,1);
lengthDistorsion = zeros(nel,1);
areaDistorsion = zeros(nel,1);
angDistorsion = zeros(nel,1);
for j = 1:nel
    nnodeEl = conn(j,1);
    connEl = conn(j,2:nnodeEl+1);

    if DEBUG
        coord(connEl,:)
        figure();
        plot([coord(connEl,1);coord(connEl(1),1)],[coord(connEl,2);coord(connEl(1),2)],'LineWidth',2)
        axis equal
    end
    
    % DISTORSIONE CALCOLATO RISPETTO ALLA LUNGHEZZA DEI LATI E ALL'AREA DEL
    % POLIGONO
    AEl = 0;
    LEl = zeros(nnodeEl,1);
    for i = 1:nnodeEl
        x1 = coord(connEl(i),1);
        y1 = coord(connEl(i), 2);
        x2 = coord(connEl(mod(i, nnodeEl) + 1), 1);  % Prossimo vertice (ciclico)
        y2 = coord(connEl(mod(i, nnodeEl) + 1), 2);
        AEl = AEl + (x1 * y2 - y1 * x2);

        p1 = [x1,y1];
        p2 = [x2,y2];
        LEl(i) = norm(p1-p2);
    end
    AEl = abs(AEl)/2;

    % DISTORSIONE RISPETTO AL LATO
    % R: del cerchio circoscritto al poligono regolare
    R = sqrt(2*AEl/(nnodeEl * sin(2*pi/nnodeEl)));
    % L: lato del poligono regolare
    L = 2*R*sin(pi/nnodeEl);

    lengthDistorsion(j) = mean(abs(LEl-L))/L;

    % DISTORSIONE RISPETTO ALL'AREA 
    L = mean(LEl);
    A = nnodeEl*L^2/(4*tan(pi/nnodeEl));

    areaDistorsion(j) = abs(AEl-A)/A;

    if lengthDistorsion(j) > 1
        indL = [indL;j];
    elseif areaDistorsion(j) > 1
        indA = [indA;j];
    end

    % DISTORSIONE CALCOLATA RISPETTO AGLI ANGOLI INTERNI
    ang = zeros(nnodeEl,1);
    for i = 1:nnodeEl
        p1 = coord(connEl(mod(i-2, nnodeEl) + 1), :); % Punto precedente
        p2 = coord(connEl(i), :);                     % Punto corrente
        p3 = coord(connEl(mod(i, nnodeEl) + 1), :);   % Punto successivo
        v1 = p1 - p2;
        v2 = p3 - p2;

        % Calcolo dell'angolo tra v1 e v2 (prodotto vettoriale e scalare in 2D)
        crossProd = v1(1)*v2(2) - v1(2)*v2(1); 
        dotProd = v1(1)*v2(1) + v1(2)*v2(2);  
        ang(i) = atan2d(abs(crossProd), dotProd);
    end

    nominalValue = (nnodeEl-2)*180/nnodeEl;
    angDistorsion(j) = mean(abs(ang - nominalValue))/nominalValue;
    
    if angDistorsion(j) > 1
        indAng = [indAng;j];
    end
end


Distorsion = [lengthDistorsion,areaDistorsion,angDistorsion];
overDistortedElements.lengthMethod = indL;
overDistortedElements.areaMethod = indA;
overDistortedElements.angleMethod = indAng;

end