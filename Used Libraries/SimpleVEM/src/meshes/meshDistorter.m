function [newPoints] = meshDistorter(Points,options)

ds = 1e-3;
ind = options.problemtype;

if ind == 1 || ind == 2 || ind == 3 || ind == 4
    % ASTA, MENSOLA E TRAVE APPOGGIATA
    l = options.data(1);
    h = options.data(2);
    polygon = [ds,ds;l-ds,ds;l-ds,h-ds;ds,h-ds];
    in = findInside(polygon(:,1),polygon(:,2),Points);
    in_ind = find(in==1);
    % plot(polygon(:,1),polygon(:,2))
    % plotNodes([Points(in==1,1),Points(in==1,2)],'red')
    nx = sum(Points(:,2)==0);
    ny = size(Points,1)/nx;
    newPoints = Points;

    for i = 1:length(in_ind)
        ind = in_ind(i);
        dist(:,1) = abs([Points(ind,1)-Points(ind-nx,1);Points(ind,1)-Points(ind+1,1);Points(ind,1)-Points(ind+nx,1);Points(ind,1)-Points(ind-1,1)]);
        dist(:,2) = abs([Points(ind,2)-Points(ind-nx,2);Points(ind,2)-Points(ind+1,2);Points(ind,2)-Points(ind+nx,2);Points(ind,2)-Points(ind-1,2)]);
        distX = dist(:,1);
        distY = dist(:,2);
        distX = distX(distX~=0);
        distY = distY(distY~=0);
        spost_maxX = min(distX);
        spost_maxY = min(distY);
        spost_max = [spost_maxX,spost_maxY];
        newPoints(ind,:) = rand(1,2).*spost_max + Points(ind,:);
    end

elseif ind == 5
    % TUBO IN PRESSIONE 
    r_int = options.data(1);
    r_ext = options.data(2);
    th1 = deg2rad(linspace(1,89,100))';
    th2 = flip(th1);
    polygon = [(r_int+ds)*cos(th1),(r_int+ds)*sin(th1);(r_ext-ds)*cos(th2),(r_ext-ds)*sin(th2)];
    in = findInside(polygon(:,1),polygon(:,2),Points);
    in_ind = find(in==1);
    % plot(polygon(:,1),polygon(:,2))
    % plotNodes([Points(in==1,1),Points(in==1,2)],'red')
    nr = sum(Points(:,1)<0.0001);
    nth = size(Points,1)/nr;
    newPoints = Points;

    for i = 1:length(in_ind)
        ind = in_ind(i);
        dist(:,1) = abs([Points(ind,1)-Points(ind-nth,1);Points(ind,1)-Points(ind+1,1);Points(ind,1)-Points(ind+nth,1);Points(ind,1)-Points(ind-1,1)]);
        dist(:,2) = abs([Points(ind,2)-Points(ind-nth,2);Points(ind,2)-Points(ind+1,2);Points(ind,2)-Points(ind+nth,2);Points(ind,2)-Points(ind-1,2)]);
        distX = dist(:,1);
        distY = dist(:,2);
        distX = distX(distX~=0);
        distY = distY(distY~=0);
        spost_maxX = min(distX);
        spost_maxY = min(distY);
        spost_max = [spost_maxX,spost_maxY];
        newPoints(ind,:) = rand(1,2).*spost_max + Points(ind,:);
    end

elseif ind == 6

    % PIASTRA FORATA
    l = options.data(1);
    h = options.data(2);
    r = options.data(3);
    th = deg2rad(linspace(1,89,100))';
    polygon = [(r+ds)*cos(th),(r+ds)*sin(th);(r+ds)*cos(th(end)),h-ds;l-ds,h-ds;l-ds,(r+ds)*sin(th(1))];
    in = findInside(polygon(:,1),polygon(:,2),Points);
    in_ind = find(in==1);
    % plot(polygon(:,1),polygon(:,2))
    % plotNodes([Points(in==1,1),Points(in==1,2)],'red')
    nr = sum(Points(:,1)<0.0001);
    nx = sum(Points(:,2)<0.0001)-nr+1;
    ny = sum(Points(:,1)==l);
    nth = (size(Points,1)-ny*(nx-1))/nr;
    newPoints = Points;

    for i = 1:length(in_ind)
        ind = in_ind(i);
        if ind<(nr*(nth-1))
            dist(:,1) = abs([Points(ind,1)-Points(ind-nth,1);Points(ind,1)-Points(ind+1,1);Points(ind,1)-Points(ind+nth,1);Points(ind,1)-Points(ind-1,1)]);
            dist(:,2) = abs([Points(ind,2)-Points(ind-nth,2);Points(ind,2)-Points(ind+1,2);Points(ind,2)-Points(ind+nth,2);Points(ind,2)-Points(ind-1,2)]);
        elseif ind>(nr*(nth-1)) && ind<(nr*nth)
            dist(:,1) = abs([Points(ind,1)-Points(ind+nth,1);Points(ind,1)-Points(ind+1,1);Points(ind,1)-Points(ind+ny,1);Points(ind,1)-Points(ind-1,1)]);
            dist(:,2) = abs([Points(ind,2)-Points(ind-nth,2);Points(ind,2)-Points(ind+1,2);Points(ind,2)-Points(ind+ny,2);Points(ind,2)-Points(ind-1,2)]);
        else
            dist(:,1) = abs([Points(ind,1)-Points(ind+ny,1);Points(ind,1)-Points(ind+1,1);Points(ind,1)-Points(ind+ny,1);Points(ind,1)-Points(ind-1,1)]);
            dist(:,2) = abs([Points(ind,2)-Points(ind-ny,2);Points(ind,2)-Points(ind+1,2);Points(ind,2)-Points(ind+ny,2);Points(ind,2)-Points(ind-1,2)]);
        end
        distX = dist(:,1);
        distY = dist(:,2);
        distX = distX(distX~=0);
        distY = distY(distY~=0);
        spost_maxX = min(distX);
        spost_maxY = min(distY);
        spost_max = [spost_maxX,spost_maxY];
        newPoints(ind,:) = rand(1,2).*spost_max + Points(ind,:);
    end

else 

    % MEMBRANA DI COOK
    H1 = options.data(1);
    H2 = options.data(2);
    H3 = options.data(3);
    L = options.data(4);
    polygon = [ds,ds;L-ds,H1+H2-H3+ds;L-ds,H1+H2-ds;ds,H1-ds];
    in = findInside(polygon(:,1),polygon(:,2),Points);
    in_ind = find(in==1);
    % plot(polygon(:,1),polygon(:,2))
    % plotNodes([Points(in==1,1),Points(in==1,2)],'red')
    ny = sum(Points(:,1)==0);
    nx = size(Points,1)/ny;
    newPoints = Points;

    for i = 1:length(in_ind)
        ind = in_ind(i);
        dist(:,1) = abs([Points(ind,1)-Points(ind-nx,1);Points(ind,1)-Points(ind+1,1);Points(ind,1)-Points(ind+nx,1);Points(ind,1)-Points(ind-1,1)]);
        dist(:,2) = abs([Points(ind,2)-Points(ind-nx,2);Points(ind,2)-Points(ind+1,2);Points(ind,2)-Points(ind+nx,2);Points(ind,2)-Points(ind-1,2)]);
        distX = dist(:,1);
        distY = dist(:,2);
        distX = distX(distX~=0);
        distY = distY(distY~=0);
        spost_maxX = min(distX);
        spost_maxY = min(distY);
        spost_max = [spost_maxX,spost_maxY];
        newPoints(ind,:) = rand(1,2).*spost_max + Points(ind,:);
    end
end

end