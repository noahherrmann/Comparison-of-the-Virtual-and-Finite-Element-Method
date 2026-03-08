function plotNodes(pnts, color)
h = plot(pnts(:,1), pnts(:,2), 'o', 'MarkerSize', 5, 'MarkerFaceColor', color); % Aggiunge i nodi come punti riempiti
set(h, 'PickableParts', 'none', 'HitTest', 'off');
end