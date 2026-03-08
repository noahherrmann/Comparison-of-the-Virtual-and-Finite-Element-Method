function varargout = Plot_skew(conn, pnts, color)

patch('Faces',conn(:,2:end),'Vertices',pnts,...
        'FaceColor', color,'EdgeColor','k', 'LineStyle', '-', 'LineWidth', 1);

axis equal;
end
