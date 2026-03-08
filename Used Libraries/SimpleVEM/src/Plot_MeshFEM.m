function varargout = Plot_MeshFEM(conn, pnts, i, linestyle, alf, varargin)

if nargin == 5
    h = NewPrettyFigure;
    varargout{1} = h;
else
    figure( varargin{1} );
    varargout{1} = [];
end

switch i
    case 1
        color = 'k';
    case 2
        color = 'b';
    case 3
        color = 'g';
    case 4
        color = 'y';
    case 5
        color = 'r';
end

patch('Faces',conn(:,2:end),'Vertices',pnts,...
        'FaceColor', [0.8 0.8 0.8],'EdgeColor',color, 'LineStyle', linestyle);

if alf
    alpha(0.)  
end
axis equal;
end
