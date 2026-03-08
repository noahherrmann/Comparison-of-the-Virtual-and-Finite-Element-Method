function varargout = plotNormals(vertices, edges, normals, varargin)

if nargin == 3
    h = NewPrettyFigure;
    varargout{1} = h;
else
    figure( varargin{1} );
    varargout{1} = [];
end

% Numero di lati
n = size(edges, 1);

% Calcolo dei punti medi e delle normali
midPoints = zeros(n, 2);

for i = 1:n
    v1 = vertices(edges(i, 1), :);
    v2 = vertices(edges(i, 2), :);
    midPoints(i, :) = (v1 + v2) / 2;
 end

quiver(midPoints(:,1), midPoints(:,2), normals(:,1), normals(:,2), 0, 'b');

end
