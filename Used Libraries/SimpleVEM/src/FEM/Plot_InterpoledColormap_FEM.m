function varargout = Plot_InterpoledColormap_FEM(connectivityMatrix,coord, value, varargin)

if nargin == 3
    h = NewPrettyFigure;
    varargout{1} = h;
else
    figure( varargin{1} );
    varargout{1} = [];
end

nnode = 4;
nelem = size(connectivityMatrix,1);

isolines = 15;
step = 5;

range = abs(max(value)-min(value));
global_min = min(value)-0.03*range;
global_max = max(value)+0.03*range;
global_levels = linspace(global_min, global_max, isolines);

for ielem=1:nelem
    [xi,eta] = meshgrid(-1:2/step:1 , -1:2/step:1);
    for i=1:length(eta)
        for j=1:length(xi)
            enne = s2424_NN([-1 -1 ; 1 -1 ; 1 1 ; -1 1],[xi(j,i),eta(j,i)]);
            for m=1:nnode
                NN(j,i,m) = enne(m);
            end
        end
    end
    VALUE = zeros(length(xi),length(eta));
    xx = zeros(length(xi),length(eta));
    yy = zeros(length(xi),length(eta));
    for i=1:size(NN,3)
        VALUE = VALUE + NN(:,:,i)*value(connectivityMatrix(ielem,i));
        xx = xx + coord(connectivityMatrix(ielem,i),1)*NN(:,:,i);
        yy = yy + coord(connectivityMatrix(ielem,i),2)*NN(:,:,i);
    end

%     minx = min(value);
%     maxx = max(value);
%     levels =  minx:(maxx-minx)/step:maxx;

    [C,h] = contourf(xx,yy,VALUE,global_levels,'EdgeColor','none');
    %     clabel(C,h) Non credo che questo sia davvero utile

end
colorbar

end
