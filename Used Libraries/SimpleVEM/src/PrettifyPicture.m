function PrettifyPicture

FontGCA = 'Nimbus Roman No9 L';
inTex = 'latex';

% axis equal;
hold all, box on, grid on
xlabel('x','interp','latex','fontname','Nimbus Roman No9 L','fontsize',10)
ylabel('y','interp','latex','fontname','Nimbus Roman No9 L','fontsize',10)
rotate3d on
set(gca,'GridLineStyle','--')

ha1 = gca;
ha1.XAxis.FontName = FontGCA;
for i = 1 : length(ha1.YAxis)
    ha1.YAxis(i).FontName = FontGCA;
    ha1.YAxis(i).TickLabelInterpreter = inTex;
end
ha1.XAxis.TickLabelInterpreter = inTex;
ha1.Color='none'; 
end