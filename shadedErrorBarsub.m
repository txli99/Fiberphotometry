function H=shadedErrorBarsub(x,m,sem,mainLineColor)
%eg. mainLineColor=[1,1,0];
%    shadedErrorBarsub(x,m,sem,mainLineColor)
H=shadedErrorBar(x,m,sem);
% mainLineColor=mainLineColor/255;
patchColor=mainLineColor+(1-mainLineColor)*0.01;
% edgeColor=mainLineColor+(1-mainLineColor)*C2;
H.mainLine.Color =mainLineColor;%中间那条线的颜色
H.mainLine.LineWidth=1;
H.patch.FaceColor  =patchColor;%填充的颜色
H.edge(1).Color='none';%上下边缘的颜色，不想要的话可以 H.edge(1).Color=‘none'
H.edge(2).Color='none';
end