function plot_inversion_params(inversion,boundary)
xObs = inversion.xObs;
dataObs = inversion.dataObs;
dataPredict = inversion.dataPredict;
xGridsCenters = inversion.xGridsCentersPlot;
zGridsCenters = inversion.zGridsCenters;
mPredict = inversion.mPredict;
xBoundary = boundary(:,1:2);
zBoundary = boundary(:,3:4);
figure('Position',[500 200 500 600]);
subplot(2,1,1);
plot(xObs,dataObs,'Color','r','LineWidth',2,'LineStyle','-','Marker','*');
hold on;
plot(xObs,dataPredict,'Color','b','LineWidth',2,'LineStyle','-','Marker','o');
legend('观测场','预测场','FontWeight','Bold','FontSize',12);
title('重力场拟合','FontWeight','bold','FontSize',30);
set(gca,'FontWeight','Bold','FontSize',12);
ylabel('Gz(mGal)','FontSize',12,'FontWeight','bold','FontName','Times New Roman');
grid on;
hold off;

subplot(2,1,2);
imagesc(xGridsCenters,zGridsCenters,mPredict);
hold on;
for i = 1:size(boundary,1)
    line(xBoundary(i,:),zBoundary(i,:),'LineStyle','-.','Color','k','Linewidth',1.5);
end
title('二维密度反演模型','FontWeight','bold','FontSize',30);
xlabel('x(m)','FontSize',30,'FontWeight','bold','FontName','Times New Roman');
ylabel('z(m)','FontSize',30,'FontWeight','bold','FontName','Times New Roman');
set(gca,'FontWeight','Bold','FontSize',12);
colormap('jet');
clim([-1 1]);
colorBar = colorbar('vert','southoutside');
colorBar.Label.String = '\rho(g/cm^3)';
colorBar.Label.FontName = 'Times New Roman';
grid on;
hold off;
end