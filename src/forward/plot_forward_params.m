function plot_forward_params(forward,dataObs,mTrue,boundary)

xBoundary = boundary(:,1:2);
zBoundary = boundary(:,3:4);
xObs = forward.xObs;
xGridsCenters = forward.xGridsCenters;
zGridsCenters = forward.zGridsCenters;

figure('Position',[500 200 500 600]);
subplot(2,1,1);
plot(xObs,dataObs,'Color','r','LineWidth',2,'LineStyle','-','Marker','o');
hold on;
title('重力场','FontWeight','bold','FontSize',30);
set(gca,'FontWeight','Bold','FontSize',12);
ylabel('Gz(mGal)','FontSize',12,'FontWeight','bold','FontName','Times New Roman');
grid on;
hold off;

subplot(2,1,2);
imagesc(xGridsCenters,zGridsCenters,mTrue);
hold on;
for i = 1:size(xBoundary,1)
    line(xBoundary(i,:),zBoundary(i,:),'LineStyle','-.','Color','k','Linewidth',2);
end
title('二维密度模型','FontWeight','bold','FontSize',30);
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