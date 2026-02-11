function plot_stablizers(inversion,boundary)
MS = inversion.minimumSupport;
New = inversion.newStablizer;
xBoundary = boundary(:,1:2);
zBoundary = boundary(:,3:4);
xGridsCenters = inversion.xGridsCentersPlot;
zGridsCenters = inversion.zGridsCenters;
figure;
subplot(2,1,1);
imagesc(xGridsCenters,zGridsCenters,MS);
hold on;
for i = 1:size(boundary,1)
    line(xBoundary(i,:),zBoundary(i,:),'LineStyle','-.','Color','k','Linewidth',1.5);
end
title('最小支撑','FontWeight','bold','FontSize',30);
xlabel('x(m)','FontSize',30,'FontWeight','bold','FontName','Times New Roman');
ylabel('z(m)','FontSize',30,'FontWeight','bold','FontName','Times New Roman');
set(gca,'FontWeight','Bold','FontSize',12);
colormap('jet');
colorbar('vert','southoutside');
grid on;
hold off;
subplot(2,1,2);
imagesc(xGridsCenters,zGridsCenters,New);
hold on;
for i = 1:size(boundary,1)
    line(xBoundary(i,:),zBoundary(i,:),'LineStyle','-.','Color','k','Linewidth',1.5);
end
title('新稳定泛函','FontWeight','bold','FontSize',30);
xlabel('x(m)','FontSize',30,'FontWeight','bold','FontName','Times New Roman');
ylabel('z(m)','FontSize',30,'FontWeight','bold','FontName','Times New Roman');
set(gca,'FontWeight','Bold','FontSize',12);
colormap('jet');
colorbar('vert','southoutside');
grid on;
hold off;
end