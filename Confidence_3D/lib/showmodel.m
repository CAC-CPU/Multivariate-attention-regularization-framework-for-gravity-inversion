function showmodel(showmodelpar)

% Show the 3-D view of the model

x = unique(showmodelpar.x);
y = unique(showmodelpar.y);
z = unique(showmodelpar.z);
xob = showmodelpar.xd;
yob = showmodelpar.yd;
topx = showmodelpar.topox;
topy = showmodelpar.topoy;
topz = showmodelpar.topoz;
comp = showmodelpar.component;
mt = showmodelpar.model;
nx = length(x); ny = length(y); nz = length(z);
dx = x(2)-x(1); dy = y(2)-y(1); dz = z(2)-z(1);
depth = z(end)-z(1)+dz;
MTS = mt(:);
color=Jetconvert(MTS);
tx = unique(topx); ty = unique(topy);
[TX TY] = meshgrid(tx,ty);
TZ=(reshape(topz,length(tx),length(ty)))';
figure(1);
view(3);
surf(TX,TY,TZ);
alpha(0.4)
colorbar;
colormap(flipud(jet));
shading flat;
hold on;
for ll=1:ny
    for mm=1:nx
        for nn=1:nz
            if MTS((ll-1)*(nx*nz)+(mm-1)*nz+nn)>0.08
                corob(x(mm),y(ll),z(nn),[dx,dy,dz],color((ll-1)*(nx*nz)+(mm-1)*nz+nn,:));
            end
        end
    end 
end
axis([min(xob) max(xob) min(yob) max(yob) min(topz) min(topz)+depth]);
grid on; 
set(gca, 'fontweight','bold');
xlabel('x.m','fontweight','bold');
ylabel('y.m','fontweight','bold');
zlabel('z.m','fontweight','bold');
set(gca, 'fontweight','bold');
set(gca,'zdir','reverse')
colormap(jet);
h=colorbar;
barLims=[min(MTS) max(MTS)];%colorbar limits
dif=diff(barLims);
Nlbl=11;
lbls=[barLims(1):dif/(Nlbl-1):barLims(2)]';
CBlim=get(h,'Ylim');
Ytcks=CBlim(1):diff(CBlim)/(Nlbl-1):CBlim(2);
set(h,'Ytick',Ytcks,'Yticklabel',num2str(lbls,3));
% set(h, 'xaxislocation','top');
set(h,'fontweight','bold')
for it =1:size(comp,1)
    temp=cutspace(comp(it,:));
    if temp(1) == 'h'
        set(get(h,'xlabel'),'string','\Delta\chi(SI)','fontsize',10,'fontweight','bold');
    elseif temp(1) == 'g'
        set(get(h,'xlabel'),'string','\Delta\rho g/cm^3','fontsize',10,'fontweight','bold');
    end
end
title('Synthetic model with topography');
axis equal;
axis([min(xob) max(xob) min(yob) max(yob) min(topz) min(topz)+depth]);
