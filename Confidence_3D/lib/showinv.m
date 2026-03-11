function showinv(showpar)

close all

m = showpar.ModelVector;
M = showpar.Model3d;
nx = showpar.numx;
ny = showpar.numy;
nz = showpar.numz;
x = showpar.x;
y = showpar.y;
z = showpar.z;
dx = x(2) -x(1);
dy = y(2) -y(1);
dz = z(2) -z(1);
dobs = showpar.observed;
dpre = showpar.predicted;
n = showpar.num_iteration;
misfit = showpar.percentmisfit;
apha = showpar.regpar;
MF = showpar.residual;
Pa = showpar.objectfun;
Inversion_component = showpar.component;
length_comp = size(Inversion_component,1);
xob = showpar.xobserved;
yob = showpar.yobserved;
topx = showpar.Xtopography;
topy = showpar.Ytopography;
topz = showpar.Ztopography;
n_data = length(xob);
cutoff = showpar.cutoff;
contint = 50;








x1=zeros(1,nx+1);
x1(1:nx)=x-0.5*dx;
x1(nx+1)=x(end)+0.5*dx;
z1=zeros(1,nz+1);
z1(1:nz)=z-0.5*dz;
z1(nz+1)=z(end)+0.5*dz;
[X,Y,Z]=meshgrid(x1,y,z1);
U=zeros(ny,nx+1,nz+1);


for i=1:ny
    U(i,1:nx,1:nz)=M(:,:,i)';
end

figure(1);
slice(X,Y,Z,U,(x(end)+0.5*dx)/2,(y(end)+0.5*dy)/2,(z(end)+0.5*dz)/2);
shading flat;
set(gca, 'fontweight','bold');
xlabel('x.m','fontweight','bold');
ylabel('y.m','fontweight','bold');
zlabel('z.m','fontweight','bold');
set(gca, 'fontweight','bold');
set(gca,'zdir','reverse');
colorbar;
h=colorbar('vert');
set(h, 'fontweight','bold');
axis equal;

for it =1:length_comp
    temp=cutspace(Inversion_component(it,:));
    if temp(1) == 'h'
        set(get(h,'xlabel'),'string','\Delta\chi','fontsize',10,'fontweight','bold');   % title
    elseif temp(1) == 'g'
        set(get(h,'xlabel'),'string','\Delta\rho g/cm^3','fontsize',10,'fontweight','bold');   % title
    end
end
color=Jetconvert(m);

figure(2);
view(3);
[X Z Y] = meshgrid(x,z,y);
xm = X(:); ym = Y(:); zm = Z(:);

hold on;

for it=1:length(xm)
    if m(it)>cutoff
       corob(xm(it),ym(it),zm(it),[dx,dy,dz],color(it,:));
    end
end
depth = z(end)-z(1)+dz;



% axis([0 200 0 100 -50 500]);
grid on;
%grid minor;
    
set(gca, 'fontweight','bold');
xlabel('x.m','fontweight','bold');
ylabel('y.m','fontweight','bold');
zlabel('z.m','fontweight','bold');
set(gca, 'fontweight','bold');
set(gca,'zdir','reverse')


colormap(jet);
h=colorbar;

barLims=[min(m) max(m)];
dif=diff(barLims);
Nlbl=11;
lbls=[barLims(1):dif/(Nlbl-1):barLims(2)]';
CBlim=get(h,'Ylim');
Ytcks=CBlim(1):diff(CBlim)/(Nlbl-1):CBlim(2);
set(h,'Ytick',Ytcks,'Yticklabel',num2str(lbls,3));
set(h,'fontweight','bold')
set(get(h,'xlabel'),'string','\Delta\chi','fontsize',8,'fontweight','bold');   % title
axis equal;
title({['3D body view of inversion result'];});
axis([min(xob) max(xob) min(yob) max(yob) min(topz) min(topz)+depth]);
% hold on;
% scatter3(topx,topy,topz,'.');

for it =1:length_comp
    temp=cutspace(Inversion_component(it,:));
    if temp(1) == 'h'
        set(get(h,'xlabel'),'string','\Delta\chi','fontsize',10,'fontweight','bold');   % title
    elseif temp(1) == 'g'
        set(get(h,'xlabel'),'string','\Delta\rho g/cm^3','fontsize',10,'fontweight','bold');   % title
    end
end


for it =1:length_comp
    temp=cutspace(Inversion_component(it,:));
    if temp(1) == 'h'
        set(get(h,'xlabel'),'string','\Delta\chi','fontsize',10,'fontweight','bold');   % title
    elseif temp(1) == 'g'
        set(get(h,'xlabel'),'string','\Delta\rho g/cm^3','fontsize',10,'fontweight','bold');   % title
    end
end



for it = 1:length_comp
    temp=cutspace(Inversion_component(it,:));
    eval([temp '=' 'dobs((it-1)*n_data+1:it*n_data);']);
    eval([temp '_pre' '=' 'dpre((it-1)*n_data+1:it*n_data);']);
end



minx=min(xob);maxx=max(xob);
miny=min(yob);maxy=max(yob);
xs=(maxx-minx)/contint;
ys=(maxy-miny)/contint;
xg=minx:xs:maxx;
yg=miny:ys:maxy;
[xgg,ygg]=meshgrid(xg,yg);



for it=1:length_comp
    temp = cutspace(Inversion_component(it,:));
    eval([temp '_inter' '=' 'griddata(xob,yob,' temp ',xgg,ygg)' ';']);
    eval([temp '_pre_inter' '=' 'griddata(xob,yob,' temp '_pre' ',xgg,ygg)' ';']);
    clear temp;
end

for it =1:length_comp
    temp=cutspace(Inversion_component(it,:));
    if temp(1) == 'h'
        set(get(h,'xlabel'),'string','\Delta\chi','fontsize',10,'fontweight','bold');   % title
    elseif temp(1) == 'g'
        set(get(h,'xlabel'),'string','\Delta\rho g/cm^3','fontsize',10,'fontweight','bold');   % title
    end
end

for it=1:length_comp
    temp=cutspace(Inversion_component(it,:));
    figure(it+2);
    subplot(2,1,1);
    set(gca,'fontweight','bold');
    eval(['contourf(xgg,ygg,' temp,'_inter,15)' ';']);
    colorbar; h=colorbar;
    title(['Observed ' temp],'fontweight','bold');
    set(h,'fontweight','bold');
    
    if temp(1) == 'h'
        if length(temp) == 2;
        set(get(h,'xlabel'),'string','nT','fontweight','bold');
        elseif length(temp) ==3
            set(get(h,'xlabel'),'string','nT/m','fontweight','bold');
        end
    elseif temp(1) == 'g'
        if length(temp) == 2
        set(get(h,'xlabel'),'string','mGal','fontweight','bold');
        elseif length(temp) ==3
            set(get(h,'xlabel'),'string','mGal','fontweight','bold');
        end
    end
    set(gca,'fontweight','bold');
    axis equal;
    xlabel('x.m','fontweight','bold');
    ylabel('y.m','fontweight','bold');
    axis equal;
    axis([minx maxx miny maxy])
    
    subplot(2,1,2);
    set(gca,'fontweight','bold');
    eval(['contourf(xgg,ygg,' temp,'_pre_inter,15)' ';']);
    colorbar; colormap(jet); h=colorbar();
    title(['Predicted ' temp],'fontweight','bold');
    set(h,'fontweight','bold');

    if temp(1) == 'h'
        if length(temp) == 2;
        set(get(h,'xlabel'),'string','nT','fontweight','bold');
        elseif length(temp) ==3
            set(get(h,'xlabel'),'string','nT/m','fontweight','bold');
        end
    elseif temp(1) == 'g'
        if length(temp) == 2
        set(get(h,'xlabel'),'string','mGal','fontweight','bold');
        elseif length(temp) ==3
            set(get(h,'xlabel'),'string','mGal','fontweight','bold');
        end
    end
    
    set(gca,'fontweight','bold');
    axis equal;
    xlabel('x.m','fontweight','bold');
    ylabel('y.m','fontweight','bold');
    axis equal;
    axis([minx maxx miny maxy]);
%     saveas(gcf,['inversion_result\' Inversion_component(it,:) '_fitting'],'jpg');
%     saveas(gcf,['inversion_result\' Inversion_component(it,:) '_fitting'],'fig');  
end


figure(4);
plot(2:n,misfit(2:n),'linewidth',2);
set(gca,'fontweight','bold');
xlabel('Iteration number');
ylabel('Normalized misfit');
set(gca,'fontweight','bold');

% % 绘制先验模型
% if isfield(showpar, 'PriorModel3d') && isfield(showpar, 'PriorModelVector')
%     disp('绘制先验模型...');
%     
%     % 提取先验模型数据
%     prior_m = showpar.PriorModelVector;
%     prior_M = showpar.PriorModel3d;
%     
%     % 显示维度信息
%     disp(['prior_M维度: ', num2str(size(prior_M))]);
%     disp(['期望维度: ', num2str([nz, nx, ny])]);
%     
%     % 检查先验模型维度
%     if size(prior_M,1) ~= nz || size(prior_M,2) ~= nx || size(prior_M,3) ~= ny
%         warning('先验模型维度 (%d×%d×%d) 与期望维度 (%d×%d×%d) 不匹配，跳过绘制', ...
%             size(prior_M,1), size(prior_M,2), size(prior_M,3), nz, nx, ny);
%         return;
%     end
%     
%     % 先验模型3D体视图 - 使用corob函数绘制方块图，参考figure2
%     figure(5);
%     title('3D View of Prior Model','fontweight','bold');
%     view(3);
%     hold on;
%     
%     % 计算先验模型的颜色
%     prior_color = Jetconvert(prior_m);
%     
%     % 绘制先验模型的3D方块图
%     for it=1:length(xm)
%         if prior_m(it)>cutoff*max(prior_m)
%             corob(xm(it),ym(it),zm(it),[dx,dy,dz],prior_color(it,:));
%         end
%     end
%     
%     grid on;
%     set(gca, 'fontweight','bold');
%     xlabel('x.m','fontweight','bold');
%     ylabel('y.m','fontweight','bold');
%     zlabel('z.m','fontweight','bold');
%     set(gca, 'fontweight','bold');
%     set(gca,'zdir','reverse');
%     
%     colormap(jet);
%     h=colorbar;
%     set(h,'fontweight','bold');
%     
%     % 设置颜色条范围和标签
%     if max(prior_m) > min(prior_m)
%         barLims=[min(prior_m) max(prior_m)];
%         dif=diff(barLims);
%         Nlbl=11;
%         lbls=[barLims(1):dif/(Nlbl-1):barLims(2)]';
%         CBlim=get(h,'Ylim');
%         Ytcks=CBlim(1):diff(CBlim)/(Nlbl-1):CBlim(2);
%         set(h,'Ytick',Ytcks,'Yticklabel',num2str(lbls,3));
%     end
%     
%     % 设置颜色条标签
%     for it =1:length_comp
%         temp=cutspace(Inversion_component(it,:));
%         if temp(1) == 'h'
%             set(get(h,'xlabel'),'string','\Delta\chi','fontsize',10,'fontweight','bold');
%         elseif temp(1) == 'g'
%             set(get(h,'xlabel'),'string','\Delta\rho g/cm^3','fontsize',10,'fontweight','bold');
%         end
%     end
%     
%     axis equal;
%     axis([min(xob) max(xob) min(yob) max(yob) min(topz) min(topz)+depth]);
% else
%     disp('先验模型数据不存在，跳过绘制');
% end

