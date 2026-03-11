function forward()

% Main part of the forward modeling process.
% This function will load file 'forwardpar' and read forward parameter from
% it to produce potential field data and show the corresponding forward
% modeling result

forwardpar = makeforwardpar();
addpath lib;
getpar;
[tpm,cut_block]=top_surf_cutblock(x,y,dz,topx,topy,topz);
MT=reshape_cut(mt,cut_block);
[XX,ZZ,YY]=meshgrid(x,z,y);
xt=reshape_cut(XX,cut_block);
yt=reshape_cut(YY,cut_block);
zt=reshape_cut(ZZ,cut_block);
num_cell_top=length(xt);
sizex = length(unique(xob));
sizey = length(unique(yob));
xp    = unique(xob);
yp    = unique(yob);
n_data=sizex*sizey;


for it=1:length_comp
    eval(['A' forward_comp(it,:) '=zeros(n_data,num_cell_top);']);
end

% Compute the forward matrix for different potential field problem
if strcompare(forward_comp,'gz')==1
    for kk=1:n_data
        for ll = 1:num_cell_top
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Agz(kk,ll)=(gam*(zt(ll)-zob(kk))/dist^3)*dx*dy*dz;
        end 
    end
    gz = Agz*MT;
    gz(:)=gz(:).*(1+2*noise*(rand(n_data,1)-0.5)/100);
end

if strcompare(forward_comp,'gxx')==1
    for kk=1:n_data
        for ll = 1:num_cell_top
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);  
            Agxx(kk,ll)=gam*(3*(xt(ll)-xob(kk))^2/dist^2-1)*dx*dy*dz/dist^3;
        end 
    end
    gxx = Agxx*MT;
    gxx(:)=gxx(:).*(1+2*noise*(rand(n_data,1)-0.5)/100);
end

if strcompare(forward_comp,'gxy')==1
    for kk=1:n_data
        for ll = 1:num_cell_top
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);  
            Agxy(kk,ll)=gam*(3*(xt(ll)-xob(kk))*(yt(ll)-yob(kk))/dist^2)*dx*dy*dz/dist^3;
        end 
    end
    gxy = Agxy*MT;
    gxy(:)=gxy(:).*(1+2*noise*(rand(n_data,1)-0.5)/100);
end

if strcompare(forward_comp,'gxz')==1
    for kk=1:n_data
        for ll = 1:num_cell_top
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);  
            Agxz(kk,ll)=gam*(3*(xt(ll)-xob(kk))*(zt(ll)-zob(kk))/dist^2)*dx*dy*dz/dist^3;
        end 
    end
    gxz = Agxz*MT;
    gxz(:)=gxz(:).*(1+2*noise*(rand(n_data,1)-0.5)/100);
end
        
if strcompare(forward_comp,'gyx')==1
    for kk=1:n_data
        for ll = 1:num_cell_top
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);  
            Agyx(kk,ll)=gam*(3*(xt(ll)-xob(kk))*(yt(ll)-yob(kk))/dist^2)*dx*dy*dz/dist^3;
        end 
    end
    gyx = Agyx*MT;
    gyx(:)=gyx(:).*(1+2*noise*(rand(n_data,1)-0.5)/100);
end    

if strcompare(forward_comp,'gyy')==1
    for kk=1:n_data
        for ll = 1:num_cell_top
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);  
            Agyy(kk,ll)=gam*(3*(yt(ll)-yob(kk))^2/dist^2-1)*dx*dy*dz/dist^3;
        end 
    end
    gyy = Agyy*MT;
    gyy(:) = gyy(:).*(1+2*noise*(rand(n_data,1)-0.5)/100);
end  
                
if strcompare(forward_comp,'gyz')==1
    for kk=1:n_data
        for ll = 1:num_cell_top
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);  
            Agyz(kk,ll)=gam*(3*(yt(ll)-yob(kk))*(zt(ll)-zob(kk))/dist^2)*dx*dy*dz/dist^3;
        end 
    end
    gyz = Agyz*MT;
    gyz(:)=gyz(:).*(1+2*noise*(rand(n_data,1)-0.5)/100);
end  

if strcompare(forward_comp,'gzx')==1
    for kk=1:n_data
        for ll = 1:num_cell_top
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);  
            Agzx(kk,ll)=gam*(3*(xt(ll)-xob(kk))*(zt(ll)-zob(kk))/dist^2)*dx*dy*dz/dist^3;
        end 
    end
    gzx = Agzx*MT;
    gzx(:)=gzx(:).*(1+2*noise*(rand(n_data,1)-0.5)/100);
end

if strcompare(forward_comp,'gzy')==1
    for kk=1:n_data
        for ll = 1:num_cell_top
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);  
            Agzy(kk,ll)=gam*(3*(yt(ll)-yob(kk))*(zt(ll)-zob(kk))/dist^2)*dx*dy*dz/dist^3;
        end 
    end
    gzy = Agzy*MT;
    gzy(:)=gzy(:).*(1+2*noise*(rand(n_data,1)-0.5)/100);
end  


if strcompare(forward_comp,'gzz')==1
    for kk=1:n_data
        for ll = 1:num_cell_top
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);  
             Agzz(kk,ll)=gam*(3*(zt(ll)-zob(kk))^2/dist^2-1)*dx*dy*dz/dist^3;
        end 
    end
    gzz = Agzz*MT;
    gzz(:)=gzz(:).*(1+2*noise*(rand(n_data,1)-0.5)/100);
end 
                
       
    
if strcompare(forward_comp,'hx')==1
    for kk=1:n_data
        for ll = 1:num_cell_top
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Ahx(kk,ll)=(-Ho/dist^3)*(lx-3*t*(xt(ll)-xob(kk))/dist^2)*dx*dy*dz;
        end 
    end
    hx = Ahx*MT;
    hx(:)=hx(:).*(1+2*noise*(rand(n_data,1)-0.5)/100);
end

if strcompare(forward_comp,'hy')==1
    for kk=1:n_data
        for ll = 1:num_cell_top
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Ahy(kk,ll)=(-Ho/dist^3)*(ly-3*t*(yt(ll)-yob(kk))/dist^2)*dx*dy*dz;
        end
    end
    hy = Ahy*MT;
    hy(:)=hy(:).*(1+2*noise*(rand(n_data,1)-0.5)/100);
end

if strcompare(forward_comp,'hz')==1
    for kk=1:n_data
        for ll = 1:num_cell_top
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Ahz(kk,ll)=(-Ho/dist^3)*(lz-3*t*(zt(ll)-zob(kk))/dist^2)*dx*dy*dz;
        end
    end
    hz = Ahz*MT;
    hz(:)=hz(:).*(1+2*noise*(rand(n_data,1)-0.5)/100);
end


if strcompare(forward_comp,'ht')==1
    for kk=1:n_data
        for ll = 1:num_cell_top
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Aht(kk,ll)=(Ho/dist^3)*(3*t^2/dist^2-1)*dx*dy*dz;
        end
    end
    ht = Aht*MT;
    ht(:)=ht(:).*(1+2*noise*(rand(n_data,1)-0.5)/100);
end

if strcompare(forward_comp,'hxx')==1
    for kk=1:n_data
        for ll = 1:num_cell_top
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Ahxx(kk,ll)=3*Ho*(((-lx*(xt(ll)-xob(kk))-t)*dist^2+2*(xob(kk)-xt(ll))^2*t)...
                /dist^7+(lx-3*t*(xt(ll)-xob(kk))/dist^2)*(xob(kk)-xt(ll))/dist^5)*dx*dy*dz;
        end
    end
    hxx = Ahxx*MT;
    hxx(:)=hxx(:).*(1+2*noise*(rand(n_data,1)-0.5)/100);
end


if strcompare(forward_comp,'hxy')==1
    for kk=1:n_data
        for ll = 1:num_cell_top
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Ahxy(kk,ll)=3*Ho*(((-ly*(xt(ll)-xob(kk)))*dist^2-2*(yob(kk)-yt(ll))*(xt(ll)-xob(kk))*t)...
                /dist^7+(lx-3*t*(xt(ll)-xob(kk))/dist^2)*(yob(kk)-yt(ll))/dist^5)*dx*dy*dz;
        end
    end
    hxy = Ahxy*MT;
    hxy(:)=hxy(:).*(1+2*noise*(rand(n_data,1)-0.5)/100);
end


if strcompare(forward_comp,'hxz')==1
    for kk=1:n_data
        for ll = 1:num_cell_top
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Ahxz(kk,ll)=3*Ho*(((-lz*(xt(ll)-xob(kk)))*dist^2-2*(zob(kk)-zt(ll))*(xt(ll)-xob(kk))*t)...
               /dist^7+(lx-3*t*(xt(ll)-xob(kk))/dist^2)*(zob(kk)-zt(ll))/dist^5)*dx*dy*dz;
        end
    end
    hxz = Ahxz*MT;
    hxz(:)=hxz(:).*(1+2*noise*(rand(n_data,1)-0.5)/100);
end

if strcompare(forward_comp,'hyy')==1
    for kk=1:n_data
        for ll = 1:num_cell_top
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Ahyy(kk,ll)=3*Ho*(((-ly*(yt(ll)-yob(kk))-t)*dist^2+2*(yob(kk)-yt(ll))^2*t)/...
                dist^7+(ly-3*t*(yt(ll)-yob(kk))/dist^2)*(yob(kk)-yt(ll))/dist^5)*dx*dy*dz;
        end
    end
    hyy = Ahyy*MT;
    hyy(:)=hyy(:).*(1+2*noise*(rand(n_data,1)-0.5)/100);
end

if strcompare(forward_comp,'hyx')==1
    for kk=1:n_data
        for ll = 1:num_cell_top
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
             Ahyx(kk,ll)=3*Ho*(((-lx*(yt(ll)-yob(kk)))*dist^2-2*(xob(kk)-xt(ll))*(yt(ll)-yob(kk))*t)/...
                 dist^7+(ly-3*t*(yt(ll)-yob(kk))/dist^2)*(xob(kk)-xt(ll))/dist^5)*dx*dy*dz;
        end
    end
    hyx = Ahyx*MT;
    hyx(:)=hyx(:).*(1+2*noise*(rand(n_data,1)-0.5)/100);
end

if strcompare(forward_comp,'hyz')==1
    for kk=1:n_data
        for ll = 1:num_cell_top
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Ahyz(kk,ll)=3*Ho*(((-lz*(yt(ll)-yob(kk)))*dist^2-2*(zob(kk)-zt(ll))*(yt(ll)-yob(kk))*t)/...
                 dist^7+(ly-3*t*(yt(ll)-yob(kk))/dist^2)*(zob(kk)-zt(ll))/dist^5)*dx*dy*dz;
        end
    end
    hyz = Ahyz*MT;
    hyz(:)=hyz(:).*(1+2*noise*(rand(n_data,1)-0.5)/100);
end


if strcompare(forward_comp,'hzz')==1
    for kk=1:n_data
        for ll = 1:num_cell_top
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Ahzz(kk,ll)=3*Ho*(((-lz*(zt(ll)-zob(kk))-t)*dist^2+2*(zob(kk)-zt(ll))^2*t)/...
                 dist^7+(lz-3*t*(zt(ll)-zob(kk))/dist^2)*(zob(kk)-zt(ll))/dist^5)*dx*dy*dz;
        end
    end
    hzz = Ahzz*MT;
    hzz(:)=hzz(:).*(1+2*noise*(rand(n_data,1)-0.5)/100);
end

if strcompare(forward_comp,'hzx')==1
    for kk=1:n_data
        for ll = 1:num_cell_top
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Ahzx(kk,ll)=3*Ho*(((-lx*(zt(ll)-zob(kk)))*dist^2-2*(xob(kk)-xt(ll))*(zt(ll)-zob(kk))*t)/...
                dist^7+(lz-3*t*(zt(ll)-zob(kk))/dist^2)*(xob(kk)-xt(ll))/dist^5)*dx*dy*dz;
        end
    end
    hzx = Ahzx*MT;
    hzx(:)=hzx(:).*(1+2*noise*(rand(n_data,1)-0.5)/100);
end

if strcompare(forward_comp,'hzy')==1
    for kk=1:n_data
        for ll = 1:num_cell_top
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Ahzy(kk,ll)=3*Ho*(((-ly*(zt(ll)-zob(kk)))*dist^2-2*(yob(kk)-yt(ll))*(zt(ll)-zob(kk))*t)/...
                dist^7+(lz-3*t*(zt(ll)-zob(kk))/dist^2)*(yob(kk)-yt(ll))/dist^5)*dx*dy*dz;
        end
    end
    hzy = Ahzy*MT;
    hzy(:)=hzy(:).*(1+2*noise*(rand(n_data,1)-0.5)/100);
end

showmodelpar.x = x; showmodelpar.y = y; showmodelpar.z = z; 
showmodelpar.model = mt; showmodelpar.xd = xob; showmodelpar.yd =yob;
showmodelpar.zd =zob; showmodelpar.topox = topx; 
showmodelpar.topoy = topy; showmodelpar.topoz = topz;
showmodelpar.component = forward_comp;

save showmodelpar.mat showmodelpar;

if forwardpar.dispflag == 1 
    showmodel(showmodelpar);
end

data=zeros(length(xob),length_comp+3);
data(:,1) = xob; data(:,2) = yob; data(:,3) = zob;
for it =1:length_comp
    temp=cutspace(forward_comp(it,:));
    data(:,it+3)=eval(temp);
    tempdata=(reshape(eval(temp),sizex,sizey))';
    
    
    if forwardpar.dispflag == 1 
        figure(it+1);
        eval(['imagesc(xp,yp,tempdata)']);
        title(['Observed ' temp],'fontweight','bold');
        set(gca,'fontweight','bold');
        colorbar;
        xlabel('x','fontweight','bold');
        ylabel('y','fontweight','bold');
        axis equal;
        axis([xp(1) xp(end) yp(1) yp(end)]);
        colorbar;h=colorbar;
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
                set(get(h,'xlabel'),'string','mGal/m','fontweight','bold');
            end
        end
    
    end  
end

save data.dat data -ascii -double;