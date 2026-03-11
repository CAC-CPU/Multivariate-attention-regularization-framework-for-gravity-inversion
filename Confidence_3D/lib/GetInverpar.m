function [A d inverpar2] = GetInverpar(inverpar,data,topo)

% generate inverpar2 for inversion
Ho = inverpar.Ho;
gam = 6.67*10^(-11)*10^8;
I = inverpar.Inclination;
D = inverpar.Declination;
Azm = inverpar.Azimuth;
CN = inverpar.MaxIteration;
ao = 0;
qq = inverpar.DecayCoef;
e = inverpar.focus;
final_misfit = inverpar.tol;
rew = inverpar.reweight;
edge_x = inverpar.edge_x;
edge_y = inverpar.edge_y;
depth = inverpar.depth;
nx = inverpar.numx;
ny = inverpar.numy;
nz = inverpar.numz;
Inversion_component = inverpar.component;
mapr_temp = inverpar.Maprior;
mini_temp = inverpar.Minitial;


if size(unique(Inversion_component)) == 2
    disp('Can not do joint inversion for magnetic and gravity field. Please check inversion component');
    return;
end


% Read data:
xob=data(:,1); yob=data(:,2); zob=data(:,3);
lengthcomp = size(Inversion_component,1);

if lengthcomp+3 ~= size(data,2)
    disp('The input data does not have the same potential field components as specified by inversion component');
    return;
end

for it = 1:size(Inversion_component,1)
    temp = cutspace(Inversion_component(it,:));
    eval([temp '=data(:,it+3);']);
    clear temp;
end


n_data = length(xob);
xp_max=max(xob); xp_min=min(xob);
yp_max=max(yob); yp_min=min(yob);
L_x=(xp_max-xp_min)+2*edge_x; 
L_y=(yp_max-yp_min)+2*edge_y; 
dx=L_x/nx; dy=L_y/ny; dz=depth/nz;

lx=cos(I*pi/180)*sin((D+Azm)*pi/180);
ly=cos(I*pi/180)*cos((D+Azm)*pi/180);
lz=sin(I*pi/180);

topx=topo(:,1); topy=topo(:,2); topz=topo(:,3); 
x=xp_min-edge_x+0.5*dx:dx:xp_max+edge_x-0.5*dx;
y=yp_min-edge_y+0.5*dy:dy:yp_max+edge_y-0.5*dy;
z=0.5*dz+min(topz):dz:depth-0.5*dz+min(topz);

[temp1,cut_block]=top_surf_cutblock(x,y,dz,topx,topy,topz);
[XX,ZZ,YY]=meshgrid(x,z,y);
xt=reshape_cut(XX,cut_block);
yt=reshape_cut(YY,cut_block);
zt=reshape_cut(ZZ,cut_block);
num_cell=length(xt);


d = zeros(n_data*lengthcomp,1);
A = zeros(n_data*lengthcomp,num_cell);
for it=1:lengthcomp
    eval(['A' Inversion_component(it,:) '=zeros(n_data,num_cell);']);
    temp = cutspace(Inversion_component(it,:));
    d(1+(it-1)*n_data:it*n_data) = eval(temp);
end


if strcompare(Inversion_component,'gz')==1
    for kk=1:n_data
        for ll = 1:num_cell
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Agz(kk,ll)=(gam*(zt(ll)-zob(kk))/dist^3)*dx*dy*dz;
        end 
    end
end

if strcompare(Inversion_component,'gxx')==1
    for kk=1:n_data
        for ll = 1:num_cell
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);  
            Agxx(kk,ll)=gam*(3*(xt(ll)-xob(kk))^2/dist^2-1)*dx*dy*dz/dist^3;
        end 
    end
end

if strcompare(Inversion_component,'gxy')==1
    for kk=1:n_data
        for ll = 1:num_cell
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);  
            Agxy(kk,ll)=gam*(3*(xt(ll)-xob(kk))*(yt(ll)-yob(kk))/dist^2)*dx*dy*dz/dist^3;
        end 
    end
end

if strcompare(Inversion_component,'gxz')==1
    for kk=1:n_data
        for ll = 1:num_cell
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);  
            Agxz(kk,ll)=gam*(3*(xt(ll)-xob(kk))*(zt(ll)-zob(kk))/dist^2)*dx*dy*dz/dist^3;
        end 
    end
end
        
if strcompare(Inversion_component,'gyx')==1
    for kk=1:n_data
        for ll = 1:num_cell
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);  
            Agyx(kk,ll)=gam*(3*(xt(ll)-xob(kk))*(yt(ll)-yob(kk))/dist^2)*dx*dy*dz/dist^3;
        end 
    end
end    

if strcompare(Inversion_component,'gyy')==1
    for kk=1:n_data
        for ll = 1:num_cell
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);  
            Agyy(kk,ll)=gam*(3*(yt(ll)-yob(kk))^2/dist^2-1)*dx*dy*dz/dist^3;
        end 
    end
end  
                
if strcompare(Inversion_component,'gyz')==1
    for kk=1:n_data
        for ll = 1:num_cell
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);  
            Agyz(kk,ll)=gam*(3*(yt(ll)-yob(kk))*(zt(ll)-zob(kk))/dist^2)*dx*dy*dz/dist^3;
        end 
    end
end  

if strcompare(Inversion_component,'gzx')==1
    for kk=1:n_data
        for ll = 1:num_cell
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);  
            Agzx(kk,ll)=gam*(3*(xt(ll)-xob(kk))*(zt(ll)-zob(kk))/dist^2)*dx*dy*dz/dist^3;
        end 
    end
end

if strcompare(Inversion_component,'gzy')==1
    for kk=1:n_data
        for ll = 1:num_cell
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);  
            Agzy(kk,ll)=gam*(3*(yt(ll)-yob(kk))*(zt(ll)-zob(kk))/dist^2)*dx*dy*dz/dist^3;
        end 
    end
end  


if strcompare(Inversion_component,'gzz')==1
    for kk=1:n_data
        for ll = 1:num_cell
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);  
             Agzz(kk,ll)=gam*(3*(zt(ll)-zob(kk))^2/dist^2-1)*dx*dy*dz/dist^3;
        end 
    end
end 
                
       
    
if strcompare(Inversion_component,'hx')==1
    for kk=1:n_data
        for ll = 1:num_cell
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Ahx(kk,ll)=(-Ho/dist^3)*(lx-3*t*(xt(ll)-xob(kk))/dist^2)*dx*dy*dz;
        end 
    end
end

if strcompare(Inversion_component,'hy')==1
    for kk=1:n_data
        for ll = 1:num_cell
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Ahy(kk,ll)=(-Ho/dist^3)*(ly-3*t*(yt(ll)-yob(kk))/dist^2)*dx*dy*dz;
        end
    end
end

if strcompare(Inversion_component,'hz')==1
    for kk=1:n_data
        for ll = 1:num_cell
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Ahz(kk,ll)=(-Ho/dist^3)*(lz-3*t*(zt(ll)-zob(kk))/dist^2)*dx*dy*dz;
        end
    end
end


if strcompare(Inversion_component,'ht')==1
    for kk=1:n_data
        for ll = 1:num_cell
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Aht(kk,ll)=(Ho/dist^3)*(3*t^2/dist^2-1)*dx*dy*dz;
        end
    end
end

if strcompare(Inversion_component,'hxx')==1
    for kk=1:n_data
        for ll = 1:num_cell
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Ahxx(kk,ll)=3*Ho*(((-lx*(xt(ll)-xob(kk))-t)*dist^2+2*(xob(kk)-xt(ll))^2*t)...
                /dist^7+(lx-3*t*(xt(ll)-xob(kk))/dist^2)*(xob(kk)-xt(ll))/dist^5)*dx*dy*dz;
        end
    end
end


if strcompare(Inversion_component,'hxy')==1
    for kk=1:n_data
        for ll = 1:num_cell
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Ahxy(kk,ll)=3*Ho*(((-ly*(xt(ll)-xob(kk)))*dist^2-2*(yob(kk)-yt(ll))*(xt(ll)-xob(kk))*t)...
                /dist^7+(lx-3*t*(xt(ll)-xob(kk))/dist^2)*(yob(kk)-yt(ll))/dist^5)*dx*dy*dz;
        end
    end
end


if strcompare(Inversion_component,'hxz')==1
    for kk=1:n_data
        for ll = 1:num_cell
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Ahxz(kk,ll)=3*Ho*(((-lz*(xt(ll)-xob(kk)))*dist^2-2*(zob(kk)-zt(ll))*(xt(ll)-xob(kk))*t)...
               /dist^7+(lx-3*t*(xt(ll)-xob(kk))/dist^2)*(zob(kk)-zt(ll))/dist^5)*dx*dy*dz;
        end
    end
end

if strcompare(Inversion_component,'hyy')==1
    for kk=1:n_data
        for ll = 1:num_cell
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Ahyy(kk,ll)=3*Ho*(((-ly*(yt(ll)-yob(kk))-t)*dist^2+2*(yob(kk)-yt(ll))^2*t)/...
                dist^7+(ly-3*t*(yt(ll)-yob(kk))/dist^2)*(yob(kk)-yt(ll))/dist^5)*dx*dy*dz;
        end
    end
end

if strcompare(Inversion_component,'hyx')==1
    for kk=1:n_data
        for ll = 1:num_cell
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
             Ahyx(kk,ll)=3*Ho*(((-lx*(yt(ll)-yob(kk)))*dist^2-2*(xob(kk)-xt(ll))*(yt(ll)-yob(kk))*t)/...
                 dist^7+(ly-3*t*(yt(ll)-yob(kk))/dist^2)*(xob(kk)-xt(ll))/dist^5)*dx*dy*dz;
        end
    end
end

if strcompare(Inversion_component,'hyz')==1
    for kk=1:n_data
        for ll = 1:num_cell
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Ahyz(kk,ll)=3*Ho*(((-lz*(yt(ll)-yob(kk)))*dist^2-2*(zob(kk)-zt(ll))*(yt(ll)-yob(kk))*t)/...
                 dist^7+(ly-3*t*(yt(ll)-yob(kk))/dist^2)*(zob(kk)-zt(ll))/dist^5)*dx*dy*dz;
        end
    end
end


if strcompare(Inversion_component,'hzz')==1
    for kk=1:n_data
        for ll = 1:num_cell
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Ahzz(kk,ll)=3*Ho*(((-lz*(zt(ll)-zob(kk))-t)*dist^2+2*(zob(kk)-zt(ll))^2*t)/...
                 dist^7+(lz-3*t*(zt(ll)-zob(kk))/dist^2)*(zob(kk)-zt(ll))/dist^5)*dx*dy*dz;
        end
    end
end

if strcompare(Inversion_component,'hzx')==1
    for kk=1:n_data
        for ll = 1:num_cell
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Ahzx(kk,ll)=3*Ho*(((-lx*(zt(ll)-zob(kk)))*dist^2-2*(xob(kk)-xt(ll))*(zt(ll)-zob(kk))*t)/...
                dist^7+(lz-3*t*(zt(ll)-zob(kk))/dist^2)*(xob(kk)-xt(ll))/dist^5)*dx*dy*dz;
        end
    end
end

if strcompare(Inversion_component,'hzy')==1
    for kk=1:n_data
        for ll = 1:num_cell
            t=lx*(xt(ll)-xob(kk))+ly*(yt(ll)-yob(kk))+lz*(zt(ll)-zob(kk));
            dist=sqrt((xt(ll)-xob(kk))^2+(yt(ll)-yob(kk))^2+(zt(ll)-zob(kk))^2);
            Ahzy(kk,ll)=3*Ho*(((-ly*(zt(ll)-zob(kk)))*dist^2-2*(yob(kk)-yt(ll))*(zt(ll)-zob(kk))*t)/...
                dist^7+(lz-3*t*(zt(ll)-zob(kk))/dist^2)*(yob(kk)-yt(ll))/dist^5)*dx*dy*dz;
        end
    end
end


for it=1:lengthcomp
    temp = cutspace(Inversion_component(it,:));
    A(1+(it-1)*n_data:it*n_data,:) = eval(['A' temp ';']);
end

inverpar2.maxiteration=CN;
inverpar2.focus = e;
inverpar2.tol = final_misfit;
inverpar2.reweight = rew;
inverpar2.InitialReg = ao;
inverpar2.DecayCoef = qq;
inverpar2.mapr = reshape_cut(mapr_temp,cut_block);
inverpar2.mstart = reshape_cut(mini_temp,cut_block);
inverpar2.cut_block = cut_block;
inverpar2.x = x;
inverpar2.y = y;
inverpar2.z = z;
inverpar2.numx = nx;
inverpar2.numy = ny;
inverpar2.numz = nz;
inverpar2.xobserved = xob;
inverpar2.yobserved = yob;
inverpar2.zobserved = zob;
inverpar2.Xtopography = topx;
inverpar2.Ytopography = topy;
inverpar2.Ztopography = topz;
inverpar2.component = inverpar.component;
inverpar2.upbound = inverpar.upbound;
inverpar2.lowbound = inverpar.lowbound;
inverpar2.cutoff = inverpar.cutoff;
inverpar2.istraditioninv = inverpar.istraditioninv;
inverpar2.apha = inverpar.apha;