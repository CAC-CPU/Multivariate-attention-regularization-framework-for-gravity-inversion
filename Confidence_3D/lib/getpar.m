% Read forward parameter from forwardpar.mat

% Ho: geomagnetic inducing field strength in the unit of nT
% gam: universal gravitational constant
% I: inclination of geomagnetic field
% D: Declination of geomagnetic field
% Azm: azimuth
% lx, ly and lz: direction consine of geomagnetic field in x,y,z direciotn
% noise: noise level of the synthetic data
% xp: x coordinate of the obsered data;
% yp: y coordinate of the oberved data;
% x: x coordinate of the model cells
% y: y coordinate of the model cells
% z: z coordinate of the model cells
% mt: 3-D synthetic model matrix
% dx: size of each cell in x dimension
% dy: size of each cell in y dimension
% dz: size of each cell in z dimension
% topx: x coordinate of topography data;
% topy: y coordinate of topography data;
% topz: z coordinate of topography data;
% length_comp: the number of potential field component to be modeled

Ho = forwardpar.Ho;
gam = 6.67*10^(-11)*10^8; %gravity constant
I = forwardpar.Inclination;
D = forwardpar.Declination;
Azm = forwardpar.Azimuth;
lx=cos(I*pi/180)*sin((D+Azm)*pi/180);
ly=cos(I*pi/180)*cos((D+Azm)*pi/180);
lz=sin(I*pi/180);
noise=forwardpar.noise;

[xtemp ytemp ztemp Mtemp] = make_model();
x = unique(xtemp);
y = unique(ytemp);
z = unique(ztemp);
mt=reshape(Mtemp,length(z),length(x),length(y));
dx = x(2) - x(1);
dy = y(2) - y(1);
dz = z(2) - z(1);
load topo.dat; % load topography data;
topx = topo(:,1);
topy = topo(:,2);
topz = topo(:,3);
z=z+min(topz);


load forwardxyz_obs.dat;
xob = forwardxyz_obs(:,1);
yob = forwardxyz_obs(:,2);
obslevel = forwardpar.obslevel;
switch obslevel
    case 0
        % The observation surface is according to topography
        zob=find_top(xob,yob,topx,topy,topz);
        if size(forwardxyz_obs,2) ~= 2
            warning('The parameter is set such that the z coordinate of the observation surface is the same as topography. The Third column of forwardxyz_obs is unnessceary');
        end
        
    case 1
        if size(forwardxyz_obs,2) ~= 3
            warning('forwardxyz_obs should has three columns in this case');
        end
        zob = forwardxyz_obs(:,3);
        
end
forward_comp=forwardpar.component;
length_comp=size(forward_comp,1);

