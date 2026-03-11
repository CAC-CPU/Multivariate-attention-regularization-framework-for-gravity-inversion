function [Xm Ym Zm M]=make_model()

% The discretization of z is from the top of the topography


% xlength: the length of domain in x dimension
% ylength: the length of domain in y dimension
% zlength: the length of domain in z dimension
% nx: number of cells in x dimension
% ny: number of cells in y dimension
% nz: number of cells in z dimension

% mture: 3-D matrix for model parameter. 
% the first dimension: model in z direction
% the second dimension: model in x direction



xlength=2000; 
ylength=1000; 
depth=1000; 
nx=2;
ny=3;
nz=4;



%==========================================================================
% Don't change these paragraphy
dx=xlength/nx;
dy=ylength/ny;
dz=depth/nz;
x=0.5*dx:dx:xlength-0.5*dx;
y=0.5*dy:dy:ylength-0.5*dy;
load topo.dat;
topz = topo(:,3);
z=0.5*dz+min(topz):dz:depth-0.5*dz+min(topz);
mture=zeros(nz,nx,ny);
%==========================================================================

% mture(9,12,3:8)=1;mture(10,11,3:8)=1;mture(11,10,3:8)=1;
% mture(12,9,3:8)=1;mture(13,8,3:8)=1;mture(14,7,3:8)=1;
mture(2,1,1)=1; mture(1,2,2)=1; mture(2,2,3)=1;


M=mture(:);
[Xmesh,Ymesh,Zmesh]=meshgrid(x,y,z);
Xm=Xmesh(:);
Ym=Ymesh(:);
Zm=Zmesh(:);