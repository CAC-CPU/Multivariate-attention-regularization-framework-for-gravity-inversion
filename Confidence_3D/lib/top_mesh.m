%This function will generate matrix TMDx, TMDy, TMDz which is the
%interpolated data for plot topography, mesh_contx,mesh_conty is the how
%many intervals of data in x direction and y direction. topx, topy topz are
%the coordinate for the topography.

function [TMDx,TMDy,TMDz]=top_mesh(mesh_contx,mesh_conty,topx,topy,topz)
% mesh_contx=100;
% mesh_conty=50;
xg=min(topx):(max(topx)-min(topx))/mesh_contx:max(topx);
yg=min(topy):(max(topy)-min(topy))/mesh_conty:max(topy);
[xgg,ygg]=meshgrid(xg,yg);
XG=reshape(xgg',length(xg)*length(yg),1);
XG=reshape(xgg',length(xg)*length(yg),1);
YG=reshape(ygg',length(xg)*length(yg),1);

INTER=griddata(topx,topy,topz,XG,YG);
TMDx=reshape(XG,length(xg),length(yg))';
TMDy=reshape(YG,length(xg),length(yg))';
TMDz=reshape(INTER,length(xg),length(yg))';