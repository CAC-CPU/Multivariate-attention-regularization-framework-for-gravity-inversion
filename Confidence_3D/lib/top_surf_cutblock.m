function [top_surf,cut_block]=top_surf_cutblock(x,y,dz,topx,topy,topz)

% This fucntion will extract topography data from topx topy and topz.
% In forwarding modeling and inversion process, we assume each cell
% indicated by x y and z has the same model parameter. The input argument
% dz is the size of the cell in z direction.
% top_surf is the discretized topography with the size of [length(y), length(x)]
% for given x and y
% cut_block is a matrix define how many cells need to be cutted from each
% column of the model according to topography
% Developed by Hongzhu Cai on March 6th, 2011.



% Last modified 02-17-2012



nx=length(x); ny=length(y);
top_surf=zeros(ny,nx);
cut_block=zeros(ny,nx);

for iny=1:ny
    for inx=1:nx       
        top_x_dif=abs(topx-x(inx));
        top_y_dif=abs(topy-y(iny));
        [indxm]=find(top_x_dif==min(top_x_dif));
        [indym]=find(top_y_dif==min(top_y_dif));
        equal_m=find_equal(indxm, indym);
        
        top_surf(iny,inx)=max(topz(equal_m));
        cut_block(iny,inx)=ceil((max(topz(equal_m))-min(topz))/dz);
        
    end
end