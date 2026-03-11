%mt is a 3D matrix denotes the model parameters;
%the size of mt is : [nz,nx,ny]
%cut_block is the number of cells need to be cutted for topography in each
%colomn, the size of cut_block is [ny, nx]
%reshape order: z faster > x >y


function m=reshape_cut(mt,cut_block)

[nz,nx,ny]=size(mt);

if [ny nx]~=size(cut_block)
    disp('Warning! Size is not correct to reshape!')
end

ind=0;
for iny=1:ny
    for inx=1:nx
        for inz=1:nz
            cut_value=cut_block(iny,inx);
            if inz>cut_value
                ind=ind+1;
                m(ind)=mt(inz,inx,iny);
%                 m(ind)=mt((iny-1)*(nx*nz)+(inx-1)*nz+inz);
            end
        end
    end
end
m=m';
