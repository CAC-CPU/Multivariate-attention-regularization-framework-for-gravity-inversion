%mt is a 3D matrix denotes the model parameters;
%the size of mt is : [nz,nx,ny]
%cut_block is the number of cells need to be cutted for topography in each
%colomn, the size of cut_block is [ny, nx]
%reshape order: z faster > x >y


function mt=reshape_back_cut(m,nx,ny,nz,cut_block)



if [ny nx]~=size(cut_block)|(length(m)+sum(sum(cut_block)))~=nx*ny*nz
    disp('Warning! Size is not correct to reshape!')
end

ind=0;
for iny=1:ny
    for inx=1:nx
        for inz=1:nz
            cut_value=cut_block(iny,inx);
            if inz<=cut_value
                mt(inz,inx,iny)=0;
            else
                ind=ind+1;
                mt(inz,inx,iny)=m(ind);
            end
        end
    end
end
m=m';
