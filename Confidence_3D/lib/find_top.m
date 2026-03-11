%find z value for each observation point;

function zob=find_top(xob,yob,topx,topy,topz)
for it=1:length(xob)
    diftopx=abs(topx-xob(it)); diftopy=abs(topy-yob(it));
    indx=find(diftopx==min(diftopx));
    indy=find(diftopy==min(diftopy));
    ind=find_equal(indx,indy);
    if length(ind)~=1
        ind=ind(1);
    end
    zob(it)=topz(ind);
end

