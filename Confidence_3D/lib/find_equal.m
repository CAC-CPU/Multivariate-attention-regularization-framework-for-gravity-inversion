function tmd=find_equal(m,n)

% Return the same elements from two vectors;

%      Example:
%      a=[1 2 3 4]; b = [3 4 5 6]
%      c = find_equal(a,b), returns [3 4]

ind=0;
for im=1:length(m)
    for in=1:length(n)
        if m(im)==n(in)
            ind=ind+1;
            tmd(ind)=m(im);
        end
    end
end