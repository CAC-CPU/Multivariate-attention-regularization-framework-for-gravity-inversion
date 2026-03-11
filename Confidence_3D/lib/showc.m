function showc(showpar)

c = showpar.c;
C = zeros(size(c,1),1);
for n=1:size(c,1)
    C(:,1)=norm(c(:,1));
end

figure
plot(C)
title('normCn')
