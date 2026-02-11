function A = generate_forward_operator(params)

dataNum = params.dataNum;
GridsNum = params.gridsNum;
xGridsVector = params.xGridsVector;
xObs = params.xObs;
zGridsVector = params.zGridsVector;
zObs = params.zObs;
dx = params.dx;
dz = params.dz;

G = 6.67*10^(-3);
A = zeros(dataNum,GridsNum);
for i = 1:dataNum
    for j = 1:GridsNum
        dist   = ...
            sqrt((xGridsVector(j)-xObs(i))^2+(zGridsVector(j)-zObs(i))^2);
        A(i,j) = ...
            (G*(zGridsVector(j) - zObs(i))/dist^3)*dx*dz;
    end
end
end