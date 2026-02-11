function smooth_model_params(inversion)
mapr = inversion.mPrior;
[xNonZero,zNonZero] = find(mapr ~= 0);
xNonZeroMiddle = ceil((max(xNonZero) - min(xNonZero))/2);
zNonZeroMiddle = ceil((max(zNonZero) - min(zNonZero))/2);
xGridsCenters = inversion.xGridsCenters;
zGridsCenters = inversion.zGridsCenters;
dx = inversion.dx;
dz = inversion.dz;
xCoordsCenters = abs(dx*(xGridsCenters - xNonZeroMiddle));
zCoordsCenters = abs(dz*(zGridsCenters - zNonZeroMiddle));

end