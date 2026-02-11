function boundary = draw_model_boundary(params,mTrue)

dx = params.dx;
dz = params.dz;

[zBoundary,xBoundary] = find(mTrue ~= 0);
zBoundaryMin = min(zBoundary) - 1;
zBoundaryMax = max(zBoundary) + 1;
xBoundaryMin = min(xBoundary) - 1;
xBoundaryMax = max(xBoundary) + 1;
mBoundary = mTrue(zBoundaryMin:zBoundaryMax,xBoundaryMin:xBoundaryMax);
zNumBoundary = size(mBoundary,1) - 2;
xNumBoundary = size(mBoundary,2) - 2;
numBoundary = zNumBoundary*xNumBoundary;
upperBoundary = zeros(numBoundary,3);
lowerBoundary = zeros(numBoundary,3);
leftBoundary = zeros(numBoundary,3);
rightBoundary = zeros(numBoundary,3);
rows = 1;
for i = (1:zNumBoundary) + 1
    for j = (1:xNumBoundary) + 1
        if mBoundary(i,j) ~= 0 && mBoundary(i-1,j) == 0
            upperBoundary(rows,:) = [j-1+xBoundaryMin,i-1+zBoundaryMin,1];
        end
        if mBoundary(i,j) ~= 0 && mBoundary(i+1,j) == 0
            lowerBoundary(rows,:) = [j-1+xBoundaryMin,i-1+zBoundaryMin,1];
        end
        if mBoundary(i,j) ~= 0 && mBoundary(i,j-1) == 0
            leftBoundary(rows,:) = [j-1+xBoundaryMin,i-1+zBoundaryMin,1];
        end
        if mBoundary(i,j) ~= 0 && mBoundary(i,j+1) == 0
            rightBoundary(rows,:) = [j-1+xBoundaryMin,i-1+zBoundaryMin,1];
        end
        rows = rows + 1;
    end
end

deleteRows = upperBoundary(:,3) == 0;
upperBoundary(deleteRows,:) = [];
upperBoundary(:,3) = [];
deleteRows = lowerBoundary(:,3) == 0;
lowerBoundary(deleteRows,:) = [];
lowerBoundary(:,3) = [];
deleteRows = leftBoundary(:,3) == 0;
leftBoundary(deleteRows,:) = [];
leftBoundary(:,3) = [];
deleteRows = rightBoundary(:,3) == 0;
rightBoundary(deleteRows,:) = [];
rightBoundary(:,3) = [];

xBoundaryEdge = [upperBoundary(:,1)+[-1 0];...
                 lowerBoundary(:,1)+[-1 0];...
                 leftBoundary(:,1)+[-1 -1];...
                 rightBoundary(:,1)+[0 0]]*dx;
zBoundaryEdge = [upperBoundary(:,2)+[-1 -1];...
                 lowerBoundary(:,2)+[0 0];...
                 leftBoundary(:,2)+[-1 0];...
                 rightBoundary(:,2)+[-1 0]].*dz;

boundary = [xBoundaryEdge zBoundaryEdge];

end