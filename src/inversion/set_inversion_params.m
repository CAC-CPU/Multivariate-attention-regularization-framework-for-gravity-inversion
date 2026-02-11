function inversion = set_inversion_params(dataObs)

[~,inversion] = config_params();

zGridsNumPlot = inversion.zGridsNum;
xGridsNumPlot = inversion.xGridsNum;
% widthPlot = inversion.width;
widthPlot = max(dataObs(:,1)) - min(dataObs(:,1));
depth = inversion.depth;
bufferLength = inversion.bufferLength;
maxIteration = inversion.maxIteration;
mPriorFileName = inversion.mPriorFileName;

% 网格坐标计算
dx = widthPlot/xGridsNumPlot;
dz = depth/zGridsNumPlot;

xGridsNum = (2*bufferLength + widthPlot)/dx;
zGridsNum = zGridsNumPlot;
gridsNum = zGridsNum*xGridsNum;
xGridsEdgeMinimun = 0 + 0.5*dx - bufferLength;
xGridsEdgeMaximun = widthPlot - 0.5*dx + bufferLength;
xGridsCenters = xGridsEdgeMinimun:dx:xGridsEdgeMaximun;

zGridsEdgeMinimun = 0 + 0.5*dz;
zGridsEdgeMaximun = depth - 0.5*dz;
zGridsCenters = zGridsEdgeMinimun:dz:zGridsEdgeMaximun;

[xGridsMesh,zGridsMesh] = meshgrid(xGridsCenters,zGridsCenters);
xGridsVector = xGridsMesh(:);
zGridsVector = zGridsMesh(:);

% 绘图
gridsNumPlot = zGridsNumPlot*xGridsNumPlot;
xGridsEdgeMinimunPlot = 0 + 0.5*dx;
xGridsEdgeMaximunPlot = widthPlot - 0.5*dx;
xGridsCentersPlot = xGridsEdgeMinimunPlot:dx:xGridsEdgeMaximunPlot;

xIndexPlotMinimum = (xGridsNum - xGridsNumPlot)/2 + 1;
xIndexPlotMaximum = xIndexPlotMinimum + xGridsNumPlot - 1;
xIndexPlot = xIndexPlotMinimum:xIndexPlotMaximum;

% 储存
inversion.dx = dx;
inversion.dz = dz;

inversion.xGridsNum = xGridsNum;
inversion.zGridsNum = zGridsNum;
inversion.gridsNum = gridsNum;
inversion.xGridsEdgeMinimun = xGridsEdgeMinimun;
inversion.xGridsEdgeMaximun = xGridsEdgeMaximun;
inversion.xGridsCenters = xGridsCenters;
inversion.zGridsEdgeMinimun = zGridsEdgeMinimun;
inversion.zGridsEdgeMaximun = zGridsEdgeMaximun;
inversion.zGridsCenters = zGridsCenters;
inversion.xGridsMesh = xGridsMesh;
inversion.zGridsMesh = zGridsMesh;
inversion.xGridsVector = xGridsVector;
inversion.zGridsVector = zGridsVector;

inversion.xGridsEdgeMinimunPlot= xGridsEdgeMinimunPlot;
inversion.xGridsEdgeMaximunPlot = xGridsEdgeMaximunPlot;
inversion.xGridsCentersPlot = xGridsCentersPlot;
inversion.gridsNumPlot = gridsNumPlot;
inversion.xGridsNumPlot = xGridsNumPlot;
inversion.zGridsNumPlot = zGridsNumPlot;
inversion.xIndexPlot = xIndexPlot;
inversion.xIndexPlotMinimum = xIndexPlotMinimum;
inversion.xIndexPlotMaximum = xIndexPlotMaximum;

% 观测数据导入
inversion.dataObs = dataObs(:,3);
inversion.dataNum = size(dataObs,1);
inversion.xObs = dataObs(:,1);
inversion.zObs = dataObs(:,2);

% 初始化模型值，先验模型
inversion.mInit = zeros(inversion.zGridsNum,inversion.xGridsNum);
inversion.mInitVector = inversion.mInit(:);
inversion.mPrior = zeros(inversion.zGridsNum,inversion.xGridsNum);
inversion.mPrior = make_model_params(inversion);
inversion.mPrior = readmatrix(mPriorFileName);
inversion.mPriorVector = inversion.mPrior(:);

inversion.misfit = zeros(maxIteration,1);
inversion.alpha = zeros(maxIteration,1);
inversion.gamma = zeros(maxIteration,1);
inversion.theta = ones(maxIteration,1);
inversion.confidenceFactor = zeros(maxIteration,gridsNum);

end