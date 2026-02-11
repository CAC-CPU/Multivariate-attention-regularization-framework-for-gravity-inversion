function forward = set_forward_params(forwardObs)

[forward,~] = config_params();

% 网格坐标计算
forward.gridsNum = forward.zGridsNum*forward.xGridsNum;
forward.dx = forward.width/forward.xGridsNum;
forward.dz = forward.depth/forward.zGridsNum;
forward.xGridsCenters = linspace(0 + 0.5*forward.dx,forward.width - 0.5*forward.dx,forward.xGridsNum)';
forward.zGridsCenters = linspace(0 + 0.5*forward.dz,forward.depth - 0.5*forward.dz,forward.zGridsNum)';
[forward.xGridsMesh,forward.zGridsMesh] = meshgrid(forward.xGridsCenters,forward.zGridsCenters);
forward.xGridsVector = forward.xGridsMesh(:);
forward.zGridsVector = forward.zGridsMesh(:);

% 观测坐标导入
forward.dataNum = size(forwardObs,1);
forward.xObs = forwardObs(:,1);
forward.zObs = forwardObs(:,2);
end