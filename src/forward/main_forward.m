function main_forward()
% 观测坐标导入
generate_observation_coords;
obsData = load('forward_obs.txt');
forwardParams = set_forward_params(obsData);
A = generate_forward_operator(forwardParams);
mTrue = make_model_params(forwardParams);
dataObs = A*mTrue(:);
data = [obsData,dataObs];
boundary = draw_model_boundary(forwardParams,mTrue);
plot_forward_params(forwardParams,dataObs,mTrue,boundary);
save('data/data.txt','data','-ascii','-double','-tabs');
disp("正演程序运行结束，模拟观测数据结果已经保存在\data\data.txt");
end