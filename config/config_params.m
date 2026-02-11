function [forward,inversion] = config_params()
% 网格参数
forward.xGridsNum = 40;
forward.zGridsNum = 20;
forward.depth = 1000;
forward.width = 2000;
forward.makeModelFlag = 1;
forward.fileName = 'm_true.txt';

inversion.xGridsNum = 40;
inversion.zGridsNum = 20;
inversion.depth = 1000;
inversion.width = 2000;
inversion.bufferLength = 0;
inversion.makeModelFlag = 0;
inversion.fileName = 'm_apr.txt';
inversion.mPriorFileName = 'm_apr.txt';
% 迭代参数
inversion.kernelMode = 5;
inversion.maxIteration = 1000;
inversion.fmisfit = 0.001;
inversion.decayCoefficient = 0.85;
inversion.focusFactor = 0.004;
inversion.reweightIndex = [1:20];
inversion.reweightGap = 5;
inversion.lambda = 1;
% inversion.eta = 0.2;
inversion.plotStablizersIndex = 51;
end