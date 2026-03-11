%% === 读取 mat 文件 ===
matfile = 'showpar.mat';   % 你的文件名
data = load(matfile);              % 载入 struct

showpar = data.showpar;

%% === 提取 Model3d 与 ModelVector ===
M3d = showpar.Model3d;      % 10×20×10
Mvec = showpar.ModelVector; % 2000×1

%% === 将 3D 矩阵展开为 1D（2000×1）===
M3d_1d = M3d(:);    % 按列展开，2000×1

%% === 保存为 txt 文件 ===
% 文件名你可以改
% save('Model3d_1D.txt',  'M3d_1d', '-ascii');
save('inversion_result.txt', 'Mvec',   '-ascii');

disp('已成功将两个变量保存为 txt 文件！');
