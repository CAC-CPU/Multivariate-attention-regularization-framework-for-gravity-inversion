% test_mapr_loading - 测试三维正则化模型的加载功能
% 
% 此脚本用于测试read_3d_mapr函数和修改后的makeinverpar.m文件
% 确保它们能够正确处理不同格式的正则化模型

% 测试1: 调用makeinverpar函数
disp('测试1: 调用makeinverpar函数...');
try
    inverpar = makeinverpar();
    disp('makeinverpar函数调用成功');
    disp(['正则化模型维度: ', num2str(size(inverpar.Maprior))]);
    disp(['期望维度: ', num2str([inverpar.numz, inverpar.numx, inverpar.numy])]);
    if all(size(inverpar.Maprior) == [inverpar.numz, inverpar.numx, inverpar.numy])
        disp('正则化模型维度正确');
    else
        disp('正则化模型维度不正确');
    end
catch ME
    disp(['makeinverpar函数调用失败: ', ME.message]);
end

% 测试2: 直接测试read_3d_mapr函数
disp('\n测试2: 直接测试read_3d_mapr函数...');
try
    addpath lib;
    numx = 20;
    numy = 10;
    numz = 20;
    Maprior = read_3d_mapr('mapr.xlsx', numx, numy, numz);
    disp('read_3d_mapr函数调用成功');
    disp(['正则化模型维度: ', num2str(size(Maprior))]);
    disp(['期望维度: ', num2str([numz, numx, numy])]);
    if all(size(Maprior) == [numz, numx, numy])
        disp('正则化模型维度正确');
    else
        disp('正则化模型维度不正确');
    end
catch ME
    disp(['read_3d_mapr函数调用失败: ', ME.message]);
end

% 测试3: 检查反演流程是否正常
disp('\n测试3: 检查反演流程是否正常...');
try
    % 加载数据和参数
    load data.dat;
    load topo.dat;
    inverpar = makeinverpar();
    
    % 检查GetInverpar函数是否能正常运行
    [A, d, inverpar2] = GetInverpar(inverpar, data, topo);
    disp('GetInverpar函数调用成功');
    disp(['正则化模型在GetInverpar后维度: ', num2str(size(inverpar2.mapr))]);
    disp(['单元格数量: ', num2str(length(inverpar2.mapr))]);
    
    % 检查linearinv函数是否能正常运行（只运行一次迭代）
    % 注意: 这里只测试函数调用，不执行完整反演
    if exist('linearinv', 'file') == 2
        disp('linearinv函数存在');
    else
        disp('linearinv函数不存在');
    end
    
    if exist('loginv', 'file') == 2
        disp('loginv函数存在');
    else
        disp('loginv函数不存在');
    end
    
catch ME
    disp(['反演流程测试失败: ', ME.message]);
end

% 测试4: 显示正则化模型的统计信息
disp('\n测试4: 显示正则化模型的统计信息...');
try
    inverpar = makeinverpar();
    disp(['正则化模型最小值: ', num2str(min(inverpar.Maprior(:)))]);
    disp(['正则化模型最大值: ', num2str(max(inverpar.Maprior(:)))]);
    disp(['正则化模型平均值: ', num2str(mean(inverpar.Maprior(:)))]);
    disp(['正则化模型非零元素数量: ', num2str(sum(inverpar.Maprior(:) ~= 0))]);
catch ME
    disp(['统计信息测试失败: ', ME.message]);
end

disp('\n测试完成!');
