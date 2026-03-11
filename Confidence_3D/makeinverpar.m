function inverpar = makeinverpar()
%This function will generate a file which contains necessary parameter for
%the inversion process


% inverpar.Ho: inducing geomagnetic IN nT
% inverpar.Inclination: inclination of geomagnetic field
% inverpar.Declination: Declination of geomagnetic field
% inverpar.Azimuth: azimuth, Set to zero except:
% when model magnetic vector/tensor data that there issome angle between 
% inverpar.MaxIteration: maximum number of iteration
% inverpar.DecayCoef: a real number between 0 and 1. Recommend: 0.6 ~ 0.95
% inverpar.focus: focus parameter(very small). Default value 1e-3. If the
% convergence is not good, try other values from 1e-1 to 1e-5.
% inverpar.tol: tolerance level(normalized misfit)
% inverpar.reweight: vector store reweighting iteration point.
% Example, inverpar.reweight=[ 3 6 10]:
% In the first two iteration, do smooth inversion and follow in iteration
% 3,6, and 10, do image focusing. Usually do reweighting in every 2~5 and
% the first 1~2 (or 1~3, 1~4. Try different value) iteration only do smoothing inversion.
% inverpar.edge_x: how much the inversion domain in x is larger than the
% data covered area. Set to 0 if they are the same.
% inverpar.edge_y: as above
% inverpar.depth: depth of the inversion domain;
% inverpar.numx: number of cells in x direction;
% inverpar.numy: number of cells in y direction;
% inverpar.numz: number of cells in z direction;

% inverpar.component: potential field component used for inversion;
% for gravity field we have one component: 'gz'
% for gravity tensor field, there are nine components and five are
% independant:
% 'gxx', 'gzy', 'gxz', 'gyx', 'gyy', 'gyz', 'gzx', 'gzy', 'gzz'
% for magnetic total field, we have one component: 'ht'
% for magnetic vector field, we have three components:'hx', 'hy', 'hz'
% for magnetic tensor field, there are nine components and five are
% independant:
% 'hxx', 'hzy', 'hxz', 'hyx', 'hyy', 'hyz', 'hzx', 'hzy', 'hzz'

% inverpar.Modelspace:
% 1 for linear space 
% 2 for logrithmic space (good for set upper and lower bound for model)
% inverpar.upbound: upper boundary for the model 
% inverpar.lowbound: lower boundary for the model
% The above boundaries only works for logrithmic space inversion and
% symmetric boundary is appreciated.
% inverpar.Maprior: a 3-D matrix for a priori model and default is zero.
%   - 三维正则化模型支持从分层Excel文件中读取
%   - Excel文件结构: 包含多个工作表，每个工作表对应一个深度层
%   - 工作表命名为'z1', 'z2', ..., 'zn'，其中n为深度层数
%   - 每个工作表的尺寸为numy行×numx列，对应水平方向的网格
%   - 当只有一个工作表时，自动扩展为三维模型（所有深度层使用相同值）
% inverpar.Minitial: a 3-D matrix for initial model and default is zero.
% inverpar.cutoff: cut off value for 3D view. It's a real number between 0 and 1
% Example: cutoff=0.6 indicates the model value which is below 0.6 of the
% maximum value of the model will not be imaged for 3D view
% 1: display inversion result 0: not display
% 
% 三维正则化模型使用指南:
% 1. 创建或修改mapr.xlsx文件
% 2. 对于三维正则化模型:
%    a. 添加多个工作表，每个工作表对应一个深度层
%    b. 工作表命名为'z1', 'z2', ..., 'zn'
%    c. 每个工作表的尺寸为numy行×numx列
%    d. 在每个工作表中填入对应深度层的正则化参数
% 3. 对于二维正则化模型:
%    a. 只使用一个工作表
%    b. 工作表尺寸为numy行×numx列
%    c. 系统会自动将其扩展为三维模型（所有深度层使用相同值）
% 4. 正则化参数值越大，表示该区域的模型变化越受到限制
% 5. 全零正则化模型表示不施加任何先验约束

inverpar.Ho = 60000;
inverpar.Inclination = 90;
inverpar.Declination = 0;
inverpar.Azimuth = 0;
inverpar.MaxIteration = 50;
inverpar.DecayCoef = 0.85;
inverpar.focus = 1e-2;
inverpar.tol = 0.05;
% inverpar.reweight = [2 3 4];
inverpar.reweight = [5:8];
inverpar.edge_x = 0;
inverpar.edge_y = 0;
inverpar.depth = 1000;
inverpar.numx = 20;
inverpar.numy = 20;
inverpar.numz =20;
inverpar.component=char('gz');  % data类型
inverpar.Modelspace = 1;        % 1 线性 2 log
inverpar.upbound = 1;           % 上限约束，没用
inverpar.lowbound = -1;         % 下限约束，没用
inverpar.istraditioninv = 0;    % 是否进行传统正则化反演
inverpar.apha = 0.01;           % 传统正则化反演apha
inverpar.isdrawsclice = 0;      % 反演后是否继续切片
inverpar.Minitial = zeros(inverpar.numz,inverpar.numx,inverpar.numy);
inverpar.cutoff = 0.2;          % fig2显示密度大于cutoff的格子
inverpar.dispflag = 1;
addpath lib;
try
    inverpar.Maprior = read_3d_mapr('mapr_2.xlsx', inverpar.numx, inverpar.numy, inverpar.numz);
catch ME
    % 如果加载失败，使用传统方法
    warning('三维正则化模型加载失败: %s，尝试使用传统方法', ME.message);
    inverpar.Maprior = readmatrix('mapr.xlsx');
    % 检查是否需要扩展为三维
    if ndims(inverpar.Maprior) == 2
        % 二维数据，扩展为三维
        [m, n] = size(inverpar.Maprior);
        if m == inverpar.numy && n == inverpar.numx
            % 尺寸匹配，扩展为三维
            temp = zeros(inverpar.numz, inverpar.numx, inverpar.numy);
            for k = 1:inverpar.numz
                temp(k, :, :) = inverpar.Maprior';
            end
            inverpar.Maprior = temp;
            disp(['注意: 正则化模型已从二维扩展为三维 (', num2str(inverpar.numz), ' 层)']);
        else
            warning('正则化模型尺寸 (%d×%d) 与期望尺寸 (%d×%d) 不匹配，使用全零模型', ...
                m, n, inverpar.numy, inverpar.numx);
            inverpar.Maprior = zeros(inverpar.numz, inverpar.numx, inverpar.numy);
        end
    elseif ndims(inverpar.Maprior) == 3
        % 三维数据，检查尺寸
        [nz, nx, ny] = size(inverpar.Maprior);
        if nz ~= inverpar.numz || nx ~= inverpar.numx || ny ~= inverpar.numy
            warning('正则化模型尺寸 (%d×%d×%d) 与期望尺寸 (%d×%d×%d) 不匹配，使用全零模型', ...
                nz, nx, ny, inverpar.numz, inverpar.numx, inverpar.numy);
            inverpar.Maprior = zeros(inverpar.numz, inverpar.numx, inverpar.numy);
        end
    else
        % 其他维度，使用全零模型
        warning('正则化模型维度 %d 不正确，使用全零模型', ndims(inverpar.Maprior));
        inverpar.Maprior = zeros(inverpar.numz, inverpar.numx, inverpar.numy);
    end
end
save inverpar.mat inverpar;