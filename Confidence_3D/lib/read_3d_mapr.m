function Maprior = read_3d_mapr(filename, numx, numy, numz)
% read_3d_mapr - 从分层Excel文件中读取三维正则化模型
%   Maprior = read_3d_mapr(filename, numx, numy, numz)
% 
% 输入参数:
%   filename - Excel文件名，包含分层存储的正则化模型
%   numx - x方向的网格数量
%   numy - y方向的网格数量
%   numz - z方向的网格数量
% 
% 输出参数:
%   Maprior - 三维正则化模型，维度为numz×numx×numy
% 
% Excel文件结构:
%   - 包含多个工作表，每个工作表对应一个深度层
%   - 工作表命名为'z1', 'z2', ..., 'zn'
%   - 每个工作表的尺寸为numy行×numx列
% 
% 注意:
%   - 当只有一个工作表时，自动扩展为三维模型（所有深度层使用相同值）
%   - 当工作表数量少于numz时，剩余层使用最后一个工作表的值
%   - 当工作表数量多于numz时，只使用前numz个工作表

% 尝试读取Excel文件的所有工作表
[~, sheets] = xlsfinfo(filename);

% 初始化三维正则化模型
Maprior = zeros(numz, numx, numy);

if isempty(sheets)
    % 没有工作表，返回全零模型
    warning('Excel文件 %s 中没有工作表，返回全零正则化模型', filename);
    return;
end

% 检查工作表数量
if length(sheets) == 1
    % 只有一个工作表，自动扩展为三维模型
    try
        % 读取第一个工作表的数据
        data = readmatrix(filename, 'Sheet', sheets{1});
        
        % 检查数据尺寸
        if size(data, 1) ~= numy || size(data, 2) ~= numx
            warning('工作表 %s 的尺寸 (%d×%d) 与期望尺寸 (%d×%d) 不匹配，返回全零正则化模型', ...
                sheets{1}, size(data, 1), size(data, 2), numy, numx);
            return;
        end
        
        % 扩展为三维模型
        for k = 1:numz
            Maprior(k, :, :) = data';  % 转置以匹配numx×numy的顺序
        end
        
        disp(['注意: 只检测到一个工作表，已自动扩展为三维正则化模型 (', num2str(numz), ' 层)']);
        
    catch ME
        warning('读取工作表时出错: %s，返回全零正则化模型', ME.message);
        return;
    end
    
else
    % 多个工作表，按顺序读取
    % 排序工作表（按z1, z2, ..., zn的顺序）
    sorted_sheets = sort_sheets(sheets);
    
    % 读取每个工作表的数据
    for k = 1:min(length(sorted_sheets), numz)
        try
            % 读取当前工作表的数据
            data = readmatrix(filename, 'Sheet', sorted_sheets{k});
            
            % 检查数据尺寸
            if size(data, 1) ~= numy || size(data, 2) ~= numx
                warning('工作表 %s 的尺寸 (%d×%d) 与期望尺寸 (%d×%d) 不匹配，跳过该工作表', ...
                    sorted_sheets{k}, size(data, 1), size(data, 2), numy, numx);
                continue;
            end
            
            % 将数据存入三维模型
            Maprior(k, :, :) = data';  % 转置以匹配numx×numy的顺序
            
        catch ME
            warning('读取工作表 %s 时出错: %s，跳过该工作表', sorted_sheets{k}, ME.message);
            continue;
        end
    end
    
    % 检查是否有足够的工作表
    if length(sorted_sheets) < numz
        warning('工作表数量 (%d) 少于期望的深度层数 (%d)，剩余层使用最后一个工作表的值', ...
            length(sorted_sheets), numz);
        
        % 用最后一个有效工作表的值填充剩余层
        if ~isempty(sorted_sheets)
            try
                last_data = readmatrix(filename, 'Sheet', sorted_sheets{end});
                if size(last_data, 1) == numy && size(last_data, 2) == numx
                    for k = length(sorted_sheets)+1:numz
                        Maprior(k, :, :) = last_data';
                    end
                end
            catch
                % 忽略错误
            end
        end
    elseif length(sorted_sheets) > numz
        warning('工作表数量 (%d) 多于期望的深度层数 (%d)，只使用前 %d 个工作表', ...
            length(sorted_sheets), numz, numz);
    end
end

end

function sorted_sheets = sort_sheets(sheets)
% sort_sheets - 按z1, z2, ..., zn的顺序排序工作表

% 提取工作表名称中的数字
z_numbers = zeros(length(sheets), 1);
valid_sheets = false(length(sheets), 1);

for i = 1:length(sheets)
    sheet_name = lower(sheets{i});
    if startsWith(sheet_name, 'z')
        % 提取数字部分
        num_str = sheet_name(2:end);
        if ~isempty(num_str) && all(isstrprop(num_str, 'digit'))
            z_numbers(i) = str2double(num_str);
            valid_sheets(i) = true;
        end
    end
end

% 只保留有效的工作表（以z开头且后面是数字）
valid_indices = find(valid_sheets);
sorted_sheets = sheets(valid_indices);
z_numbers = z_numbers(valid_indices);

% 按数字排序
if ~isempty(z_numbers)
    [~, sort_order] = sort(z_numbers);
    sorted_sheets = sorted_sheets(sort_order);
end

% 如果没有有效的工作表，返回所有工作表
if isempty(sorted_sheets)
    sorted_sheets = sheets;
end

end