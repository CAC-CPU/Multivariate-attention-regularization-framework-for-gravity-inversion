s1_mat_get


%% ========== 设置 3D 模型的物理尺寸和网格参数 ==========
% 设置模型，就和反演空间一样
x_length = 4000;   % X 方向长度（米）
y_width = 4000;    % Y 方向宽度（米）  
z_height = 1000;   % Z 方向高度（米）
nx = 20;   % X 方向网格点数
ny = 20;   % Y 方向网格点数
nz = 20;   % Z 方向网格点数

% 定义平面方程: Ax + By + Cz = D
A = 0;  % x 系数
B = 1;  % y 系数
C = 0;  % z 系数
D = 2100; % 常数项

% 选择绘图方式: 'pcolor', 'contourf' 或 'imagesc'
plot_type = 'imagesc';


%% ========== 读取数据并处理 ==========
data = importdata('inversion_result.txt');
data = data(:,1);

dx = x_length / (nx - 1);   % X 方向间距
dy = y_width / (ny - 1);    % Y 方向间距
dz = z_height / (nz - 1);   % Z 方向间距

[x, y, z] = meshgrid(0:dx:x_length, 0:dy:y_width, 0:dz:z_height);

new_data = reshape(data, [nx, ny, nz]);
new_data = permute(new_data, [2, 3, 1]);  % 调整维度顺序

num_points = 100;

if B ~= 0
    [x_plane, z_plane] = meshgrid(linspace(min(x(:)), max(x(:)), num_points), ...
                                  linspace(min(z(:)), max(z(:)), num_points));
    y_plane = (D - A * x_plane - C * z_plane) / B;
    valid_mask = (y_plane >= min(y(:))) & (y_plane <= max(y(:)));
elseif A ~= 0
    [y_plane, z_plane] = meshgrid(linspace(min(y(:)), max(y(:)), num_points), ...
                                  linspace(min(z(:)), max(z(:)), num_points));
    x_plane = (D - C * z_plane) / A;
    valid_mask = (x_plane >= min(x(:))) & (x_plane <= max(x(:)));
elseif C ~= 0
    [x_plane, y_plane] = meshgrid(linspace(min(x(:)), max(x(:)), num_points), ...
                                  linspace(min(y(:)), max(y(:)), num_points));
    z_plane = ones(size(x_plane)) * (D / C);
    valid_mask = true(size(x_plane));
else
    error('无效的平面参数: A, B, C 不能同时为 0');
end

x_plane = x_plane(valid_mask);
y_plane = y_plane(valid_mask);
z_plane = z_plane(valid_mask);

v_plane = interp3(x, y, z, new_data, x_plane, y_plane, z_plane, 'linear');

if B ~= 0
    x_range = max(x_plane) - min(x_plane);
    z_range = max(z_plane) - min(z_plane);
    if x_range > 0 && z_range > 0
        aspect_ratio = x_range / z_range;
        screen_ratio = 0.5; 
        if aspect_ratio >= 1
            width = screen_ratio;
            height = screen_ratio / aspect_ratio;
        else
            height = screen_ratio;
            width = screen_ratio * aspect_ratio;
        end
        left = (1 - width) / 2;
        bottom = (1 - height) / 2;
        figure('Units', 'normalized', 'Position', [left bottom width height]);
    else
        figure;
    end
    
    F = scatteredInterpolant(x_plane(:), z_plane(:), v_plane(:), 'linear', 'none');
    [xq, zq] = meshgrid(linspace(min(x_plane), max(x_plane), num_points), ...
                         linspace(min(z_plane), max(z_plane), num_points));
    vq = F(xq, zq);
    if strcmp(plot_type, 'pcolor')
        pcolor(xq, zq, vq);
    elseif strcmp(plot_type, 'contourf')
        contourf(xq, zq, vq, 10, 'LineStyle', 'none');
    elseif strcmp(plot_type, 'imagesc')
        imagesc(xq(1,:), zq(:,1), vq);
    else
        error('不支持的绘图类型，请选择: pcolor, contourf, imagesc');
    end
    xlabel('X (m)', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Z (m)', 'FontSize', 10, 'FontWeight', 'bold');
    set(gca, 'YDir', 'reverse');
    axis equal tight;
    
elseif A ~= 0
    y_range = max(y_plane) - min(y_plane);
    z_range = max(z_plane) - min(z_plane);
    if y_range > 0 && z_range > 0
        aspect_ratio = y_range / z_range;
        screen_ratio = 0.5; 
        if aspect_ratio >= 1
            width = screen_ratio;
            height = screen_ratio / aspect_ratio;
        else
            height = screen_ratio;
            width = screen_ratio * aspect_ratio;
        end
        left = (1 - width) / 2;
        bottom = (1 - height) / 2;
        figure('Units', 'normalized', 'Position', [left bottom width height]);
    else
        figure;
    end
    
    F = scatteredInterpolant(y_plane(:), z_plane(:), v_plane(:), 'linear', 'none');
    [yq, zq] = meshgrid(linspace(min(y_plane), max(y_plane), num_points), ...
                         linspace(min(z_plane), max(z_plane), num_points));
    vq = F(yq, zq);
    if strcmp(plot_type, 'pcolor')
        pcolor(yq, zq, vq);
    elseif strcmp(plot_type, 'contourf')
        contourf(yq, zq, vq, 10, 'LineStyle', 'none');
    elseif strcmp(plot_type, 'imagesc')
        imagesc(yq(1,:), zq(:,1), vq);
    else
        error('不支持的绘图类型，请选择: pcolor, contourf, imagesc');
    end
    xlabel('Y (m)', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Z (m)', 'FontSize', 10, 'FontWeight', 'bold');
    set(gca, 'YDir', 'reverse');
    axis equal tight;
    
else
    x_range = max(x_plane) - min(x_plane);
    y_range = max(y_plane) - min(y_plane);
    if x_range > 0 && y_range > 0
        aspect_ratio = x_range / y_range;
        screen_ratio = 0.5; 
        if aspect_ratio >= 1
            width = screen_ratio;
            height = screen_ratio / aspect_ratio;
        else
            height = screen_ratio;
            width = screen_ratio * aspect_ratio;
        end
        left = (1 - width) / 2;
        bottom = (1 - height) / 2;
        figure('Units', 'normalized', 'Position', [left bottom width height]);
    else
        figure;
    end
    
    F = scatteredInterpolant(x_plane(:), y_plane(:), v_plane(:), 'linear', 'none');
    [xq, yq] = meshgrid(linspace(min(x_plane), max(x_plane), num_points), ...
                         linspace(min(y_plane), max(y_plane), num_points));
    vq = F(xq, yq);
    if strcmp(plot_type, 'pcolor')
        pcolor(xq, yq, vq);
    elseif strcmp(plot_type, 'contourf')
        contourf(xq, yq, vq, 10, 'LineStyle', 'none');
    elseif strcmp(plot_type, 'imagesc')
        imagesc(xq(1,:), yq(:,1), vq);
        set(gca, 'YDir', 'normal');
    else
        error('不支持的绘图类型，请选择: pcolor, contourf, imagesc');
    end
    xlabel('X (m)', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Y (m)', 'FontSize', 10, 'FontWeight', 'bold');
    axis equal tight;
end

shading interp;
colormap(jet);
colorbar;
caxis([min(data) max(data)]);

title('2D slice', 'FontSize', 10, 'FontWeight', 'bold');
set(gca, 'FontSize', 10, 'FontWeight', 'bold');

c = colorbar;
set(c, 'FontSize', 10, 'FontWeight', 'bold');
set(get(c, 'label'), 'string', 'ρ g/cm^3');

drawnow;
clear