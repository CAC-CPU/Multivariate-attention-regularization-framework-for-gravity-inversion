function ycolor=Jetconvert(m)

colmap = jet(64);
ncolors = size(colmap, 1);

data_min = min(m);
data_max = max(m);
data_range = data_max - data_min;

size_m = length(m);
ycolor = zeros(size_m, 3);

if data_range > 0
    for in = 1:size_m
        normalized_val = (m(in) - data_min) / data_range;
        
        t = 63 * normalized_val;
        
        if t == 0
            t = 0.1;
        end
        
        t = min(t, 63);
        
        color_idx = floor(t) + 1;
        
        if color_idx < 1
            color_idx = 1;
        elseif color_idx > ncolors
            color_idx = ncolors;
        end
        
        ycolor(in, :) = colmap(color_idx, :);
    end
else
    mid_color = colmap(round(ncolors/2), :);
    for in = 1:size_m
        ycolor(in, :) = mid_color;
    end
end

