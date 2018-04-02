function map_mat = recoveryMat(metrics_cells,yr,missing_value,low_high) 

if nargin<4
    low_high = [2 98];
end

[rows,cols] = size(metrics_cells);
map_mat_raw = zeros(rows,cols);

for i = 1:rows
    for j = 1:cols
        if yr == 2
            cell_idx = 10;
        elseif yr == 4
            cell_idx = 11;
        end
       map_mat_raw(i,j) = metrics_cells{i,j}{cell_idx};
    end
end

map_mat_raw(isnan(map_mat_raw)) = missing_value;

%contrast stretch the ouput image
map_mat = vi_stretch(map_mat_raw,low_high);
end