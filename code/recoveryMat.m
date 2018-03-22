function map_mat = recoveryMat(metrics_cells,yr,missing_value) 

[rows,cols] = size(metrics_cells);
map_mat = zeros(rows,cols);

for i = 1:rows
    for j = 1:cols
        if yr == 2
            cell_idx = 10;
        elseif yr == 4
            cell_idx = 11;
        end
       map_mat(i,j) = metrics_cells{i,j}{cell_idx};
    end
end

map_mat(isnan(map_mat)) = missing_value;
end