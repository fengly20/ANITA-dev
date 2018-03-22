function map_mat = dateThresMat(metrics_cells,missing_value)

[rows,cols] = size(metrics_cells);
map_mat = zeros(rows,cols);

for i = 1:rows
    for j = 1:cols
       map_mat(i,j) = metrics_cells{i,j}{12};
    end
end

map_mat(isnan(map_mat)) = missing_value;
end