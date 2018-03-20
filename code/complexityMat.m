function map_mat = complexityMat(results_cells, missing_value)

[rows,cols] = size(results_cells);
map_mat = zeros(rows,cols);

for i = 1:rows
    for j = 1:cols
       map_mat(i,j) = results_cells{i,j}{1};
    end
end

map_mat(isnan(map_mat)) = missing_value;
end
