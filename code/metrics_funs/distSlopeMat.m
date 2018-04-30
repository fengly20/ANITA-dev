function map_mat = distSlopeMat(metrics_cells) 

[rows,cols] = size(metrics_cells);
map_mat = zeros(rows,cols);

for i = 1:rows
    for j = 1:cols
       metrics_cell = metrics_cells{i,j};
       map_mat(i,j) = metrics_cell{6};
    end
end

end