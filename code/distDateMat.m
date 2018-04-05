function map_mat = distDateMat(metrics_cells,option) 

% option -- see function calDistDate

[rows,cols] = size(metrics_cells);
map_mat = zeros(rows,cols);

for i = 1:rows
    for j = 1:cols
       metrics_cell = metrics_cells{i,j};
       map_mat(i,j) = calDistDate(metrics_cell,option);
    end
end

end