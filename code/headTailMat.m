function [head_mat,tail_mat] = headTailMat(metrics_cells)
[rows,cols] = size(metrics_cells);
head_mat = zeros(rows,cols);
tail_mat = zeros(rows,cols);

for i = 1:rows
    for j = 1:cols
        metrics_cell = metrics_cells{i,j};
        head_mat(i,j) =  metrics_cell{10}(1,2);
        tail_mat(i,j) =  metrics_cell{10}(end,2);
    end
end
end