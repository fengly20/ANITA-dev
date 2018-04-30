function [head_mat,tail_mat] = headTailMat(metrics_cells)
[rows,cols] = size(metrics_cells);
head_mat = zeros(rows,cols);
tail_mat = zeros(rows,cols);

for i = 1:rows
    for j = 1:cols
        metrics_cell = metrics_cells{i,j};
        if size(metrics_cell{10},2)>1
            head_mat(i,j) =  metrics_cell{10}(1,2);
            tail_mat(i,j) =  metrics_cell{10}(end,2);
        else
            head_mat(i,j) =  NaN;
            tail_mat(i,j) =  NaN;
        end
    end
    
end
end