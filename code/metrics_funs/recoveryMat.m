function map_mat = recoveryMat(duration,option,metrics_cells, progress_bar)

[rows,cols] = size(metrics_cells);
map_mat = zeros(rows,cols);

delta = duration*1000;

for i = 1:rows
    for j = 1:cols
        metrics_cell = metrics_cells{i,j};
        start_date = metrics_cell{4};
        end_date = start_date+delta;
        map_mat(i,j) = valueChange(start_date,end_date,option,metrics_cell);
    end
    
    if nargin>3 
        if mod(i,10)==0
            i
        end
    end
end
end