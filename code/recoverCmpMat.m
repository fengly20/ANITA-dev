function map_mat = recoverCmpMat(duration,metrics_cells,progress_bar)

[rows,cols] = size(metrics_cells);
map_mat = zeros(rows,cols);

delta = duration*1000;

for i = 1:rows
    parfor j = 1:cols
        metrics_cell = metrics_cells{i,j};
        coeff_before = metrics_cell{12};
        start_date = metrics_cell{3};
        end_date = start_date+delta;
        map_mat(i,j) = valueChange(start_date,end_date,'diff',metrics_cell)/coeff_before;
    end
    if nargin>2 
        if mod(i,10)==0
            i
        end
    end
end
end