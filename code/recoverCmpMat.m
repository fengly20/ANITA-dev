function map_mat = recoverCmpMat(duration,metrics_cells)

[rows,cols] = size(metrics_cells);
map_mat = zeros(rows,cols);

delta = duration*1000;

for i = 1:rows
    for j = 1:cols
        if mod(i,10)==0
            i
        end
        
        metrics_cell = metrics_cells{i,j};
        coeff_before = metrics_cell{12};
        start_date = metrics_cell{3};
        end_date = start_date+delta;
        map_mat(i,j) = valueChange(start_date,end_date,'diff',metrics_cell)/coeff_before;
    end
end
end