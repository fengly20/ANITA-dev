function map_mat = valueChangeMat(start_date,end_time,option,metrics_cells)

% duration can be 9999, or distributed date

[rows,cols] = size(metrics_cells);
map_mat = zeros(rows,cols);

for i = 1:rows
    for j = 1:cols
        
        if start_date == -9999
            metrics_cell = metrics_cells{i,j};
            first_date = metrics_cell{10}(1,1);
            start_date = first_date;
        end
        
        if end_time == 9999
            last_date = metrics_cell{10}(end,1);
            end_date = last_date; 
        end
            
        map_mat(i,j) = valueChange(start_date,end_date,option,metrics_cell);

    end
end
end