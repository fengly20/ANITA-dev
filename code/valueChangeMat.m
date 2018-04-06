function map_mat = valueChangeMat(start_date,end_date,option,metrics_cells,progress_bar)

[rows,cols] = size(metrics_cells);
map_mat = zeros(rows,cols);

for i = 1:rows
    for j = 1:cols
        metrics_cell = metrics_cells{i,j};
        
        if start_date == -9999            
            first_date = metrics_cell{10}(1,1);
            start_date_run = first_date;
        else 
            start_date_run = start_date;
        end
        
        if end_date == 9999
            last_date = metrics_cell{10}(end,1);
            end_date_run = last_date; 
        else 
            end_date_run = end_date;
        end
            
        map_mat(i,j) = valueChange(start_date_run,end_date_run,option,metrics_cell);
        %disp([i j map_mat(i,j)])
    end
    
    if nargin>4
        if mod(i,10)==0
            i
        end
    end
end
end