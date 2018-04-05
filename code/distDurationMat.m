function map_mat = distDurationMat(metrics_cells) 

[rows,cols] = size(metrics_cells);
map_mat = zeros(rows,cols);

for i = 1:rows
    for j = 1:cols
       metrics_cell = metrics_cells{i,j};
       run =  metrics_cell{5};
       map_mat(i,j) = floor(run/1000)*365+mod(run,1000)/1000*365;
       %map_mat(i,j) = run;
    end
end

end